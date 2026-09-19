# run_batch_gene -------------------------------------------------------------
# New params vs original: save_knockoff, load_knockoff, knockoff_dir,
#   knockoff_sample_ids, stage1_only, read_mid_exist, export_switch,
#   plink_threads, temp_dir.
#
# knockoff_dir is the chr-level subdirectory (e.g. <knockoff_root>/chr1/).
# One RDS file per gene is written there for the gene_buffer knockoff only.
# Enhancer knockoffs are always generated fresh (fast + many; not saved).
# ----------------------------------------------------------------------------
run_batch_gene <- function(
  genes,
  b,
  batch_index,
  geno.file,
  obj_nullmodel,
  window_length,
  plink_prefix,
  M,
  genome_build,
  Gsub.id,
  bim_metadata,
  plink_keep_file     = NULL,
  reference_id        = NULL,
  seed                = NULL,
  use_glmm            = FALSE,
  abc_df,
  gh_df,
  sparseSigma         = NULL,
  ratio               = NULL,
  user_cores          = 1,
  save_knockoff       = FALSE,
  load_knockoff       = FALSE,
  knockoff_dir        = NULL,       # chr-level subdir, e.g. <root>/chr1
  knockoff_sample_ids = NULL,
  stage1_only         = FALSE,
  read_mid_exist      = TRUE,
  export_switch       = NULL,
  plink_threads       = NULL,
  temp_dir            = NULL
) {
  kk_vec <- batch_index[[b]]
  tmpdir <- if (is.null(temp_dir)) tempdir() else temp_dir
  if (!dir.exists(tmpdir) &&
      !dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE))
    stop("Unable to create temporary directory: ", tmpdir)
  chr    <- as.numeric(gsub("chr", "", genes[kk_vec[1], chr]))

  gene_buffer_extension <- 5000 + 50000
  start_all    <- min(genes[kk_vec, start]) - gene_buffer_extension
  end_all      <- max(genes[kk_vec, end])   + gene_buffer_extension
  batch_bim <- bim_metadata[
    bim_metadata$pos >= start_all & bim_metadata$pos <= end_all, , drop = FALSE
  ]
  # A gene batch can legitimately have no variants in the input dataset.
  # Skip that empty analysis unit before PLINK turns it into a no-output error.
  if (nrow(batch_bim) == 0L) return(NULL)

  batch_prefix <- tempfile(
    sprintf("KnockoffPipeline_chr%d_batch_%d_%d_",
            chr, min(kk_vec), max(kk_vec)),
    tmpdir = tmpdir
  )

  keep_arg <- if (is.null(plink_keep_file)) "" else
    paste("--keep", shQuote(plink_keep_file))
  status <- .run_plink_additive_export(
    plink_prefix = plink_prefix, geno_file = geno.file, chr = chr,
    start = start_all, stop = end_all, keep_arg = keep_arg,
    out_prefix = batch_prefix, export_switch = export_switch,
    plink_threads = plink_threads
  )
  if (!identical(status, 0L))
    stop("PLINK failed while exporting chr", chr, ":", start_all, "-", end_all, ".")

  raw_file <- paste0(batch_prefix, ".raw")
  if (!file.exists(raw_file)) return(NULL)
  raw <- data.table::fread(
    raw_file, data.table = FALSE, check.names = FALSE,
    keepLeadingZeros = TRUE
  )
  unlink(paste0(batch_prefix, c(".raw", ".log", ".nosex")), force = TRUE)
  if (ncol(raw) <= 6) return(NULL)

  message("  Batch ", b, " / ", length(batch_index),
              " (snp ", start_all, "-", end_all, ")")

  prepared <- .prepare_raw_genotypes(
    raw = raw, target_ids = Gsub.id, bim_metadata = batch_bim
  )
  G_batch <- prepared$geno
  variants_batch <- as.numeric(prepared$variant_metadata$pos)
  variant_metadata_batch <- prepared$variant_metadata
  # Downstream upstream-method functions expect numeric positions in the
  # genotype column names.  These values now come from .bim, never from rsID
  # digits or another filename convention.
  colnames(G_batch) <- as.character(variants_batch)
  rm(raw); gc()

  ## ===== Per-gene function =====
  safe_fun <- function(kk) {
    tryCatch({
      gene_start <- genes[kk, start]
      gene_end   <- genes[kk, end]
      gene_id    <- genes[kk, id]
      knockoff_seed <- .derive_unit_seed(
        seed, "Gene_Centric", chr, as.character(gene_id)
      )

      # Knockoff file path for this gene (gene_buffer knockoff only)
      ko_file <- if (!is.null(knockoff_dir))
        file.path(knockoff_dir,
                  paste0("gene_", gsub("[^a-zA-Z0-9._-]", "_", gene_id), "_ko.rds"))
      else NULL

      load_this_gene <- isTRUE(load_knockoff) ||
        (isTRUE(stage1_only) && isTRUE(read_mid_exist) &&
           !is.null(ko_file) && file.exists(ko_file))
      save_this_gene <- isTRUE(save_knockoff) && !load_this_gene

      # Gene buffer SNPs (±5kb around gene body)
      idx_gene_buffer <- which(variants_batch >= gene_start-5000 & variants_batch <= gene_end+5000)
      idx_gene_surround <- which(variants_batch >= gene_start-gene_buffer_extension & variants_batch <= gene_end+gene_buffer_extension)
      
      if (length(idx_gene_buffer) <= 1) return(NULL)
      G_gene          <- G_batch[, idx_gene_surround, drop = FALSE]
      variant_metadata_gene <- variant_metadata_batch[
        idx_gene_surround, , drop = FALSE
      ]
      gene_buffer.pos <- c(min(variants_batch[idx_gene_buffer]),
                           max(variants_batch[idx_gene_buffer]))

      # Enhancer regions
      abc_enhancers <- abc_df[TargetGene == gene_id, .(start, end)]
      gh_enhancers  <- gh_df[gene == gene_id,
                             .(start = GH_start, end = GH_end)]
      enhancers     <- unique(rbind(abc_enhancers, gh_enhancers))

      G_EnhancerAll_surround        <- NULL
      variants_EnhancerAll_surround <- NULL
      Enhancer.pos                  <- NULL
      p_EnhancerAll_surround        <- NULL
      p_EnhancerAll                 <- NULL
      R <- 0

      if (nrow(enhancers) > 0) {
        for (r in seq_len(nrow(enhancers))) {
          e_start <- enhancers$start[r]; e_end <- enhancers$end[r]
          idx_e_surrond <- which(variants_batch >= e_start-5000 & variants_batch <= e_end+5000)
          idx_e <- which(variants_batch >= e_start & variants_batch <= e_end)

          if (length(idx_e) > 5) {
            G_EnhancerAll_surround        <- cbind(G_EnhancerAll_surround,
                                                   G_batch[, idx_e_surrond, drop = FALSE])
            Enhancer.pos                  <- rbind(Enhancer.pos, c(e_start, e_end))
            variants_EnhancerAll_surround <- c(variants_EnhancerAll_surround,
                                               variants_batch[idx_e_surrond])
            p_EnhancerAll_surround        <- c(p_EnhancerAll_surround, length(idx_e_surrond))
            p_EnhancerAll                 <- c(p_EnhancerAll, length(idx_e))
            R <- R + 1
          }
        }
      }

      # Dispatch to analysis function
      if (!isTRUE(use_glmm)) {
        full_results <- GeneScan3D.KnockoffGeneration(
          G_gene_buffer_surround        = G_gene,
          variants_gene_buffer_surround = variants_batch[idx_gene_surround],
          gene_buffer.pos               = gene_buffer.pos,
          R                             = R,
          G_EnhancerAll_surround        = G_EnhancerAll_surround,
          variants_EnhancerAll_surround = variants_EnhancerAll_surround,
          p_EnhancerAll_surround        = p_EnhancerAll_surround,
          Enhancer.pos                  = Enhancer.pos,
          p.EnhancerAll                 = p_EnhancerAll,
          window.size                   = window_length,
          result.null.model             = obj_nullmodel,
          M                             = M,
          Gsub.id                       = Gsub.id,
          variant_metadata_gene_buffer_surround = variant_metadata_gene,
          genome_build                  = genome_build,
          reference_id                  = reference_id,
          knockoff_seed                 = knockoff_seed,
          save_knockoff                 = save_this_gene,
          load_knockoff                 = load_this_gene,
          knockoff_file                 = ko_file,
          knockoff_sample_ids           = knockoff_sample_ids,
          stage1_only                   = stage1_only
        )
      } else {
        full_results <- GeneScan3D.UKB.GLMM.KnockoffGeneration(
          G_gene_buffer_surround        = G_gene,
          variants_gene_buffer_surround = variants_batch[idx_gene_surround],
          gene_buffer.pos               = gene_buffer.pos,
          R                             = R,
          G_EnhancerAll_surround        = G_EnhancerAll_surround,
          variants_EnhancerAll_surround = variants_EnhancerAll_surround,
          p_EnhancerAll_surround        = p_EnhancerAll_surround,
          Enhancer.pos                  = Enhancer.pos,
          p.EnhancerAll                 = p_EnhancerAll,
          window.size                   = window_length,
          result.null.model             = obj_nullmodel,
          M                             = M,
          Gsub.id                       = Gsub.id,
          sparseSigma                   = sparseSigma,
          ratio                         = ratio,
          variant_metadata_gene_buffer_surround = variant_metadata_gene,
          genome_build                  = genome_build,
          reference_id                  = reference_id,
          knockoff_seed                 = knockoff_seed,
          save_knockoff                 = save_this_gene,
          load_knockoff                 = load_this_gene,
          knockoff_file                 = ko_file,
          knockoff_sample_ids           = knockoff_sample_ids,
          stage1_only                   = stage1_only
        )
      }

      if (isTRUE(stage1_only) || is.null(full_results)) return(invisible(NULL))

      results <- data.frame(
        chr        = chr,
        gene_id    = gene_id,
        gene_start = gene_start,
        gene_end   = gene_end,
        GeneScan3D.Cauchy = full_results$GeneScan3D.Cauchy[1],
        t(full_results$GeneScan3D.Cauchy_knockoff[, 1, drop = FALSE]),
        stringsAsFactors = FALSE
      )
      colnames(results)[6:ncol(results)] <-
        paste0("GeneScan3D.Cauchy_knockoff_", seq_len(M))
      return(results)

    }, error = function(e) {
      if (isTRUE(load_knockoff) || isTRUE(stage1_only)) stop(e)
      message("  !! Gene ", genes[kk, id], " (chr ", chr, ") failed: ",
              conditionMessage(e))
      NULL
    })
  }

  out <- parallel::mclapply(kk_vec, safe_fun, mc.cores = user_cores)
  failed <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("Gene batch ", b, " failed while validating or loading saved knockoffs: ",
         paste(as.character(out[failed]), collapse = "; "))
  }
  out <- Filter(Negate(is.null), out)
  rm(G_batch); gc()

  if (length(out) == 0) return(NULL)
  # FIX (original bug): was rbindlist(as.data.table(result_list)) — wrong nesting
  return(data.table::rbindlist(out, fill = TRUE))
}


# .gene_ko_load_or_gen -------------------------------------------------------
# Load gene_buffer knockoff from RDS, or generate fresh.
# Validates test type, M, build, reference identifier, and the complete ordered
# variant/allele fingerprint. Any mismatch fails closed.
# The saved and current character-IID sets must match exactly; row order may
# differ and is aligned explicitly.
#
# Knockoff RDS format:
#   $G_gene_buffer_knockoff  — array [M × n × p_gene_buffer]
#   $sample_ids              — character, length n (matched_ids from stage1)
#   $snp_pos                 — numeric, length p_gene_buffer
# ---------------------------------------------------------------------------
.gene_ko_load_or_gen <- function(
  load_knockoff,
  save_knockoff,
  knockoff_file,
  matched_ids,
  p_expected,    # number of SNPs in gene buffer after QC (current run)
  gen_fun,       # zero-arg function that returns the knockoff array
  snp_pos,       # current SNP positions for saving
  context
) {
  need_generate <- TRUE

  if (isTRUE(load_knockoff)) {
    if (is.null(knockoff_file) || !file.exists(knockoff_file))
      stop("Required saved knockoff file not found: ", knockoff_file)
    ko_obj <- readRDS(knockoff_file)
    .assert_knockoff_context(ko_obj$context, context, knockoff_file)
    arr <- ko_obj$G_gene_buffer_knockoff
    if (length(dim(arr)) != 3L || dim(arr)[1L] != context$M ||
        dim(arr)[3L] != p_expected)
      stop("Saved gene knockoff array has incompatible dimensions: ", knockoff_file)

    current_ids <- as.character(matched_ids)
    saved_ids <- as.character(ko_obj$sample_ids)
    if (.need_regenerate_samples(current_ids, saved_ids)) {
      stop(
        "Saved knockoff cannot be reused because its sample-ID set differs ",
        "from the current run: ", knockoff_file, ". Regenerate knockoffs for ",
        "the exact analysis sample set."
      )
    }
    row_map <- match(current_ids, saved_ids)
    arr <- arr[, row_map, , drop = FALSE]
    need_generate <- FALSE
  }

  if (need_generate) {
    arr <- gen_fun()
    if (is.null(arr)) return(NULL)
    if (length(dim(arr)) != 3L || dim(arr)[1L] != context$M ||
        dim(arr)[2L] != length(matched_ids) || dim(arr)[3L] != p_expected)
      stop("Generated gene knockoff array has incompatible dimensions.")

    if (isTRUE(save_knockoff) && !is.null(knockoff_file)) {
      dir.create(dirname(knockoff_file), recursive = TRUE, showWarnings = FALSE)
      .atomic_save_rds(
        list(
          G_gene_buffer_knockoff = arr,
          sample_ids             = matched_ids,
          snp_pos                = snp_pos,
          context                = context
        ),
        path = knockoff_file
      )
    }
  }

  arr
}


# ---- Summary helpers (unchanged from original) ----------------------------

extract_position_universal <- function(col_names) {
  positions <- sapply(col_names, function(col) {
    numbers <- regmatches(col, gregexpr("\\d+", col))[[1]]
    if (length(numbers) == 0) return(NA)
    if (length(numbers) == 1) return(as.numeric(numbers))
    candidate <- numbers[nchar(numbers) >= 3 & nchar(numbers) <= 9]
    if (length(candidate) > 0) return(as.numeric(utils::tail(candidate, 1)))
    return(as.numeric(max(numbers)))
  })
  return(as.numeric(positions))
}

preprocess_for_GeneScan3DKnock <- function(p0, p_ko, M) {
  p0   <- as.numeric(p0);  p0[is.na(p0)]     <- 0.5
  p_ko <- as.matrix(p_ko); p_ko[is.na(p_ko)] <- 0.5
  list(p0 = p0, p_ko = p_ko)
}

GeneScan3DKnock_Summary <- function(result, M, fdr = 0.1) {
  result <- as.data.frame(result)
  result <- result[order(result[, 4]), ]
  result <- result[order(result[, 3]), ]
  p0  <- as.numeric(result[, 5])
  pk  <- as.matrix(result[, 6:(5 + M), drop = FALSE])
  pre <- preprocess_for_GeneScan3DKnock(p0, pk, M)
  res <- GeneScan3DKnock(M = M, p0 = pre$p0, p_ko = pre$p_ko,
                         gene_id = result[, 2], fdr = fdr)
  result$W           <- res$W
  result$W_Threshold <- rep(res$W.threshold, nrow(result))
  result$Qvalue      <- res$Qvalue
  result$indicator   <- res$Qvalue <= fdr
  return(result)
}
