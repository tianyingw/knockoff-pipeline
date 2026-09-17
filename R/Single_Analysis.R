# run_single_block -----------------------------------------------------------
# New params vs original: save_knockoff, load_knockoff, knockoff_file,
#   knockoff_sample_ids, stage1_only.
# When stage1_only=TRUE: saves knockoff then returns NULL immediately.
# ----------------------------------------------------------------------------
run_single_block <- function(
  blocks,
  kk,
  geno.file,
  obj_nullmodel,
  window_length,
  plink_prefix,
  impute.method,
  M,
  Gsub.id,
  bim_metadata,
  genome_build,
  reference_id        = NULL,
  plink_keep_file     = NULL,
  knockoff_seed       = NULL,
  save_knockoff       = FALSE,
  load_knockoff       = FALSE,
  knockoff_file       = NULL,
  knockoff_sample_ids = NULL,
  stage1_only         = FALSE
) {
  chr   <- blocks[kk, chr]
  start <- blocks[kk, start]
  stop  <- blocks[kk, stop]

  # Some reference LD blocks contain no variants in the input dataset.  PLINK
  # reports that situation as a failed/no-output export, but it is an empty
  # analysis unit rather than a pipeline error.  Detect it from the chromosome
  # BIM metadata before invoking PLINK so genuine export failures remain fatal.
  block_bim <- bim_metadata[
    bim_metadata$pos >= start & bim_metadata$pos <= stop, , drop = FALSE
  ]
  if (nrow(block_bim) == 0L) return(NULL)

  tmpdir       <- tempdir()
  block_prefix <- file.path(tmpdir, sprintf("temp_chr%d_block%d", chr, kk))
  keep_arg <- if (is.null(plink_keep_file)) "" else
    paste("--keep", shQuote(plink_keep_file))
  status <- .run_plink_additive_export(
    plink_prefix = plink_prefix, geno_file = geno.file, chr = chr,
    start = start, stop = stop, keep_arg = keep_arg,
    out_prefix = block_prefix
  )
  if (!identical(status, 0L))
    stop("PLINK failed while exporting chr", chr, ":", start, "-", stop, ".")

  raw_file <- paste0(block_prefix, ".raw")
  if (!file.exists(raw_file)) return(NULL)

  raw <- data.table::fread(
    raw_file, data.table = FALSE, check.names = FALSE,
    keepLeadingZeros = TRUE
  )
  unlink(paste0(block_prefix, c(".raw", ".log", ".nosex")), force = TRUE)
  if (ncol(raw) <= 6) return(NULL)

  prepared <- .prepare_raw_genotypes(
    raw = raw, target_ids = Gsub.id, bim_metadata = block_bim
  )
  df <- prepared$geno
  variant_metadata <- prepared$variant_metadata
  rm(raw)
  cat(sprintf("chr: %s, start: %d, end: %d, snp count: %d\n", chr, start, stop, ncol(df)))

  results <- Single_Window_Analysis(
      nullobj             = obj_nullmodel,
      geno                = df,
      chr                 = chr,
      window_length       = window_length,
      M                   = M,
      impute.method       = impute.method,
      Gsub.id             = Gsub.id,
      variant_metadata    = variant_metadata,
      genome_build        = genome_build,
      reference_id        = reference_id,
      knockoff_seed       = knockoff_seed,
      save_knockoff       = save_knockoff,
      load_knockoff       = load_knockoff,
      knockoff_file       = knockoff_file,
      knockoff_sample_ids = knockoff_sample_ids,
      stage1_only         = stage1_only
    )
  rm(df); gc()
  return(results)
}


# Single_Window_Analysis -----------------------------------------------------
# New params: save_knockoff, load_knockoff, knockoff_file,
#   knockoff_sample_ids, stage1_only.
#
# Knockoff save format (RDS manifest + file-backed matrices):
#   $descriptors — serializable descriptors for M file-backed matrices [n × p]
#   $sample_ids   — character IIDs in row order (length n)
#   $snp_pos    — numeric SNP positions (length p); used for column validation
#
# Load validation:
#   1. The complete ordered variant/build/reference/construction context must
#      match; an obsolete or incompatible manifest fails closed.
#   2. The saved and current character-IID sets must match exactly; row order
#      may differ and is aligned explicitly.
# ----------------------------------------------------------------------------
Single_Window_Analysis <- function(
  nullobj,
  geno,
  chr,
  window_length       = NULL,
  M                   = 5,
  thres.single        = 0.01,
  thres.ultrarare     = 25,
  thres.missing       = 0.10,
  midout.dir          = NULL,
  jobtitle            = NULL,
  impute.method       = "fixed",
  Gsub.id             = NULL,
  variant_metadata,
  genome_build,
  reference_id        = NULL,
  knockoff_seed       = NULL,
  bigmemory           = TRUE,
  leveraging          = TRUE,
  LD.filter           = NULL,
  save_knockoff       = FALSE,
  load_knockoff       = FALSE,
  knockoff_file       = NULL,
  knockoff_sample_ids = NULL,
  stage1_only         = FALSE
) {
  imputation_seed <- .derive_unit_seed(knockoff_seed, "imputation")
  preprocess <- .with_local_seed(
    imputation_seed,
    function() Preprocess(
      geno = geno, chr = chr, window = window_length,
      impute.method = impute.method, variant_metadata = variant_metadata
    )
  )
  if (is.null(preprocess)) return(NULL)

  G          <- preprocess$G
  pos        <- preprocess$pos
  variant_metadata <- preprocess$variant_metadata
  window.bed <- preprocess$window.bed
  current_context <- .make_knockoff_context(
    test_type = "Single_Window", M = M, genome_build = genome_build,
    variant_metadata = variant_metadata, reference_id = reference_id,
    construction_id = paste0(
      "KnockoffScreen-SCIP-v1;impute=", impute.method,
      ";imputation_seed=", if (is.null(imputation_seed)) "NULL" else imputation_seed,
      ";corr_max=0.75;maxBP=100000;thres_ultrarare=25;R2=1;method=shrinkage"
    ),
    random_seed = knockoff_seed
  )

  # ---- Knockoff: load or generate -----------------------------------------
  G_k          <- NULL
  need_generate <- TRUE

  if (isTRUE(load_knockoff)) {
    if (is.null(knockoff_file) || !file.exists(knockoff_file))
      stop("Required saved knockoff file not found: ", knockoff_file)
    ko_obj <- readRDS(knockoff_file)
    .assert_knockoff_context(ko_obj$context, current_context, knockoff_file)

    if (identical(ko_obj$storage, "bigmemory_filebacked")) {
      if (length(ko_obj$descriptor_files) != M)
        stop("Saved knockoff manifest has the wrong number of matrix descriptors: ", knockoff_file)
      ko_obj$G_k <- lapply(ko_obj$descriptor_files, function(desc) {
        desc_path <- file.path(dirname(knockoff_file), desc)
        if (!file.exists(desc_path))
          stop("Saved knockoff descriptor not found: ", desc_path)
        bigmemory::attach.big.matrix(
          basename(desc_path), path = dirname(desc_path)
        )
      })
    } else if (identical(ko_obj$storage, "r_matrix")) {
      if (length(ko_obj$G_k) != M)
        stop("Saved knockoff manifest has the wrong number of matrices: ", knockoff_file)
    } else {
      stop("Saved knockoff storage metadata are missing or unsupported: ", knockoff_file,
           ". Regenerate knockoffs with the current package version.")
    }

    target_ids <- if (!is.null(Gsub.id)) Gsub.id else knockoff_sample_ids
    saved_ids <- ko_obj$sample_ids
    if (.need_regenerate_samples(target_ids, saved_ids)) {
      stop(
        "Saved knockoff cannot be reused because its sample-ID set differs ",
        "from the current run: ", knockoff_file, ". Regenerate knockoffs for ",
        "the exact analysis sample set."
      )
    }
    ko_obj <- .align_knockoff_samples(ko_obj, target_ids)
    G_k <- ko_obj$G_k
    need_generate <- FALSE
  }

  if (need_generate) {
    backing_path <- if (isTRUE(save_knockoff) && isTRUE(bigmemory) &&
                       !is.null(knockoff_file)) dirname(knockoff_file) else NULL
    backing_prefix <- if (is.null(backing_path)) NULL else basename(tempfile(
      pattern = paste0(
        tools::file_path_sans_ext(basename(knockoff_file)), "_matrix-"
      ),
      tmpdir = backing_path
    ))
    descriptor_files <- if (is.null(backing_prefix)) NULL else
      paste0(backing_prefix, "_", seq_len(M), ".desc")
    new_backing_files <- if (is.null(descriptor_files)) character(0) else
      file.path(
        backing_path,
        c(descriptor_files, sub("\\.desc$", ".bin", descriptor_files))
      )
    checkpoint_committed <- FALSE
    on.exit({
      if (!checkpoint_committed && length(new_backing_files) > 0L) {
        # The manifest may already have been atomically installed immediately
        # before an interrupt.  Never remove files referenced by the currently
        # installed manifest, even if the in-memory flag was not yet updated.
        unlink(
          .uncommitted_backing_cleanup(new_backing_files, knockoff_file),
          force = TRUE
        )
      }
    }, add = TRUE)
    G_k <- .with_local_seed(
      .derive_unit_seed(knockoff_seed, "knockoff"),
      function() create.KS(
        X = G, pos = pos, M = M, bigmemory = bigmemory,
        backing_path = backing_path, backing_prefix = backing_prefix
      )
    )

    if (isTRUE(save_knockoff) && !is.null(knockoff_file)) {
      dir.create(dirname(knockoff_file), recursive = TRUE, showWarnings = FALSE)
      storage <- if (isTRUE(bigmemory)) "bigmemory_filebacked" else "r_matrix"
      if (!isTRUE(bigmemory)) descriptor_files <- NULL
      previous_manifest <- if (file.exists(knockoff_file)) {
        tryCatch(readRDS(knockoff_file), error = function(e) NULL)
      } else NULL
      manifest <- list(
          storage    = storage,
          G_k        = if (identical(storage, "r_matrix")) G_k else NULL,
          descriptor_files = descriptor_files,
          sample_ids = if (!is.null(Gsub.id)) Gsub.id else knockoff_sample_ids,
          snp_pos    = pos,
          context    = current_context
      )
      .atomic_save_rds(manifest, knockoff_file)
      checkpoint_committed <- TRUE

      # Only after the new manifest is installed do we retire backing files
      # referenced by the previous generation.
      old_files <- setdiff(
        .manifest_backing_files(previous_manifest, knockoff_file),
        new_backing_files
      )
      if (length(old_files) > 0L) unlink(old_files, force = TRUE)
    }
  }

  # ---- Stage 1: return after saving knockoff (no association test) --------
  if (isTRUE(stage1_only)) return(invisible(NULL))

  # ---- Association test ---------------------------------------------------
  fit <- KS.chr(
    result.prelim = nullobj,
    input.X       = G,
    window.bed    = window.bed,
    input.G_k     = G_k,
    M             = M,
    thres.single  = thres.single,
    Gsub.id       = Gsub.id
  )
  rm(G_k, G, pos, window.bed); gc()
  return(fit)
}


# Preprocess -----------------------------------------------------------------
Preprocess <- function(geno, chr, window = NULL, thres.maf = 0,
                       thres.missing = 0.1, impute.method = "fixed",
                       variant_metadata = NULL) {

  G <- as.matrix(geno)
  if (is.null(variant_metadata) || nrow(variant_metadata) != ncol(G))
    stop("variant_metadata must contain one row per genotype column.")
  G <- 2 - G

  if (length(G) == 0 || ncol(G) == 0) {
    warning("Number of variants in the specified range is 0", call. = FALSE)
    return(NULL)   # FIX: was `next`
  }
  if (ncol(G) == 1) {
    warning("Number of variants in the specified range is 1", call. = FALSE)
    return(NULL)   # FIX: was `next`
  }

  variant_key <- paste(
    variant_metadata$chr, variant_metadata$variant_id, variant_metadata$pos,
    variant_metadata$a1, variant_metadata$a2, sep = ":"
  )
  unique_variant <- match(unique(variant_key), variant_key)
  G <- G[, unique_variant, drop = FALSE]
  variant_metadata <- variant_metadata[unique_variant, , drop = FALSE]

  # Missing imputation
  G[G < 0 | G > 2] <- NA
  G <- Impute(G, impute.method)

  # Filter constant variants
  s <- apply(G, 2, sd)
  keep_variable <- !is.na(s) & s != 0
  G <- G[, keep_variable, drop = FALSE]
  variant_metadata <- variant_metadata[keep_variable, , drop = FALSE]
  if (ncol(G) < 2) return(NULL)

  # Reorder by position
  pos <- as.numeric(variant_metadata$pos)
  pos_order <- order(pos, variant_metadata$variant_id)
  G <- G[, pos_order, drop = FALSE]
  variant_metadata <- variant_metadata[pos_order, , drop = FALSE]
  pos <- pos[pos_order]
  MAF <- colMeans(G) / 2
  G   <- as.matrix(G)
  flip_to_minor <- MAF > 0.5 & !is.na(MAF)
  G[, flip_to_minor] <- 2 - G[, flip_to_minor, drop = FALSE]
  variant_metadata$coded_allele <- ifelse(
    flip_to_minor, variant_metadata$counted_allele,
    ifelse(variant_metadata$counted_allele == variant_metadata$a1,
           variant_metadata$a2, variant_metadata$a1)
  )
  MAF <- colMeans(G) / 2
  MAC <- colSums(G)
  G   <- Matrix::Matrix(G, sparse = TRUE)

  colnames(G) <- pos
  start <- min(pos); end <- max(pos)

  # Clustering to remove highly correlated SNPs
  cor.X      <- sparse.cor(Matrix::Matrix(G))$cor
  Sigma.dist <- as.dist(1 - abs(cor.X))
  fit_clust  <- hclust(Sigma.dist, method = "single")
  clusters   <- cutree(fit_clust, h = 1 - 0.75)

  cluster.idx <- match(unique(clusters), clusters)
  G   <- G[, cluster.idx, drop = FALSE]
  MAF <- MAF[cluster.idx]; MAC <- MAC[cluster.idx]; pos <- pos[cluster.idx]
  variant_metadata <- variant_metadata[cluster.idx, , drop = FALSE]

  unique.idx <- match(unique(pos), pos)
  G   <- G[, unique.idx, drop = FALSE]
  MAF <- MAF[unique.idx]; MAC <- MAC[unique.idx]; pos <- pos[unique.idx]
  variant_metadata <- variant_metadata[unique.idx, , drop = FALSE]

  if (ncol(G) < 2) return(NULL)

  # Window bed
  if (length(window) != 0) {
    window.bed <- c()
    for (size in window) {
      pos.tag    <- seq(start, end, by = size * 0.5)
      window.bed <- rbind(window.bed, cbind(chr, pos.tag, pos.tag + size))
    }
    window.bed <- window.bed[order(as.numeric(window.bed[, 2])), ]
    return(list(G = G, chr = chr, pos = pos, window.bed = window.bed,
                variant_metadata = variant_metadata))
  } else {
    return(list(G = G, chr = chr, pos = pos, window.bed = NULL,
                variant_metadata = variant_metadata))
  }
}
