#' Run Knockoff Genome-wide Pipeline
#'
#' Main entry function for the KnockoffPipeline package. Supports both
#' single-window and gene-centric association testing with knockoff-based
#' FDR control.
#'
#' @param outdir        Character. Output directory (created recursively if absent).
#' @param test_type     Character. One of \code{"Single_Window"} or
#'   \code{"Gene_Centric"}.
#' @param pheno_file    Character. Path to the phenotype file (tab- or
#'   comma-separated).
#' @param geno_file     Character. PLINK genotype file prefix (no extension).
#' @param phenotype     Character. Column name of the outcome variable.
#' @param pheno_id      Character or \code{NULL}. Column name of the sample ID
#'   in the phenotype file. \code{NULL} assumes rows are already aligned with
#'   the PLINK \code{.fam} file.
#' @param covar_cols    Character vector of continuous covariate column names,
#'   or \code{NULL}.
#' @param cat_covar_cols Character vector of binary/categorical covariate column
#'   names, or \code{NULL}.
#' @param user_cores    Integer. Number of parallel cores for within-chromosome
#'   processing. Default \code{1}.
#' @param sliding_window_length Integer or \code{NULL}. Sliding window size (bp).
#' @param geno_missing_imputation Character. Genotype imputation method.
#'   Default \code{"fixed"}.
#' @param plink_path    Character. Path to the PLINK executable. Default
#'   \code{"plink"} (assumes it is on \code{PATH}).
#' @param M             Integer. Number of knockoff copies. Default \code{5}.
#' @param genome_build  Character. One of \code{"hg19"} or \code{"hg38"}.
#' @param sample_uncorrelated Logical. \code{TRUE} fits a standard (non-mixed)
#'   null model; \code{FALSE} fits a GLMM via SAIGE.
#' @param grm_file      Character or \code{NULL}. Path to the sparse GRM file
#'   (required when \code{sample_uncorrelated = FALSE}).
#' @param grm_id_file   Character or \code{NULL}. Path to the sparse GRM ID
#'   file (required when \code{sample_uncorrelated = FALSE}).
#' @param fdr           Numeric in \code{(0, 1)}. Target FDR level. Default
#'   \code{0.1}.
#' @param chromosomes   Integer vector of autosomes to analyse. Defaults to
#'   \code{1:22}.
#' @param batch_size    Integer. Number of genes per batch in gene-centric mode.
#'   Default \code{20}.
#' @param read_mid_exist Logical. Skip chromosomes whose intermediate result
#'   files already exist. Default \code{TRUE}.
#'
#' @return Invisibly returns \code{TRUE} on success.
#'
#' @import SKAT
#' @import Matrix
#' @import WGScan
#' @import SPAtest
#' @import CompQuadForm
#' @import irlba
#' @import bigmemory
#' @import data.table
#' @import dplyr
#' @import parallel
#' @import qqman
#' @import abind
#' @import SAIGE
#' @export
run_pipeline <- function(
  outdir,
  test_type,
  pheno_file,
  geno_file,
  phenotype,
  pheno_id = NULL,
  covar_cols         = NULL,
  cat_covar_cols     = NULL,
  user_cores = 1L,
  sliding_window_length = NULL,
  geno_missing_imputation = "fixed",
  plink_path = "plink",
  M = 5L,
  genome_build = "hg19",
  sample_uncorrelated = TRUE,
  grm_file = NULL,
  grm_id_file = NULL,
  fdr = 0.1,
  chromosomes = 1:22,
  batch_size = 20L,
  read_mid_exist = TRUE
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  stopifnot(
    "outdir must be a single non-empty string"         = is.character(outdir) && length(outdir) == 1L && nzchar(outdir),
    "pheno_file must be a single non-empty string"     = is.character(pheno_file) && length(pheno_file) == 1L && nzchar(pheno_file),
    "geno_file must be a single non-empty string"      = is.character(geno_file) && length(geno_file) == 1L && nzchar(geno_file),
    "phenotype must be a single non-empty string"      = is.character(phenotype) && length(phenotype) == 1L && nzchar(phenotype),
    "M must be a positive integer"                     = is.numeric(M) && length(M) == 1L && M >= 1L,
    "fdr must be a numeric value in (0, 1)"            = is.numeric(fdr) && length(fdr) == 1L && fdr > 0 && fdr < 1,
    "user_cores must be a positive integer"            = is.numeric(user_cores) && length(user_cores) == 1L && user_cores >= 1L,
    "batch_size must be a positive integer"            = is.numeric(batch_size) && length(batch_size) == 1L && batch_size >= 1L,
    "sample_uncorrelated must be logical"              = is.logical(sample_uncorrelated) && length(sample_uncorrelated) == 1L,
    "read_mid_exist must be logical"                   = is.logical(read_mid_exist) && length(read_mid_exist) == 1L
  )

  if (!test_type %in% c("Single_Window", "Gene_Centric")) {
    stop("'test_type' must be one of \"Single_Window\" or \"Gene_Centric\".")
  }

  if (!genome_build %in% c("hg19", "hg38")) {
    stop("'genome_build' must be one of \"hg19\" or \"hg38\".")
  }

  if (!is.null(covar_cols) && !is.character(covar_cols)) {
    stop("'covar_cols' must be a character vector or NULL.")
  }

  if (!is.null(cat_covar_cols) && !is.character(cat_covar_cols)) {
    stop("'cat_covar_cols' must be a character vector or NULL.")
  }

  if (!is.null(pheno_id) && (!is.character(pheno_id) || length(pheno_id) != 1L)) {
    stop("'pheno_id' must be a single string or NULL.")
  }

  # Validate chromosomes: accept integer-coercible values only
  chr_numeric <- suppressWarnings(as.integer(chromosomes))
  if (any(is.na(chr_numeric))) {
    stop("'chromosomes' contains values that cannot be coerced to integers: ",
         paste(chromosomes[is.na(chr_numeric)], collapse = ", "))
  }
  chr_vector <- intersect(chr_numeric, 1L:22L)
  if (length(chr_vector) == 0L) {
    stop("No valid autosomes (1-22) found in 'chromosomes'.")
  }

  if (!file.exists(pheno_file)) {
    stop("Phenotype file not found: ", pheno_file)
  }

  plink_fam <- paste0(geno_file, ".fam")
  if (!file.exists(plink_fam)) {
    stop("PLINK .fam file not found: ", plink_fam,
         "\nCheck that 'geno_file' is the correct prefix.")
  }

  # ---------------------------------------------------------------------------
  # Setup
  # ---------------------------------------------------------------------------

  if (!dir.exists(outdir)) {
    message("Creating output directory: ", outdir)
    dir.create(outdir, recursive = TRUE)
  }

  mid_dir <- file.path(outdir, "mid")
  if (!dir.exists(mid_dir)) {
    dir.create(mid_dir)
  }

  # Limit BLAS thread contention when using mclapply
  if (user_cores > 1L) {
    Sys.setenv(MKL_NUM_THREADS = 1)
  }

  # ---------------------------------------------------------------------------
  # Load phenotype and detect outcome type
  # ---------------------------------------------------------------------------
  
  message("Reading phenotype file: ", pheno_file)
  pheno <- data.table::fread(pheno_file)

  if (!phenotype %in% colnames(pheno)) {
    stop("Phenotype column \"", phenotype, "\" not found in phenotype file.")
  }

  all_covar_cols <- c(covar_cols, cat_covar_cols)
  missing_covar  <- setdiff(all_covar_cols, colnames(pheno))
  if (length(missing_covar) > 0L) {
    stop("The following covariate column(s) are missing from the phenotype file: ",
         paste(missing_covar, collapse = ", "))
  }

  if (!is.null(pheno_id) && !pheno_id %in% colnames(pheno)) {
    stop("Sample ID column \"", pheno_id, "\" not found in phenotype file.")
  }

  pheno_values <- unique(pheno[[phenotype]])
  is_binary    <- length(pheno_values) == 2L && all(sort(pheno_values) == c(0, 1))
  out_type     <- if (is_binary) "D" else "C"
  message("Outcome type detected: ", if (is_binary) "binary (D)" else "continuous (C)")

  # ---------------------------------------------------------------------------
  # Build covariate matrix (may be NULL when no covariates are specified)
  # ---------------------------------------------------------------------------

  .build_covar_matrix <- function(dt, cols) {
    if (length(cols) == 0L) return(NULL)
    as.matrix(dt[, cols, with = FALSE])
  }

  # ---------------------------------------------------------------------------
  # Align phenotype file to PLINK .fam by sample ID (when pheno_id provided)
  # ---------------------------------------------------------------------------
  # Gsub.id semantics: a character/integer vector of IIDs present in the .fam
  # file that should be used in analysis, in .fam order.  When pheno_id is
  # NULL the caller guarantees pheno rows are already in .fam order.
 
  fam <- data.table::fread(plink_fam, header = FALSE,
                           col.names = c("FID","IID","PAT","MAT","SEX","PHENO"))
  if (!is.null(pheno_id)) {
 
    pheno_iid <- as.numeric(pheno[[pheno_id]])
    fam_iid   <- as.numeric(fam$IID)
 
    # Samples present in both phenotype file and .fam
    shared_iid <- intersect(fam_iid, pheno_iid)
 
    if (length(shared_iid) == 0L) {
      stop("No samples matched between the phenotype file (column \"", pheno_id,
           "\") and the PLINK .fam file.\n",
           "  Example pheno IID : ", paste(head(pheno_iid, 3L), collapse = ", "), "\n",
           "  Example .fam  IID : ", paste(head(fam_iid,   3L), collapse = ", "), "\n",
           "Check that both files use the same sample ID format.")
    }
 
    n_pheno_only <- length(setdiff(pheno_iid, fam_iid))
    n_fam_only   <- length(setdiff(fam_iid, pheno_iid))
    if (n_pheno_only > 0L)
      message("  ", n_pheno_only, " sample(s) in phenotype file not found in .fam — excluded.")
    if (n_fam_only > 0L)
      message("  ", n_fam_only,   " sample(s) in .fam not found in phenotype file — excluded.")
    message("  ", length(shared_iid), " samples matched and will be used.")
 
    # Re-order pheno to match .fam order for the shared samples
    pheno <- pheno[match(shared_iid, pheno_iid)]
 
    # Rebuild covariate matrix on the aligned subset
    covar_matrix <- .build_covar_matrix(pheno, all_covar_cols)
 
    # Gsub.id: IIDs of the matched samples, in .fam order
    Gsub.id <- shared_iid
 
  } else {
 
    # No ID column: assume rows are already aligned; pass all .fam IIDs
    if (nrow(pheno) != nrow(fam)) {
      stop("pheno_id is NULL but the phenotype file has ", nrow(pheno),
           " rows while the .fam file has ", nrow(fam), " rows. ",
           "Either supply pheno_id or ensure the two files have the same number ",
           "of rows in the same sample order.")
    }
    Gsub.id <- NULL
  }
 
  rm(fam)

  # ---------------------------------------------------------------------------
  # Fit null model
  # ---------------------------------------------------------------------------

  message("Fitting null model (sample_uncorrelated = ", sample_uncorrelated, ") ...")

  if (sample_uncorrelated) {

    null_model_args <- list(
      Y        = pheno[[phenotype]],
      X        = covar_matrix,
      out_type = out_type
    )
    if (!is.null(pheno_id)) {
      null_model_args$id <- pheno[[pheno_id]]
    }
    nullobj <- do.call(Fit_null_model, null_model_args)

  } else {

    nullobj <- Fit_null_model_GLMM(
      geno_file,
      pheno_file,
      phenotype,
      plink_path,
      outcome_type = out_type,
      sample_id_col = pheno_id,
      covar_cols = covar_cols,
      cat_covar_cols = cat_covar_cols,
      sparse_grm_file = grm_file,
      sparse_grm_id_file = grm_id_file
    )
  }

  rm(pheno); gc()

  # ===========================================================================
  # Single-Window analysis
  # ===========================================================================

  if (test_type == "Single_Window") {

    block_filename <- if (genome_build == "hg19") {
      "LAVA_s2500_m25_f1_w200.blocks"
    } else {
      "deCODE_EUR_LD_blocks.bed"
    }

    block_file <- file.path(
      system.file("extdata", package = "KnockoffPipeline"),
      block_filename
    )

    if (!file.exists(block_file)) {
      stop("LD block reference file not found: ", block_file)
    }

    blocks <- data.table::fread(block_file)
    unique_chr <- intersect(sort(unique(blocks$chr)), chr_vector)

    if (length(unique_chr) == 0L) {
      stop("No chromosomes remain after intersecting block file with requested chromosomes.")
    }

    for (c in unique_chr) {

      message("--- Processing chromosome ", c, " (Single_Window) ---")

      single_mid_file <- file.path(mid_dir, paste0("Single_mid_results_chr", c, ".txt"))
      window_mid_file <- file.path(mid_dir, paste0("Window_mid_results_chr", c, ".txt"))

      if (read_mid_exist && file.exists(single_mid_file) && file.exists(window_mid_file)) {
        message("Intermediate files found; skipping chromosome ", c)
        next
      }

      block_chr <- blocks[blocks$chr==c]
      sub_seq_id <- seq_len(nrow(block_chr))
      
      safe_run_block <- function(kk) {
        tryCatch(
          run_single_block(
            blocks = block_chr,
            kk = kk,
            geno.file = geno_file,
            obj_nullmodel = nullobj,
            window_length = sliding_window_length,
            plink_prefix = plink_path,
            impute.method = geno_missing_imputation,
            M = M,
            Gsub.id = Gsub.id
          ),
          error = function(e) {
            NULL
          }
        )
      }

      out <- parallel::mclapply(sub_seq_id, safe_run_block, mc.cores=user_cores)
      out <- Filter(Negate(is.null), out)

      if (length(out) == 0L) {
        warning("All blocks failed for chromosome ", c, "; skipping.")
        next
      }
      print(str(out))

      single_chr <- data.table::rbindlist(
        lapply(out, function(x) data.table::as.data.table(x$result.single)),
        fill = TRUE
      )
      window_chr <- data.table::rbindlist(
        lapply(out, function(x) data.table::as.data.table(x$result.window)),
        fill = TRUE
      )
      data.table::fwrite(single_chr, single_mid_file, sep = "\t")
      data.table::fwrite(window_chr, window_mid_file, sep = "\t")

      rm(out,single_chr,window_chr); gc()
    }

    # ---- Merge and summarise ------------------------------------------------

    message("Merging intermediate results across chromosomes ...")

    .read_mid <- function(prefix, chr) {
      f <- file.path(mid_dir, paste0(prefix, "_mid_results_chr", chr, ".txt"))
      if (!file.exists(f)) {
        warning("Intermediate file missing, skipping: ", f)
        return(NULL)
      }
      data.table::fread(f)
    }

    result.single.all <- data.table::rbindlist(
      lapply(unique_chr, .read_mid, prefix = "Single"), fill = TRUE
    )
    result.window.all <- data.table::rbindlist(
      lapply(unique_chr, .read_mid, prefix = "Window"), fill = TRUE
    )

    if (nrow(result.single.all) == 0L || nrow(result.window.all) == 0L) {
      stop("No results were produced across all chromosomes. ",
           "Check intermediate files in: ", mid_dir)
    }

    summary <- KS_summary(
      as.matrix(result.window.all),
      as.matrix(result.single.all),
      M,
      fdr=fdr
    )

    keep_cols  <- c("chr", "start", "end", "Qvalue", "W_KS", "W_Threshold", "detect")
    result_all <- summary_res[, keep_cols, with = FALSE]

    out_file <- file.path(outdir, "Single_Window_results.csv")
    data.table::fwrite(result_all, out_file)
    message("Single_Window results written to: ", out_file)

    plot_manhattan(result_all, outdir, "manhattan_plot_single.png")

  }

  # ===========================================================================
  # Gene-Centric analysis
  # ===========================================================================

  if (test_type == "Gene_Centric") {

    gene_file <- file.path(
      system.file("extdata",package="KnockoffPipeline"),
      genome_build,
      paste0("coding.genes.TSS.",genome_build,".tsv")
    )

    if (!file.exists(gene_file)) {
      stop("Gene annotation file not found: ", gene_file)
    }

    genes_info <- data.table::fread(gene_file)
    genes_info$chr <- as.numeric(gsub("[^0-9]","",genes_info$chr))
    genes_info      <- genes_info[!is.na(chr)]
    unique_chr <- intersect(sort(unique(genes_info$chr)), chr_vector)

    if (length(unique_chr) == 0L) {
      stop("No chromosomes remain after intersecting gene file with requested chromosomes.")
    }

    for (c in unique_chr) {

      message("--- Processing chromosome ", c, " (Gene_Centric) ---")

      mid_file_chr <- file.path(mid_dir, paste0("GeneCentric_mid_results_chr", c, ".txt"))

      if (read_mid_exist && file.exists(mid_file_chr)) {
        message("Intermediate file found; skipping chromosome ", c)
        next
      }

      chr_genes <- genes_info[chr==c]
      n_gene <- nrow(chr_genes)

      # ---- Load enhancer annotation for this chromosome --------------------

      abc_file_chr <- file.path(
        system.file("extdata",package="KnockoffPipeline"),
        genome_build,
        paste0("ABC_combined_chr", c, ".csv")
      )
      gh_file_chr <- file.path(
        system.file("extdata",package="KnockoffPipeline"),
        genome_build,
        paste0("GH.data_chr", c, ".csv")
      )
      if (!file.exists(abc_file_chr)) {
        stop("ABC enhancer file not found for chr", c, ": ", abc_file_chr)
      }
      if (!file.exists(gh_file_chr)) {
        stop("GH enhancer file not found for chr", c, ": ", gh_file_chr)
      }
      abc_df <- data.table::fread(abc_file_chr)
      gh_df  <- data.table::fread(gh_file_chr)
      
      # ---- Batch processing ------------------------------------------------

      batch_index <- split(seq_len(n_gene), ceiling(seq_len(n_gene) / batch_size))
      result_list_chr <- vector("list", length(batch_index))

      for (b in seq_along(batch_index)) {

        message("  Batch ", b, " / ", length(batch_index),
                " (genes ", batch_index[[b]][1L], "-",
                batch_index[[b]][length(batch_index[[b]])], ")")

        result_list_chr[[b]] <- run_batch_gene(
          genes = chr_genes,
          kk_vec = batch_index[[b]],
          geno.file = geno_file,
          obj_nullmodel = if(!sample_uncorrelated) nullobj$result.null.model.GLMM else nullobj,
          window_length = sliding_window_length,
          plink_prefix = plink_path,
          M = M,
          genome_build = genome_build,
          Gsub.id = Gsub.id,
          sparseSigma = if(!sample_uncorrelated) nullobj$sparseSigma else NULL,
          ratio = if(!sample_uncorrelated) nullobj$ratio else NULL,
          user_cores = user_cores
        )
        gc()
      }

      result_chr <- data.table::rbindlist(data.table::as.data.table(result_list_chr), fill=TRUE)
      data.table::fwrite(result_chr, mid_file_chr, sep="\t")

      rm(result_list_chr, result_chr); gc()
    }

    # ---- Merge and summarise ------------------------------------------------

    result.all <- data.table::rbindlist(
      lapply(unique_chr, function(c) {
        f <- file.path(mid_dir, paste0("GeneCentric_mid_results_chr", c, ".txt"))
        if (!file.exists(f)) {
          warning("Intermediate file missing, skipping: ", f)
          return(NULL)
        }
        data.table::fread(f)
      }),
      fill = TRUE
    )

    if (nrow(result.all) == 0L) {
      stop("No gene-centric results were produced. ",
           "Check intermediate files in: ", mid_dir)
    }

    summary_res <- GeneScan3DKnock_Summary(result.all, M = M, fdr = fdr)

    keep_cols  <- c("chr", "gene_id", "gene_start", "gene_end",
                    "Qvalue", "W", "W_Threshold", "detect")
    result_all <- summary_res[, keep_cols, with = FALSE]

    data.table::setnames(
      result_all,
      old = c("gene_start", "gene_end", "W"),
      new = c("start",      "end",      "W_KS")
    )

    out_file <- file.path(outdir, "GeneCentric_results.csv")
    data.table::fwrite(result_all, out_file)
    message("Gene_Centric results written to: ", out_file)

    plot_manhattan(result_all, outdir, "manhattan_plot_gene.png")
  }

  invisible(TRUE)
}
