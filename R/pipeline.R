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
#' @param phenotype     Character \strong{vector} of phenotype column name(s).
#'   When multiple phenotypes are provided, samples with missing values in
#'   \emph{any} phenotype or covariate are removed once before any analysis,
#'   knockoffs are generated on the first pass and automatically reused for
#'   subsequent phenotypes (requires \code{save_knockoff = TRUE} or
#'   \code{NULL}).
#' @param pheno_id      Character or \code{NULL}. Column name of the sample ID
#'   in the phenotype file. \code{NULL} assumes rows are already aligned with
#'   the PLINK \code{.fam} file.
#' @param covar_cols    Character vector of continuous covariate column names,
#'   or \code{NULL}.
#' @param cat_covar_cols Character vector of binary/categorical covariate column
#'   names, or \code{NULL}.
#' @param user_cores    Integer. Number of parallel cores. Default \code{1}.
#' @param sliding_window_length Integer or \code{NULL}. Sliding window size (bp).
#' @param geno_missing_imputation Character. Genotype imputation method.
#'   Default \code{"fixed"}.
#' @param plink_path    Character. Path to PLINK executable. Default
#'   \code{"plink"}.
#' @param M             Integer. Number of knockoff copies. Default \code{5}.
#' @param genome_build  Character. One of \code{"hg19"} or \code{"hg38"}.
#' @param sample_uncorrelated Logical. \code{TRUE} fits a standard null model;
#'   \code{FALSE} fits a GLMM via SAIGE.
#' @param grm_file      Character or \code{NULL}. Path to the sparse GRM file.
#' @param grm_id_file   Character or \code{NULL}. Path to the sparse GRM ID
#'   file.
#' @param fdr           Numeric in \code{(0, 1)}. Target FDR level. Default
#'   \code{0.1}.
#' @param chromosomes   Integer vector of autosomes to analyse. Default
#'   \code{1:22}.
#' @param batch_size    Integer. Genes per batch in gene-centric mode. Default
#'   \code{20}.
#' @param read_mid_exist Logical. Skip chromosomes with existing intermediate
#'   files. Default \code{TRUE}.
#' @param pipeline_stage Character. Controls the run mode:
#'   \describe{
#'     \item{\code{"full"}}{(Default) Complete end-to-end pipeline.}
#'     \item{\code{"stage1_knockoff"}}{Generate and save knockoffs only; no
#'       association testing. Writes a sample-list file for reproducibility.
#'       Implies \code{save_knockoff = TRUE}.}
#'     \item{\code{"stage2_analysis"}}{Load pre-generated knockoffs and run
#'       association testing only. Requires \code{knockoff_dir}.
#'       Sample overlap between the saved knockoffs and the current
#'       genotype/phenotype data is validated automatically, with warnings
#'       for any discrepancies.}
#'   }
#' @param save_knockoff Logical or \code{NULL}. Whether to save generated
#'   knockoffs to \code{knockoff_dir}.
#'   \code{NULL} (default) sets this automatically: \code{TRUE} when
#'   \code{pipeline_stage = "stage1_knockoff"} or when multiple phenotypes
#'   are provided; \code{FALSE} otherwise.
#' @param knockoff_dir  Character or \code{NULL}. Directory for knockoff
#'   \code{.rds} files.  Defaults to \code{<outdir>/knockoffs}.  Must be
#'   provided (and populated) when \code{pipeline_stage =
#'   "stage2_analysis"}.
#'
#' @details
#' \strong{Knockoff file format.}
#' Each knockoff file is an RDS list with fields:
#' \code{G_k} (knockoff matrices), \code{sample_ids} (IIDs in row order),
#' and \code{snp_pos} (SNP positions, for validation).
#'
#' \strong{Two-stage workflow.}
#' Run stage 1 to pre-generate knockoffs, then run stage 2 with any
#' phenotype.  Stage 2 reads the sample list produced by stage 1, computes
#' the intersection with the current phenotype/genotype data, reindexes
#' knockoff rows accordingly, and issues warnings for any lost samples.
#'
#' \strong{Multiple phenotypes.}
#' Supply a character vector.  Samples with missing values in \emph{any}
#' phenotype or covariate are removed once.  Per-phenotype output goes to
#' \code{<outdir>/<phenotype_name>/}.
#'
#' @return Invisibly returns \code{TRUE} on success.
#'
#' @import SKAT Matrix WGScan SPAtest CompQuadForm irlba bigmemory
#' @import data.table dplyr parallel qqman abind SAIGE
#' @export
run_pipeline <- function(
  outdir,
  test_type,
  pheno_file,
  geno_file,
  phenotype,
  pheno_id                = NULL,
  covar_cols              = NULL,
  cat_covar_cols          = NULL,
  user_cores              = 1L,
  sliding_window_length   = NULL,
  geno_missing_imputation = "fixed",
  plink_path              = "plink",
  M                       = 5L,
  genome_build            = "hg19",
  sample_uncorrelated     = TRUE,
  grm_file                = NULL,
  grm_id_file             = NULL,
  fdr                     = 0.1,
  chromosomes             = 1:22,
  batch_size              = 20L,
  read_mid_exist          = TRUE,
  pipeline_stage          = "full",
  save_knockoff           = NULL,
  knockoff_dir            = NULL
) {

  # ---------------------------------------------------------------------------
  # 1.  Input validation
  # ---------------------------------------------------------------------------

  stopifnot(
    "outdir must be a single non-empty string"        = is.character(outdir)     && length(outdir)     == 1L && nzchar(outdir),
    "pheno_file must be a single non-empty string"    = is.character(pheno_file) && length(pheno_file) == 1L && nzchar(pheno_file),
    "geno_file must be a single non-empty string"     = is.character(geno_file)  && length(geno_file)  == 1L && nzchar(geno_file),
    "phenotype must be a non-empty character vector"  = is.character(phenotype)  && length(phenotype)  >= 1L && all(nzchar(phenotype)),
    "M must be a positive integer"                    = is.numeric(M)            && length(M)          == 1L && M >= 1L,
    "fdr must be numeric in (0, 1)"                   = is.numeric(fdr)          && length(fdr)        == 1L && fdr > 0 && fdr < 1,
    "user_cores must be a positive integer"           = is.numeric(user_cores)   && length(user_cores) == 1L && user_cores >= 1L,
    "batch_size must be a positive integer"           = is.numeric(batch_size)   && length(batch_size) == 1L && batch_size >= 1L,
    "sample_uncorrelated must be logical"             = is.logical(sample_uncorrelated) && length(sample_uncorrelated) == 1L,
    "read_mid_exist must be logical"                  = is.logical(read_mid_exist) && length(read_mid_exist) == 1L
  )

  if (!test_type %in% c("Single_Window", "Gene_Centric"))
    stop("'test_type' must be \"Single_Window\" or \"Gene_Centric\".")
  if (!genome_build %in% c("hg19", "hg38"))
    stop("'genome_build' must be \"hg19\" or \"hg38\".")
  if (!pipeline_stage %in% c("full", "stage1_knockoff", "stage2_analysis"))
    stop("'pipeline_stage' must be one of \"full\", \"stage1_knockoff\", \"stage2_analysis\".")
  if (!is.null(covar_cols)     && !is.character(covar_cols))     stop("'covar_cols' must be a character vector or NULL.")
  if (!is.null(cat_covar_cols) && !is.character(cat_covar_cols)) stop("'cat_covar_cols' must be a character vector or NULL.")
  if (!is.null(pheno_id) && (!is.character(pheno_id) || length(pheno_id) != 1L))
    stop("'pheno_id' must be a single string or NULL.")
  if (!is.null(save_knockoff) && !is.logical(save_knockoff))
    stop("'save_knockoff' must be TRUE, FALSE, or NULL.")

  chr_numeric <- suppressWarnings(as.integer(chromosomes))
  if (any(is.na(chr_numeric)))
    stop("'chromosomes' contains non-integer values: ", paste(chromosomes[is.na(chr_numeric)], collapse = ", "))
  chr_vector  <- intersect(chr_numeric, 1L:22L)
  if (length(chr_vector) == 0L) stop("No valid autosomes (1-22) in 'chromosomes'.")

  if (!file.exists(pheno_file))
    stop("Phenotype file not found: ", pheno_file)
  plink_fam <- paste0(geno_file, ".fam")
  if (!file.exists(plink_fam))
    stop("PLINK .fam not found: ", plink_fam, "\nCheck that 'geno_file' is the correct prefix.")

  # ---------------------------------------------------------------------------
  # 2.  Resolve save/load flags
  # ---------------------------------------------------------------------------

  multi_pheno <- length(phenotype) > 1L

  # Auto-resolve save_knockoff
  if (is.null(save_knockoff))
    save_knockoff <- (pipeline_stage == "stage1_knockoff") || multi_pheno

  if (pipeline_stage == "stage1_knockoff" && !isTRUE(save_knockoff))
    stop("pipeline_stage = 'stage1_knockoff' requires save_knockoff = TRUE (or NULL).")
  if (pipeline_stage == "stage2_analysis" && is.null(knockoff_dir))
    stop("pipeline_stage = 'stage2_analysis' requires 'knockoff_dir' to be specified.")
  if (multi_pheno && !isTRUE(save_knockoff))
    warning("Multiple phenotypes provided but save_knockoff = FALSE. ",
            "Knockoffs will be regenerated for each phenotype independently, ",
            "which is less efficient and breaks exchangeability across phenotypes.")

  # Resolve knockoff directory
  if (is.null(knockoff_dir)) knockoff_dir <- file.path(outdir, "knockoffs")
  knockoff_sample_file <- file.path(knockoff_dir, "knockoff_sample_list.txt")

  # ---------------------------------------------------------------------------
  # 3.  Create directories
  # ---------------------------------------------------------------------------

  if (!dir.exists(outdir)) { message("Creating output directory: ", outdir); dir.create(outdir, recursive = TRUE) }
  if (isTRUE(save_knockoff) && !dir.exists(knockoff_dir)) {
    message("Creating knockoff directory: ", knockoff_dir); dir.create(knockoff_dir, recursive = TRUE)
  }
  if (user_cores > 1L) Sys.setenv(MKL_NUM_THREADS = 1)

  # ---------------------------------------------------------------------------
  # 4.  Load phenotype file; remove samples missing in ANY phenotype / covariate
  # ---------------------------------------------------------------------------

  message("Reading phenotype file: ", pheno_file)
  pheno <- data.table::fread(pheno_file)

  missing_pheno <- setdiff(phenotype, colnames(pheno))
  if (length(missing_pheno) > 0L)
    stop("Phenotype column(s) not found: ", paste(missing_pheno, collapse = ", "))

  all_covar_cols <- c(covar_cols, cat_covar_cols)
  missing_covar  <- setdiff(all_covar_cols, colnames(pheno))
  if (length(missing_covar) > 0L)
    stop("Covariate column(s) not found: ", paste(missing_covar, collapse = ", "))
  if (!is.null(pheno_id) && !pheno_id %in% colnames(pheno))
    stop("Sample ID column \"", pheno_id, "\" not found in phenotype file.")

  # Drop rows missing in ANY phenotype or covariate (single consistent sample set)
  check_cols    <- unique(c(phenotype, all_covar_cols))
  complete_mask <- complete.cases(pheno[, check_cols, with = FALSE])
  n_incomplete  <- sum(!complete_mask)
  if (n_incomplete > 0L) {
    message(sprintf(
      "%d sample(s) removed: missing in at least one of [%s].",
      n_incomplete, paste(check_cols, collapse = ", ")
    ))
    pheno <- pheno[complete_mask]
  }
  message(nrow(pheno), " sample(s) retained after missing-value filtering.")

  # ---------------------------------------------------------------------------
  # 5.  Align phenotype file to PLINK .fam
  # ---------------------------------------------------------------------------

  fam <- data.table::fread(plink_fam, header = FALSE,
                           col.names = c("FID","IID","PAT","MAT","SEX","PHENO"))

  if (!is.null(pheno_id)) {
    pheno_iid  <- as.numeric(pheno[[pheno_id]])
    fam_iid    <- as.numeric(fam$IID)
    shared_iid <- intersect(fam_iid, pheno_iid)   # keeps .fam order

    if (length(shared_iid) == 0L)
      stop("No samples matched between phenotype (column \"", pheno_id,
           "\") and PLINK .fam.\n",
           "  Example pheno IID : ", paste(head(pheno_iid, 3L), collapse = ", "), "\n",
           "  Example .fam  IID : ", paste(head(fam_iid,   3L), collapse = ", "))

    n_pheno_only <- length(setdiff(pheno_iid, fam_iid))
    n_fam_only   <- length(setdiff(fam_iid,   pheno_iid))
    if (n_pheno_only > 0L) message("  ", n_pheno_only, " sample(s) in phenotype not in .fam — excluded.")
    if (n_fam_only   > 0L) message("  ", n_fam_only,   " sample(s) in .fam not in phenotype — excluded.")
    message("  ", length(shared_iid), " sample(s) matched.")

    pheno   <- pheno[match(shared_iid, pheno_iid)]
    Gsub.id <- shared_iid

  } else {
    if (nrow(pheno) != nrow(fam))
      stop("pheno_id is NULL but phenotype has ", nrow(pheno),
           " rows while .fam has ", nrow(fam), " rows.")
    Gsub.id <- NULL
  }
  rm(fam); gc()

  # ---------------------------------------------------------------------------
  # 6.  Stage 2: validate saved knockoff sample list; restrict sample set
  # ---------------------------------------------------------------------------

  # The 'canonical' IDs used for knockoff row-ordering are what we call
  # knockoff_sample_ids.  This is written during stage1 / multi-pheno pass 1
  # and read back during stage2 / multi-pheno pass ≥2.

  if (pipeline_stage == "stage2_analysis") {
    if (!file.exists(knockoff_sample_file))
      stop("Knockoff sample list not found: ", knockoff_sample_file,
           "\nRun pipeline_stage = 'stage1_knockoff' first.")

    saved_ids   <- as.numeric(readLines(knockoff_sample_file))
    current_ids <- if (!is.null(Gsub.id)) Gsub.id else seq_len(nrow(pheno))
    common_ids  <- intersect(saved_ids, current_ids)   # preserves saved order

    n_only_saved   <- length(setdiff(saved_ids,   current_ids))
    n_only_current <- length(setdiff(current_ids, saved_ids))

    if (n_only_saved > 0L)
      warning(n_only_saved, " sample(s) in saved knockoffs are absent from the ",
              "current genotype/phenotype — those knockoff rows will be ignored.")
    if (n_only_current > 0L)
      warning(n_only_current, " sample(s) present in the current data are absent ",
              "from the saved knockoffs — they will be EXCLUDED from this analysis.")
    if (length(common_ids) == 0L)
      stop("No samples shared between saved knockoffs and current dataset.")

    message(length(common_ids), " sample(s) shared between saved knockoffs and current data.")

    # Restrict pheno to the common set in saved-ids order for row alignment
    if (!is.null(Gsub.id)) {
      pheno   <- pheno[match(common_ids, as.numeric(pheno[[pheno_id]]))]
      Gsub.id <- common_ids
    }
    knockoff_sample_ids <- saved_ids   # full saved list (needed inside load helpers)

  } else {
    # stage1 or full: determine and (optionally) write the canonical IDs
    knockoff_sample_ids <- if (!is.null(Gsub.id)) Gsub.id else seq_len(nrow(pheno))
    if (isTRUE(save_knockoff)) {
      writeLines(as.character(knockoff_sample_ids), knockoff_sample_file)
      message("Knockoff sample list written to: ", knockoff_sample_file)
    }
  }

  # ---------------------------------------------------------------------------
  # 7.  Stage 1: knockoff generation only (no null model, no tests)
  # ---------------------------------------------------------------------------

  if (pipeline_stage == "stage1_knockoff") {
    message("\n======= Stage 1: Knockoff Generation =======")
    .run_knockoff_generation(
      test_type               = test_type,
      geno_file               = geno_file,
      Gsub.id                 = Gsub.id,
      knockoff_dir            = knockoff_dir,
      chr_vector              = chr_vector,
      M                       = M,
      genome_build            = genome_build,
      sliding_window_length   = sliding_window_length,
      geno_missing_imputation = geno_missing_imputation,
      plink_path              = plink_path,
      batch_size              = batch_size,
      sample_uncorrelated     = sample_uncorrelated,
      user_cores              = user_cores,
      read_mid_exist          = read_mid_exist
    )
    message("\nStage 1 complete.")
    message("  Knockoffs : ", knockoff_dir)
    message("  Sample list: ", knockoff_sample_file)
    return(invisible(TRUE))
  }

  # ===========================================================================
  # 8.  Stages "full" / "stage2_analysis": loop over phenotype(s)
  # ===========================================================================

  # For multi-phenotype "full" runs:
  #   pass 1 (first phenotype)  → save_knockoff = TRUE  (generate + save)
  #   pass ≥2 (later phenotypes) → load_knockoff = TRUE  (load saved)
  # For "stage2_analysis": always load.

  knockoffs_ready <- (pipeline_stage == "stage2_analysis")

  for (pheno_name in phenotype) {

    message("\n======= Phenotype: \"", pheno_name, "\" =======")

    # Per-phenotype subdirectory (only created when multiple phenotypes)
    p_outdir  <- if (multi_pheno) file.path(outdir, pheno_name) else outdir
    p_mid_dir <- file.path(p_outdir, "mid")
    for (d in c(p_outdir, p_mid_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

    # ---- Outcome type -------------------------------------------------------
    pv       <- unique(pheno[[pheno_name]])
    is_bin   <- length(pv) == 2L && all(sort(as.numeric(pv)) == c(0, 1))
    out_type <- if (is_bin) "D" else "C"
    message("Outcome type: ", if (is_bin) "binary (D)" else "continuous (C)")

    # ---- Covariate matrix ---------------------------------------------------
    covar_matrix <- if (length(all_covar_cols) > 0L)
      as.matrix(pheno[, all_covar_cols, with = FALSE]) else NULL

    # ---- Fit null model -----------------------------------------------------
    message("Fitting null model (sample_uncorrelated = ", sample_uncorrelated, ") ...")

    if (sample_uncorrelated) {
      nm_args <- list(Y = pheno[[pheno_name]], X = covar_matrix, out_type = out_type)
      if (!is.null(pheno_id)) nm_args$id <- pheno[[pheno_id]]
      nullobj <- do.call(Fit_null_model, nm_args)
    } else {
      nullobj <- Fit_null_model_GLMM(
        geno_file, pheno_file, pheno_name, plink_path,
        outcome_type       = out_type,
        sample_id_col      = pheno_id,
        covar_cols         = covar_cols,
        cat_covar_cols     = cat_covar_cols,
        sparse_grm_file    = grm_file,
        sparse_grm_id_file = grm_id_file
      )
    }

    # ---- Knockoff flags for this pass ---------------------------------------
    ko_save <- isTRUE(save_knockoff) && !knockoffs_ready
    ko_load <- isTRUE(save_knockoff) &&  knockoffs_ready

    # ---- Branch by test type ------------------------------------------------
    if (test_type == "Single_Window") {

      .run_single_window(
        outdir                  = p_outdir,
        mid_dir                 = p_mid_dir,
        geno_file               = geno_file,
        nullobj                 = nullobj,
        Gsub.id                 = Gsub.id,
        chr_vector              = chr_vector,
        genome_build            = genome_build,
        sliding_window_length   = sliding_window_length,
        geno_missing_imputation = geno_missing_imputation,
        plink_path              = plink_path,
        M                       = M,
        user_cores              = user_cores,
        read_mid_exist          = read_mid_exist,
        fdr                     = fdr,
        save_knockoff           = ko_save,
        load_knockoff           = ko_load,
        knockoff_dir            = knockoff_dir,
        knockoff_sample_ids     = knockoff_sample_ids
      )

    } else {

      .run_gene_centric(
        outdir                  = p_outdir,
        mid_dir                 = p_mid_dir,
        geno_file               = geno_file,
        nullobj                 = nullobj,
        Gsub.id                 = Gsub.id,
        chr_vector              = chr_vector,
        genome_build            = genome_build,
        sliding_window_length   = sliding_window_length,
        plink_path              = plink_path,
        M                       = M,
        user_cores              = user_cores,
        read_mid_exist          = read_mid_exist,
        fdr                     = fdr,
        batch_size              = batch_size,
        sample_uncorrelated     = sample_uncorrelated,
        save_knockoff           = ko_save,
        load_knockoff           = ko_load,
        knockoff_dir            = knockoff_dir,
        knockoff_sample_ids     = knockoff_sample_ids,
        sparseSigma             = if (!sample_uncorrelated) nullobj$sparseSigma else NULL,
        ratio                   = if (!sample_uncorrelated) nullobj$ratio       else NULL
      )
    }

    knockoffs_ready <- TRUE   # knockoffs now exist on disk for subsequent phenotypes
    gc()
  }

  invisible(TRUE)
}


# =============================================================================
# Internal: Stage-1 knockoff generation only (no association testing)
# Dispatches to Single_Window or Gene_Centric generation helpers.
# =============================================================================

.run_knockoff_generation <- function(
  test_type, geno_file, Gsub.id, knockoff_dir, chr_vector,
  M, genome_build, sliding_window_length, geno_missing_imputation,
  plink_path, batch_size, sample_uncorrelated, user_cores, read_mid_exist
) {
  # A minimal "null object" is not needed here: run_single_block /
  # run_batch_gene accept save_knockoff = TRUE without running tests.
  # We pass nullobj = NULL and the analysis branches will short-circuit.

  if (test_type == "Single_Window") {
    block_filename <- if (genome_build == "hg19") "LAVA_s2500_m25_f1_w200.blocks" else "deCODE_EUR_LD_blocks.bed"
    block_file     <- file.path(system.file("extdata", package = "KnockoffPipeline"), block_filename)
    if (!file.exists(block_file)) stop("LD block reference not found: ", block_file)

    blocks     <- data.table::fread(block_file)
    unique_chr <- intersect(sort(unique(blocks$chr)), chr_vector)
    if (length(unique_chr) == 0L) stop("No chromosomes after intersecting block file with requested chromosomes.")

    for (c in unique_chr) {
      message("--- chr ", c, " ---")
      chr_ko_dir <- file.path(knockoff_dir, paste0("chr", c))
      if (!dir.exists(chr_ko_dir)) dir.create(chr_ko_dir)
      block_chr  <- blocks[blocks$chr == c]

      parallel::mclapply(seq_len(nrow(block_chr)), function(kk) {
        ko_file <- .ko_file_single(chr_ko_dir, kk)
        if (read_mid_exist && file.exists(ko_file)) return(invisible(NULL))
        run_single_block(
          blocks                  = block_chr,
          kk                      = kk,
          geno.file               = geno_file,
          obj_nullmodel           = NULL,   # no null model needed
          window_length           = sliding_window_length,
          plink_prefix            = plink_path,
          impute.method           = geno_missing_imputation,
          M                       = M,
          Gsub.id                 = Gsub.id,
          save_knockoff           = TRUE,
          load_knockoff           = FALSE,
          knockoff_file           = ko_file,
          knockoff_sample_ids     = Gsub.id,
          stage1_only             = TRUE    # skip association test
        )
      }, mc.cores = user_cores)
    }

  } else {
    # Gene_Centric
    gene_file <- file.path(system.file("extdata", package = "KnockoffPipeline"),
                           genome_build, paste0("coding.genes.TSS.", genome_build, ".tsv"))
    if (!file.exists(gene_file)) stop("Gene annotation file not found: ", gene_file)

    genes_info <- data.table::fread(gene_file)
    genes_info$chr <- as.numeric(gsub("[^0-9]", "", genes_info$chr))
    genes_info      <- genes_info[!is.na(chr)]
    unique_chr      <- intersect(sort(unique(genes_info$chr)), chr_vector)
    if (length(unique_chr) == 0L) stop("No chromosomes remain.")

    for (c in unique_chr) {
      message("--- chr ", c, " ---")
      chr_ko_dir <- file.path(knockoff_dir, paste0("chr", c))
      if (!dir.exists(chr_ko_dir)) dir.create(chr_ko_dir)

      chr_genes <- genes_info[chr == c]
      abc_df    <- data.table::fread(.extdata_path(genome_build, paste0("ABC_combined_chr",   c, ".csv")))
      gh_df     <- data.table::fread(.extdata_path(genome_build, paste0("GH.data_chr",        c, ".csv")))

      batch_index <- split(seq_len(nrow(chr_genes)), ceiling(seq_len(nrow(chr_genes)) / batch_size))
      for (b in seq_along(batch_index)) {
        message("  Batch ", b, " / ", length(batch_index))
        run_batch_gene(
          genes               = chr_genes,
          kk_vec              = batch_index[[b]],
          geno.file           = geno_file,
          obj_nullmodel       = NULL,
          window_length       = sliding_window_length,
          plink_prefix        = plink_path,
          M                   = M,
          genome_build        = genome_build,
          Gsub.id             = Gsub.id,
          abc_df              = abc_df,
          gh_df               = gh_df,
          sparseSigma         = NULL,
          ratio               = NULL,
          user_cores          = user_cores,
          save_knockoff       = TRUE,
          load_knockoff       = FALSE,
          knockoff_dir        = chr_ko_dir,
          knockoff_sample_ids = Gsub.id,
          stage1_only         = TRUE,
          read_mid_exist      = read_mid_exist
        )
        gc()
      }
    }
  }
}


# =============================================================================
# Internal: Single_Window analysis for one phenotype
# =============================================================================

.run_single_window <- function(
  outdir, mid_dir, geno_file, nullobj, Gsub.id, chr_vector, genome_build,
  sliding_window_length, geno_missing_imputation, plink_path, M,
  user_cores, read_mid_exist, fdr,
  save_knockoff, load_knockoff, knockoff_dir, knockoff_sample_ids
) {
  block_filename <- if (genome_build == "hg19") "LAVA_s2500_m25_f1_w200.blocks" else "deCODE_EUR_LD_blocks.bed"
  block_file     <- file.path(system.file("extdata", package = "KnockoffPipeline"), block_filename)
  if (!file.exists(block_file)) stop("LD block reference not found: ", block_file)

  blocks     <- data.table::fread(block_file)
  unique_chr <- intersect(sort(unique(blocks$chr)), chr_vector)
  if (length(unique_chr) == 0L) stop("No chromosomes remain after intersecting block file.")

  for (c in unique_chr) {
    message("--- chr ", c, " (Single_Window) ---")
    single_mid_file <- file.path(mid_dir, paste0("Single_mid_results_chr", c, ".txt"))
    window_mid_file <- file.path(mid_dir, paste0("Window_mid_results_chr", c, ".txt"))

    if (read_mid_exist && file.exists(single_mid_file) && file.exists(window_mid_file)) {
      message("  Existing intermediate files found — skipping chr ", c); next
    }

    block_chr  <- blocks[blocks$chr == c]
    chr_ko_dir <- file.path(knockoff_dir, paste0("chr", c))

    out <- parallel::mclapply(seq_len(nrow(block_chr)), function(kk) {
      tryCatch(
        run_single_block(
          blocks                  = block_chr,
          kk                      = kk,
          geno.file               = geno_file,
          obj_nullmodel           = nullobj,
          window_length           = sliding_window_length,
          plink_prefix            = plink_path,
          impute.method           = geno_missing_imputation,
          M                       = M,
          Gsub.id                 = Gsub.id,
          save_knockoff           = save_knockoff,
          load_knockoff           = load_knockoff,
          knockoff_file           = .ko_file_single(chr_ko_dir, kk),
          knockoff_sample_ids     = knockoff_sample_ids,
          stage1_only             = FALSE
        ),
        error = function(e) { warning("Block ", kk, " chr ", c, " failed: ", conditionMessage(e)); NULL }
      )
    }, mc.cores = user_cores)

    out <- Filter(Negate(is.null), out)
    if (length(out) == 0L) { warning("All blocks failed for chr ", c, " — skipping."); next }

    single_chr <- data.table::rbindlist(lapply(out, function(x) data.table::as.data.table(x$result.single)), fill = TRUE)
    window_chr <- data.table::rbindlist(lapply(out, function(x) data.table::as.data.table(x$result.window)), fill = TRUE)
    data.table::fwrite(single_chr, single_mid_file, sep = "\t")
    data.table::fwrite(window_chr, window_mid_file, sep = "\t")
    rm(out, single_chr, window_chr); gc()
  }

  # ---- Merge and summarise -------------------------------------------------
  message("Merging intermediate results ...")

  read_mid <- function(prefix, chr) {
    f <- file.path(mid_dir, paste0(prefix, "_mid_results_chr", chr, ".txt"))
    if (!file.exists(f)) { warning("Missing intermediate file: ", f); return(NULL) }
    data.table::fread(f)
  }

  result.single.all <- data.table::rbindlist(lapply(unique_chr, read_mid, prefix = "Single"), fill = TRUE)
  result.window.all <- data.table::rbindlist(lapply(unique_chr, read_mid, prefix = "Window"), fill = TRUE)

  if (nrow(result.single.all) == 0L || nrow(result.window.all) == 0L)
    stop("No results across all chromosomes. Check intermediate files in: ", mid_dir)

  # FIX: was `summary <- ...` (name clash + undefined summary_res below)
  summary_res <- KS_summary(
    as.matrix(result.window.all),
    as.matrix(result.single.all),
    M, fdr = fdr
  )

  keep_cols  <- c("chr", "start", "end", "Qvalue", "W_KS", "W_Threshold", "detect")
  result_all <- summary_res[, keep_cols]

  out_file <- file.path(outdir, "Single_Window_results.csv")
  data.table::fwrite(result_all, out_file)
  message("Single_Window results written to: ", out_file)
  plot_manhattan(result_all, outdir, "manhattan_plot_single.png")
}


# =============================================================================
# Internal: Gene_Centric analysis for one phenotype
# =============================================================================

.run_gene_centric <- function(
  outdir, mid_dir, geno_file, nullobj, Gsub.id, chr_vector, genome_build,
  sliding_window_length, plink_path, M, user_cores, read_mid_exist, fdr,
  batch_size, sample_uncorrelated,
  save_knockoff, load_knockoff, knockoff_dir, knockoff_sample_ids,
  sparseSigma, ratio
) {
  gene_file <- file.path(system.file("extdata", package = "KnockoffPipeline"),
                         genome_build, paste0("coding.genes.TSS.", genome_build, ".tsv"))
  if (!file.exists(gene_file)) stop("Gene annotation file not found: ", gene_file)

  genes_info <- data.table::fread(gene_file)
  genes_info$chr <- as.numeric(gsub("[^0-9]", "", genes_info$chr))
  genes_info      <- genes_info[!is.na(chr)]
  unique_chr      <- intersect(sort(unique(genes_info$chr)), chr_vector)
  if (length(unique_chr) == 0L) stop("No chromosomes remain.")

  for (c in unique_chr) {
    message("--- chr ", c, " (Gene_Centric) ---")
    mid_file_chr <- file.path(mid_dir, paste0("GeneCentric_mid_results_chr", c, ".txt"))

    if (read_mid_exist && file.exists(mid_file_chr)) {
      message("  Existing intermediate file found — skipping chr ", c); next
    }

    chr_genes  <- genes_info[chr == c]
    chr_genes <- chr_genes[order(chr_genes$start), ]
    chr_ko_dir <- file.path(knockoff_dir, paste0("chr", c))
    abc_df     <- data.table::fread(.extdata_path(genome_build, paste0("ABC_combined_chr", c, ".csv")))
    gh_df      <- data.table::fread(.extdata_path(genome_build, paste0("GH.data_chr",      c, ".csv")))

    batch_index     <- split(seq_len(nrow(chr_genes)), ceiling(seq_len(nrow(chr_genes)) / batch_size))
    result_list_chr <- vector("list", length(batch_index))

    for (b in seq_along(batch_index)) {

      result_list_chr[[b]] <- run_batch_gene(
        genes               = chr_genes,
        b                   = b,
        batch_index         = batch_index,
        geno.file           = geno_file,
        obj_nullmodel       = if (!sample_uncorrelated) nullobj$result.null.model.GLMM else nullobj,
        window_length       = sliding_window_length,
        plink_prefix        = plink_path,
        M                   = M,
        genome_build        = genome_build,
        Gsub.id             = Gsub.id,
        abc_df              = abc_df,
        gh_df               = gh_df,
        sparseSigma         = sparseSigma,
        ratio               = ratio,
        user_cores          = user_cores,
        save_knockoff       = save_knockoff,
        load_knockoff       = load_knockoff,
        knockoff_dir        = chr_ko_dir,
        knockoff_sample_ids = knockoff_sample_ids,
        stage1_only         = FALSE,
        read_mid_exist      = read_mid_exist
      )
      gc()
    }

    result_chr <- data.table::rbindlist(Filter(Negate(is.null), result_list_chr), fill = TRUE)
    data.table::fwrite(result_chr, mid_file_chr, sep = "\t")
    rm(result_list_chr, result_chr); gc()
  }

  # ---- Merge and summarise -------------------------------------------------
  result.all <- data.table::rbindlist(
    lapply(unique_chr, function(c) {
      f <- file.path(mid_dir, paste0("GeneCentric_mid_results_chr", c, ".txt"))
      if (!file.exists(f)) { warning("Missing intermediate file: ", f); return(NULL) }
      data.table::fread(f)
    }),
    fill = TRUE
  )
  if (nrow(result.all) == 0L)
    stop("No gene-centric results. Check intermediate files in: ", mid_dir)

  summary_res <- GeneScan3DKnock_Summary(result.all, M = M, fdr = fdr)

  keep_cols  <- c("chr", "gene_id", "gene_start", "gene_end", "Qvalue", "W", "W_Threshold", "detect")
  result_all <- summary_res[, keep_cols]
  data.table::setnames(result_all, old = c("gene_start","gene_end","W"), new = c("start","end","W_KS"))

  out_file <- file.path(outdir, "GeneCentric_results.csv")
  data.table::fwrite(result_all, out_file)
  message("Gene_Centric results written to: ", out_file)
  plot_manhattan(result_all, outdir, "manhattan_plot_gene.png")
}


# =============================================================================
# Knockoff file path helpers
# =============================================================================

#' @keywords internal
.ko_file_single <- function(chr_ko_dir, kk)
  file.path(chr_ko_dir, sprintf("block_%04d_knockoff.rds", kk))

#' @keywords internal
.ko_file_gene <- function(chr_ko_dir, gene_id)
  file.path(chr_ko_dir, paste0("gene_", gsub("[^a-zA-Z0-9_.-]", "_", gene_id), "_knockoff.rds"))

#' @keywords internal
.extdata_path <- function(...) {
  p <- file.path(system.file("extdata", package = "KnockoffPipeline"), ...)
  if (!file.exists(p)) stop("Required data file not found: ", p)
  p
}

#' Align knockoff rows to a new sample order
#'
#' @param ko_obj    List returned by \code{saveRDS}: must contain
#'   \code{sample_ids} and either \code{G_k} (Single_Window) or
#'   \code{G_gene_buffer_knockoff} (Gene_Centric).
#' @param target_ids Numeric/character vector of IDs in the desired row order.
#' @return The same list with knockoff matrices reindexed to \code{target_ids}.
#' @keywords internal
.align_knockoff_samples <- function(ko_obj, target_ids) {
  saved_ids <- ko_obj$sample_ids
  row_map   <- match(as.numeric(target_ids), as.numeric(saved_ids))

  if (any(is.na(row_map)))
    warning(sum(is.na(row_map)), " target sample(s) not found in saved knockoff — ",
            "they will be NA in the aligned knockoff matrices.")

  # Single_Window: G_k is a list of M sparse matrices (n × p)
  if (!is.null(ko_obj$G_k)) {
    ko_obj$G_k <- lapply(ko_obj$G_k, function(m) m[row_map, , drop = FALSE])
  }

  # Gene_Centric: knockoffs are M × n × p arrays
  if (!is.null(ko_obj$G_gene_buffer_knockoff)) {
    ko_obj$G_gene_buffer_knockoff <- ko_obj$G_gene_buffer_knockoff[, row_map, , drop = FALSE]
  }
  if (!is.null(ko_obj$G_EnhancerAll_knockoff) && length(dim(ko_obj$G_EnhancerAll_knockoff)) == 3L) {
    ko_obj$G_EnhancerAll_knockoff <- ko_obj$G_EnhancerAll_knockoff[, row_map, , drop = FALSE]
  }

  ko_obj$sample_ids <- target_ids
  ko_obj
}