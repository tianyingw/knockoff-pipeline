#' Run Knockoff Genome-wide Pipeline
#'
#' Main entry function for the KnockoffPipeline package. Supports both
#' single-window and gene-centric association testing with knockoff-based
#' FDR control.
#'
#' @param outdir        Character. Output directory (created recursively if absent).
#' @param test_type     Character. One of \code{"Single_Window"} or
#'   \code{"Gene_Centric"}.
#' @param pheno_file    Character. Path to the phenotype file
#'   (tab- or comma-separated). Required for every pipeline stage so stage 1
#'   and downstream analysis use the same complete-case sample set.
#' @param geno_file     Character. PLINK genotype file prefix (no extension).
#' @param phenotype     Character \strong{vector} of phenotype column name(s).
#'   When multiple phenotypes are provided, samples with missing values in
#'   \emph{any} phenotype or covariate are removed once before any analysis,
#'   reusable knockoffs are generated and persisted internally on the first
#'   pass. SNP/window and gene-buffer knockoffs are reused; gene-centric
#'   enhancer knockoffs are regenerated for each phenotype.
#' @param pheno_id      Character or \code{NULL}. Column name of the sample ID
#'   in the phenotype file. \code{NULL} assumes rows are already aligned with
#'   the PLINK \code{.fam} file.
#' @param covar_cols    Character vector of continuous covariate column names,
#'   or \code{NULL}.
#' @param cat_covar_cols Character vector of binary/categorical covariate column
#'   names, or \code{NULL}.
#' @param user_cores    Integer. Number of parallel cores. Default \code{1}.
#' @param sliding_window_length Integer vector. Sliding window sizes (bp).
#'   Default \code{c(1000, 5000, 10000)}.
#' @param geno_missing_imputation Character. Genotype imputation method for
#'   \code{"Single_Window"}: \code{"fixed"}, \code{"random"}, or
#'   \code{"bestguess"}. Default \code{"fixed"}. Gene-centric analysis
#'   currently uses fixed imputation.
#' @param plink_path    Character. Path to a PLINK 1.9 or PLINK 2 executable.
#'   The additive-export syntax is selected from the executable version.
#'   Default \code{"plink2"}.
#' @param M             Integer. Number of knockoff copies. Default \code{5}.
#' @param seed          Integer or \code{NULL}. Optional base seed. When set,
#'   each chromosome/block or gene receives a deterministic derived seed so
#'   parallel and resumed runs regenerate the same knockoffs. Required when
#'   \code{geno_missing_imputation = "random"}. Supply a seed for
#'   reproducible restart, stage-1/stage-2, or multi-phenotype runs; with
#'   \code{NULL}, independent process restarts are not reproducible.
#' @param genome_build  Character. One of \code{"hg19"} or \code{"hg38"}.
#' @param ld_block_file Character or \code{NULL}. Optional ancestry-matched LD
#'   block file for \code{"Single_Window"}; it must contain \code{chr},
#'   \code{start}, and \code{stop} columns. The bundled European-ancestry
#'   resource for \code{genome_build} is used when \code{NULL}.
#' @param sample_uncorrelated Logical model-path selector. \code{TRUE} fits a
#'   standard GLM and does not test or prune sample relatedness; \code{FALSE}
#'   selects the BIGKnock/SAIGE GLMM path for \code{"Gene_Centric"}.
#'   \code{"Single_Window"} currently accepts \code{TRUE} only.
#' @param grm_file      Character or \code{NULL}. Path to the sparse GRM file.
#' @param grm_id_file   Character or \code{NULL}. Path to the sparse GRM ID
#'   file.
#' @param relatedness_cutoff Numeric. SAIGE sparse-GRM relatedness cutoff.
#'   Default \code{0.125}.
#' @param n_markers_grm Integer. Number of randomly selected markers used when
#'   constructing a sparse GRM. Default \code{1000}.
#' @param fdr           Numeric in \code{(0, 1)}. Target FDR level. Default
#'   \code{0.1}.
#' @param chromosomes   Integer vector of autosomes to analyse. Default
#'   \code{1:22}.
#' @param batch_size    Integer. Genes per batch in gene-centric mode. Default
#'   \code{20}.
#' @param read_mid_exist Logical. Resume from compatible block-, batch-, or
#'   chromosome-level intermediate files when \code{TRUE}. A versioned manifest
#'   must match the current input identities and settings. \code{FALSE} removes
#'   and recreates only the dedicated \code{mid} directory for a fresh run.
#'   Default \code{TRUE}.
#' @param pipeline_stage Character. Controls the run mode:
#'   \describe{
#'     \item{\code{"full"}}{(Default) Complete end-to-end pipeline.}
#'     \item{\code{"stage1_knockoff"}}{Generate and save knockoffs only; no
#'       association testing. First applies the same phenotype/covariate
#'       complete-case and PLINK-ID matching rules as a full run, then writes
#'       a sample-list file for reproducibility. Implies
#'       \code{save_knockoff = TRUE}.}
#'     \item{\code{"stage2_analysis"}}{Load pre-generated knockoffs and run
#'       association testing only. Requires \code{knockoff_dir}.
#'       The saved and current character-IID sets must match exactly (row order
#'       may differ); incompatible sample, variant, or run metadata error.}
#'   }
#' @param save_knockoff Logical or \code{NULL}. Controls whether the knockoff
#'   directory is \strong{retained} after the run completes.
#'   \code{NULL} (default) sets this automatically to \code{TRUE} when
#'   \code{pipeline_stage = "stage1_knockoff"}, and \code{FALSE} otherwise.
#'   For multi-phenotype runs knockoffs are always written to disk internally
#'   during the run (so they can be reused across phenotypes). A temporary
#'   knockoff directory created by the current full run is deleted at the end
#'   when \code{save_knockoff = FALSE}; a pre-existing directory and every
#'   directory used by \code{pipeline_stage = "stage2_analysis"} are never
#'   recursively deleted. Use a dedicated knockoff directory because its
#'   checkpoint files are managed as one unit.
#' @param knockoff_dir  Character or \code{NULL}. Directory for knockoff
#'   files. Defaults to \code{<outdir>/knockoffs}. Must be
#'   provided (and populated) when \code{pipeline_stage =
#'   "stage2_analysis"}.
#' @param temp_dir Character or \code{NULL}. Directory for temporary PLINK
#'   exports. When \code{NULL}, a writable \code{SLURM_TMPDIR} is preferred and
#'   R's \code{tempdir()} is used otherwise. On clusters, choose node-local
#'   scratch rather than a shared output directory to avoid parallel I/O
#'   contention.
#'
#' @details
#' \strong{Knockoff file format.}
#' Each saved unit has an RDS object recording sample IDs, SNP positions,
#' matrix-storage metadata, and a versioned compatibility context.
#' For \code{"Single_Window"}, the RDS is a manifest for one file-backed
#' \pkg{bigmemory} matrix per knockoff copy; its adjacent \code{.desc} and
#' \code{.bin} files must be kept with the manifest. Gene-centric arrays are
#' stored inside their per-gene RDS files.
#'
#' \strong{Saved-knockoff compatibility.}
#' Reuse requires the same character-IID set (row order may differ) and
#' compares an ordered fingerprint containing BIM
#' chromosome, variant ID, base position, both alleles, and the coded allele.
#' The analysis path, \code{M}, genome build, construction settings, and the
#' identity of the LD-block definition file (Single_Window) or gene-annotation
#' file (Gene_Centric) are also compared. Missing or obsolete
#' metadata and any mismatch fail closed with an error.
#'
#' \strong{Analysis checkpoints.}
#' Before intermediate scores are reused, a manifest validates the phenotype,
#' character-IID sample order, PLINK and reference identities, model path, seed,
#' and analysis settings. Legacy intermediates without a manifest and any
#' mismatch fail closed. Setting \code{read_mid_exist = FALSE} starts fresh by
#' replacing only the analysis-specific \code{mid} directory.
#'
#' \strong{Two-stage workflow.}
#' Run stage 1 with the phenotype(s) and covariates that define the analysis
#' complete cases, then run stage 2 with the same resulting sample set. Stage 1
#' uses these columns only to establish the sample set and does not fit a null
#' model. Stage 2 reads the saved sample list, requires an exact character-IID
#' set match, and reindexes rows if their order differs.
#'
#' \strong{Multiple phenotypes.}
#' Supply a character vector. Samples with missing values in \emph{any}
#' phenotype or covariate are removed once. SNP/window and gene-buffer
#' knockoffs are reused, while gene-centric enhancer knockoffs are regenerated
#' for each phenotype. Per-phenotype output goes to
#' \code{<outdir>/<phenotype_name>/}.
#'
#' @return Invisibly returns \code{TRUE} on success.
#'
#' @examples
#' \dontrun{
#' demo_bed <- system.file(
#'   "examples", "input", "demo.bed",
#'   package = "KnockoffPipeline", mustWork = TRUE
#' )
#' demo_geno <- tools::file_path_sans_ext(demo_bed)
#' demo_pheno <- system.file(
#'   "examples", "input", "phenotype.csv",
#'   package = "KnockoffPipeline", mustWork = TRUE
#' )
#' run_pipeline(
#'   outdir = tempfile("KnockoffPipeline-SNP-Window-"),
#'   test_type = "Single_Window", geno_file = demo_geno,
#'   pheno_file = demo_pheno, phenotype = "Y", pheno_id = "IID",
#'   covar_cols = "X1", chromosomes = 1, seed = 20260915L,
#'   read_mid_exist = FALSE
#' )
#' }
#'
#' @import SKAT Matrix WGScan CompQuadForm irlba bigmemory
#' @import data.table parallel qqman abind SAIGE
#' @importFrom graphics abline
#' @importFrom grDevices dev.off png
#' @importFrom stats as.dist binomial complete.cases cutree dbeta end gaussian
#'   glm hclust median pcauchy pchisq rbinom sd start var
#' @importFrom utils capture.output read.table write.table
#' @export
run_pipeline <- function(
  outdir,
  test_type,
  geno_file,
  pheno_file              = NULL,
  phenotype               = NULL,
  pheno_id                = NULL,
  covar_cols              = NULL,
  cat_covar_cols          = NULL,
  user_cores              = 1L,
  sliding_window_length   = c(1000, 5000, 10000),
  geno_missing_imputation = "fixed",
  plink_path              = "plink2",
  M                       = 5L,
  seed                    = NULL,
  genome_build            = "hg19",
  ld_block_file           = NULL,
  sample_uncorrelated     = TRUE,
  grm_file                = NULL,
  grm_id_file             = NULL,
  relatedness_cutoff      = 0.125,
  n_markers_grm           = 1000L,
  fdr                     = 0.1,
  chromosomes             = 1:22,
  batch_size              = 20L,
  read_mid_exist          = TRUE,
  pipeline_stage          = "full",
  save_knockoff           = NULL,
  knockoff_dir            = NULL,
  temp_dir                = NULL
) {

  # ---------------------------------------------------------------------------
  # 1.  Input validation
  # ---------------------------------------------------------------------------

  stopifnot(
    "outdir must be a single non-empty string"        = is.character(outdir)     && length(outdir)     == 1L && nzchar(outdir),
    "geno_file must be a single non-empty string"     = is.character(geno_file)  && length(geno_file)  == 1L && nzchar(geno_file),
    "M must be a positive integer"                    = is.numeric(M)            && length(M)          == 1L && is.finite(M) && M >= 1L && M == as.integer(M),
    "fdr must be numeric in (0, 1)"                   = is.numeric(fdr)          && length(fdr)        == 1L && fdr > 0 && fdr < 1,
    "user_cores must be a positive integer"           = is.numeric(user_cores)   && length(user_cores) == 1L && is.finite(user_cores) && user_cores >= 1L && user_cores == as.integer(user_cores),
    "batch_size must be a positive integer"           = is.numeric(batch_size)   && length(batch_size) == 1L && is.finite(batch_size) && batch_size >= 1L && batch_size == as.integer(batch_size),
    "n_markers_grm must be a positive integer"         = is.numeric(n_markers_grm) && length(n_markers_grm) == 1L && is.finite(n_markers_grm) && n_markers_grm >= 1L && n_markers_grm == as.integer(n_markers_grm),
    "relatedness_cutoff must be in (0, 1]"              = is.numeric(relatedness_cutoff) && length(relatedness_cutoff) == 1L && is.finite(relatedness_cutoff) && relatedness_cutoff > 0 && relatedness_cutoff <= 1,
    "sample_uncorrelated must be TRUE or FALSE"        = is.logical(sample_uncorrelated) && length(sample_uncorrelated) == 1L && !is.na(sample_uncorrelated),
    "read_mid_exist must be TRUE or FALSE"             = is.logical(read_mid_exist) && length(read_mid_exist) == 1L && !is.na(read_mid_exist)
  )

  if (!test_type %in% c("Single_Window", "Gene_Centric"))
    stop("'test_type' must be \"Single_Window\" or \"Gene_Centric\".")
  if (identical(test_type, "Single_Window") && !isTRUE(sample_uncorrelated))
    stop("'Single_Window' currently supports only the standard GLM path (sample_uncorrelated = TRUE), which does not explicitly adjust for sample relatedness. The BIGKnock/SAIGE GLMM path is available only for 'Gene_Centric'.")
  if (!genome_build %in% c("hg19", "hg38"))
    stop("'genome_build' must be \"hg19\" or \"hg38\".")
  if (!pipeline_stage %in% c("full", "stage1_knockoff", "stage2_analysis"))
    stop("'pipeline_stage' must be one of \"full\", \"stage1_knockoff\", \"stage2_analysis\".")
  if (!is.null(covar_cols)     && !is.character(covar_cols))     stop("'covar_cols' must be a character vector or NULL.")
  if (!is.null(cat_covar_cols) && !is.character(cat_covar_cols)) stop("'cat_covar_cols' must be a character vector or NULL.")
  if (length(intersect(covar_cols, cat_covar_cols)) > 0L)
    stop("A covariate cannot appear in both 'covar_cols' and 'cat_covar_cols'.")
  if (!is.null(pheno_id) && (!is.character(pheno_id) || length(pheno_id) != 1L))
    stop("'pheno_id' must be a single string or NULL.")
  if (!is.null(save_knockoff) &&
      (!is.logical(save_knockoff) || length(save_knockoff) != 1L || is.na(save_knockoff)))
    stop("'save_knockoff' must be TRUE, FALSE, or NULL.")
  if (!is.null(temp_dir) &&
      (!is.character(temp_dir) || length(temp_dir) != 1L ||
       is.na(temp_dir) || !nzchar(temp_dir)))
    stop("'temp_dir' must be NULL or one non-empty path.")
  if (!is.character(geno_missing_imputation) ||
      length(geno_missing_imputation) != 1L ||
      !geno_missing_imputation %in% c("fixed", "random", "bestguess"))
    stop("'geno_missing_imputation' must be one of 'fixed', 'random', or 'bestguess'.")
  if (identical(test_type, "Gene_Centric") &&
      !identical(geno_missing_imputation, "fixed"))
    stop("'Gene_Centric' currently supports geno_missing_imputation = 'fixed' only.")
  if (xor(is.null(grm_file), is.null(grm_id_file)))
    stop("'grm_file' and 'grm_id_file' must be supplied together.")
  if (!is.null(seed)) .derive_unit_seed(seed, "validation")
  if (identical(test_type, "Single_Window") &&
      identical(geno_missing_imputation, "random") && is.null(seed))
    stop("'seed' is required when geno_missing_imputation = 'random'.")
  if (!is.null(ld_block_file)) {
    if (!is.character(ld_block_file) || length(ld_block_file) != 1L ||
        !nzchar(ld_block_file) || !file.exists(ld_block_file))
      stop("'ld_block_file' must be NULL or the path to an existing file.")
  }

  chr_numeric <- suppressWarnings(as.integer(chromosomes))
  if (any(is.na(chr_numeric)))
    stop("'chromosomes' contains non-integer values: ", paste(chromosomes[is.na(chr_numeric)], collapse = ", "))
  chr_vector  <- intersect(chr_numeric, 1L:22L)
  if (length(chr_vector) == 0L) stop("No valid autosomes (1-22) in 'chromosomes'.")

  if (is.null(pheno_file) || !is.character(pheno_file) ||
      length(pheno_file) != 1L || !nzchar(pheno_file))
    stop("'pheno_file' is required for every pipeline stage.")
  if (!file.exists(pheno_file))
    stop("Phenotype file not found: ", pheno_file)
  if (is.null(phenotype) || !is.character(phenotype) ||
      length(phenotype) == 0L || !all(nzchar(phenotype)))
    stop("'phenotype' must be a non-empty character vector for every pipeline stage.")
  plink_fam <- paste0(geno_file, ".fam")
  if (!file.exists(plink_fam))
    stop("PLINK .fam not found: ", plink_fam, "\nCheck that 'geno_file' is the correct prefix.")

  # ---------------------------------------------------------------------------
  # 2.  Resolve save/load flags
  # ---------------------------------------------------------------------------

  multi_pheno <- !is.null(phenotype) && length(phenotype) > 1L

  # save_knockoff controls retention after the run, not whether temporary
  # on-disk persistence may be needed internally during the run. Keep all
  # stage/save/load decisions in one pure helper so the state machine can be
  # tested without running PLINK or the association methods.
  ko_plan <- .resolve_knockoff_persistence(
    pipeline_stage = pipeline_stage,
    save_knockoff  = save_knockoff,
    multi_pheno    = multi_pheno
  )
  save_knockoff <- ko_plan$save_knockoff

  if (pipeline_stage == "stage2_analysis" && is.null(knockoff_dir))
    stop("pipeline_stage = 'stage2_analysis' requires 'knockoff_dir' to be specified.")

  # internal_persist: whether the first pass must write knockoffs to disk.
  internal_persist <- ko_plan$write_first_pass

  # Resolve knockoff directory
  if (is.null(knockoff_dir)) knockoff_dir <- file.path(outdir, "knockoffs")
  knockoff_sample_file <- file.path(knockoff_dir, "knockoff_sample_list.txt")

  # ---------------------------------------------------------------------------
  # 3.  Create directories
  # ---------------------------------------------------------------------------

  if (!dir.exists(outdir)) { message("Creating output directory: ", outdir); dir.create(outdir, recursive = TRUE) }
  if (is.null(temp_dir)) {
    slurm_tmp <- Sys.getenv("SLURM_TMPDIR", unset = "")
    temp_dir <- if (nzchar(slurm_tmp) && dir.exists(slurm_tmp) &&
                    file.access(slurm_tmp, mode = 2L) == 0L) {
      slurm_tmp
    } else {
      tempdir()
    }
  }
  if (!dir.exists(temp_dir) &&
      !dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE))
    stop("Unable to create temporary directory: ", temp_dir)
  temp_dir <- normalizePath(temp_dir, winslash = "/", mustWork = TRUE)
  if (file.access(temp_dir, mode = 2L) != 0L)
    stop("Temporary directory is not writable: ", temp_dir)
  knockoff_dir_created <- .prepare_knockoff_directory(
    knockoff_dir,
    create           = ko_plan$write_first_pass,
    require_existing = ko_plan$load_first_pass
  )
  if (knockoff_dir_created) message("Creating knockoff directory: ", knockoff_dir)
  if (user_cores > 1L) {
    thread_vars <- c(
      "MKL_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
      "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"
    )
    previous_threads <- Sys.getenv(thread_vars, unset = NA_character_)
    on.exit({
      present <- !is.na(previous_threads)
      if (any(present))
        do.call(Sys.setenv, as.list(stats::setNames(
          previous_threads[present], thread_vars[present]
        )))
      if (any(!present)) Sys.unsetenv(thread_vars[!present])
    }, add = TRUE)
    do.call(Sys.setenv, as.list(stats::setNames(
      rep("1", length(thread_vars)), thread_vars
    )))
  }

  # ---------------------------------------------------------------------------
  # 4.  Determine sample set
  #     Every stage uses the same phenotype/covariate complete-case rule and
  #     aligns the retained samples to the PLINK .fam order.
  # ---------------------------------------------------------------------------

  sample_data <- .prepare_analysis_samples(
    pheno_file = pheno_file,
    phenotype = phenotype,
    pheno_id = pheno_id,
    covar_cols = covar_cols,
    cat_covar_cols = cat_covar_cols,
    plink_fam = plink_fam
  )
  pheno <- sample_data$pheno
  Gsub.id <- sample_data$sample_ids
  plink_keep_fam <- sample_data$plink_keep_fam
  all_covar_cols <- sample_data$all_covar_cols
  all_fam_samples_retained <- sample_data$all_fam_samples_retained
  rm(sample_data); gc()

  # ---------------------------------------------------------------------------
  # 6.  Stage 2: validate the exact saved sample set and restore its row order
  # ---------------------------------------------------------------------------

  # The 'canonical' IDs used for knockoff row-ordering are what we call
  # knockoff_sample_ids.  This is written during stage1 / multi-pheno pass 1
  # and read back during stage2 / multi-pheno pass >=2.

  if (pipeline_stage == "stage2_analysis" || internal_persist) {
    # A pre-existing directory is a checkpoint unit.  Validate its sample list
    # before changing anything, then preserve its canonical row order.  This
    # keeps a failed or interrupted restart from corrupting reusable knockoffs.
    sample_state <- .reconcile_knockoff_sample_file(
      current_ids = Gsub.id,
      sample_file = knockoff_sample_file,
      create       = internal_persist
    )
    current_order <- sample_state$order
    plink_keep_fam <- plink_keep_fam[current_order]
    if (!is.null(pheno)) pheno <- pheno[current_order]
    Gsub.id <- sample_state$sample_ids
    knockoff_sample_ids <- sample_state$sample_ids

    if (sample_state$created) {
      message("Knockoff sample list written to: ", knockoff_sample_file)
    } else {
      message(length(Gsub.id),
              " sample(s) matched the saved knockoff sample list.")
    }
  } else {
    knockoff_sample_ids <- Gsub.id
  }

  # Write one PLINK keep file for every genotype export. PLINK does not promise
  # to follow keep-file order, so each .raw file is also reordered explicitly
  # by IID in .prepare_raw_genotypes().
  if (isTRUE(all_fam_samples_retained)) {
    plink_keep_file <- NULL
  } else {
    plink_keep_file <- tempfile("KnockoffPipeline_keep_", fileext = ".fam")
    data.table::fwrite(plink_keep_fam, plink_keep_file, sep = "\t",
                       col.names = FALSE)
    on.exit(unlink(plink_keep_file, force = TRUE), add = TRUE)
  }

  # Detect the additive-export syntax once.  The block/batch workers reuse the
  # resolved switch instead of spawning one `plink --version` process per
  # export.  PLINK itself stays single-threaded because parallelism is managed
  # by the outer R workers.
  plink_export_switch <- .plink_additive_export_switch(plink_path)
  plink_threads <- 1L

  # ---------------------------------------------------------------------------
  # 7.  Stage 1: knockoff generation only (no null model, no tests)
  # ---------------------------------------------------------------------------

  if (pipeline_stage == "stage1_knockoff") {
    message("\n======= Stage 1: Knockoff Generation =======")
    rm(pheno, all_covar_cols)
    gc()
    .run_knockoff_generation(
      test_type               = test_type,
      geno_file               = geno_file,
      Gsub.id                 = Gsub.id,
      knockoff_dir            = knockoff_dir,
      chr_vector              = chr_vector,
      M                       = M,
      seed                    = seed,
      genome_build            = genome_build,
      sliding_window_length   = sliding_window_length,
      geno_missing_imputation = geno_missing_imputation,
      plink_path              = plink_path,
      batch_size              = batch_size,
      sample_uncorrelated     = sample_uncorrelated,
      user_cores              = user_cores,
      read_mid_exist          = read_mid_exist,
      plink_keep_file         = plink_keep_file,
      ld_block_file           = ld_block_file,
      export_switch           = plink_export_switch,
      plink_threads           = plink_threads,
      temp_dir                = temp_dir
    )
    message("\nStage 1 complete.")
    message("  Knockoffs : ", knockoff_dir)
    message("  Sample list: ", knockoff_sample_file)
    return(invisible(TRUE))
  }

  # Record compact identities of the external inputs once.  These are used to
  # prevent a restart from silently mixing intermediate results produced from
  # another phenotype file, PLINK dataset, reference, or model configuration.
  checkpoint_sources <- list(
    genotype = .plink_dataset_identity(geno_file),
    phenotype = .reference_file_id(pheno_file),
    reference = .analysis_reference_identity(
      test_type = test_type, genome_build = genome_build,
      chromosomes = chr_vector, ld_block_file = ld_block_file
    ),
    grm = .optional_file_pair_identity(grm_file, grm_id_file)
  )

  # ===========================================================================
  # 8.  Stages "full" / "stage2_analysis": loop over phenotype(s)
  # ===========================================================================

  # For multi-phenotype "full" runs:
  #   pass 1 (first phenotype) -> save_knockoff = TRUE (generate + save)
  #   pass >=2 (later phenotypes) -> load_knockoff = TRUE (load saved)
  # For "stage2_analysis": always load.

  knockoffs_ready <- ko_plan$load_first_pass

  for (pheno_name in phenotype) {

    message("\n======= Phenotype: \"", pheno_name, "\" =======")

    # Per-phenotype subdirectory (only created when multiple phenotypes)
    p_outdir  <- if (multi_pheno) file.path(outdir, pheno_name) else outdir
    p_mid_dir <- file.path(p_outdir, "mid")
    if (!dir.exists(p_outdir)) dir.create(p_outdir, recursive = TRUE)

    # ---- Outcome type -------------------------------------------------------
    pv       <- unique(pheno[[pheno_name]])
    is_bin   <- length(pv) == 2L && all(sort(as.numeric(pv)) == c(0, 1))
    out_type <- if (is_bin) "D" else "C"
    message("Outcome type: ", if (is_bin) "binary (D)" else "continuous (C)")

    checkpoint_context <- .make_analysis_checkpoint_context(
      test_type = test_type, phenotype = pheno_name,
      sample_ids = Gsub.id, pheno_id = pheno_id,
      covar_cols = covar_cols, cat_covar_cols = cat_covar_cols,
      out_type = out_type, M = M, seed = seed,
      genome_build = genome_build,
      sliding_window_length = sliding_window_length,
      geno_missing_imputation = geno_missing_imputation,
      sample_uncorrelated = sample_uncorrelated,
      relatedness_cutoff = relatedness_cutoff,
      n_markers_grm = n_markers_grm, fdr = fdr,
      sources = checkpoint_sources
    )
    .prepare_analysis_checkpoint(
      mid_dir = p_mid_dir, context = checkpoint_context,
      read_mid_exist = read_mid_exist
    )

    # ---- Covariate matrix ---------------------------------------------------
    covar_matrix <- .build_glm_covariates(
      pheno, covar_cols = covar_cols,
      cat_covar_cols = cat_covar_cols
    )

    # ---- Fit null model -----------------------------------------------------
    message("Fitting null model (sample_uncorrelated = ", sample_uncorrelated, ") ...")

    if (sample_uncorrelated) {
      nm_args <- list(Y = pheno[[pheno_name]], X = covar_matrix,
                      id = Gsub.id, out_type = out_type)
      nullobj <- do.call(Fit_null_model, nm_args)
    } else {
      saige_input_dir <- file.path(p_outdir, "saige_input")
      if (!dir.exists(saige_input_dir)) dir.create(saige_input_dir, recursive = TRUE)

      pheno_saige_file <- file.path(saige_input_dir, paste0("pheno_", pheno_name, ".tsv"))
      pheno_saige <- data.table::copy(pheno)
      saige_sample_id_col <- pheno_id
      if (is.null(saige_sample_id_col)) {
        saige_sample_id_col <- ".KnockoffPipeline_IID"
        pheno_saige[, (saige_sample_id_col) := Gsub.id]
      }
      data.table::fwrite(pheno_saige, pheno_saige_file, sep = "\t")

      nullobj <- .with_local_seed(
        .derive_unit_seed(seed, "SAIGE_null", pheno_name),
        function() Fit_null_model_GLMM(
          geno_file, pheno_saige_file, pheno_name, plink_path,
          outcome_type       = out_type,
          sample_id_col      = saige_sample_id_col,
          covar_cols         = covar_cols,
          cat_covar_cols     = cat_covar_cols,
          output_prefix      = file.path(p_outdir, "saige_output"),
          sparse_grm_file    = grm_file,
          sparse_grm_id_file = grm_id_file,
          n_threads          = user_cores,
          num_random_marker_for_sparse_kin = as.integer(n_markers_grm),
          relatedness_cutoff = relatedness_cutoff,
          random_seed        = .derive_unit_seed(
            seed, "SAIGE_PLINK_thinning", pheno_name
          )
        )
      )
    }

    # ---- Knockoff flags for this pass ---------------------------------------
    pass_flags <- .knockoff_pass_flags(ko_plan, knockoffs_ready)
    ko_save <- pass_flags$save_knockoff
    ko_load <- pass_flags$load_knockoff

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
        ld_block_file           = ld_block_file,
        plink_keep_file         = plink_keep_file,
        sliding_window_length   = sliding_window_length,
        geno_missing_imputation = geno_missing_imputation,
        plink_path              = plink_path,
        M                       = M,
        seed                    = seed,
        user_cores              = user_cores,
        read_mid_exist          = read_mid_exist,
        fdr                     = fdr,
        save_knockoff           = ko_save,
        load_knockoff           = ko_load,
        knockoff_dir            = knockoff_dir,
        knockoff_sample_ids     = knockoff_sample_ids,
        export_switch           = plink_export_switch,
        plink_threads           = plink_threads,
        temp_dir                = temp_dir
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
        plink_keep_file         = plink_keep_file,
        sliding_window_length   = sliding_window_length,
        plink_path              = plink_path,
        M                       = M,
        seed                    = seed,
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
        ratio                   = if (!sample_uncorrelated) nullobj$ratio       else NULL,
        export_switch           = plink_export_switch,
        plink_threads           = plink_threads,
        temp_dir                = temp_dir
      )
    }

    if (ko_plan$reuse_later_passes)
      knockoffs_ready <- TRUE   # knockoffs now exist for subsequent phenotypes
    gc()
  }

  removed_knockoff_dir <- .remove_owned_knockoff_directory(
    knockoff_dir,
    cleanup        = ko_plan$cleanup_after_run,
    created_by_run = knockoff_dir_created
  )
  if (removed_knockoff_dir) {
    message("Knockoff directory removed (save_knockoff = FALSE).")
  } else if (ko_plan$cleanup_after_run && dir.exists(knockoff_dir)) {
    message("Pre-existing knockoff directory retained; it is not owned by this run: ",
            knockoff_dir)
  }

  invisible(TRUE)
}


# =============================================================================
# Internal: persistence state and directory ownership helpers
# =============================================================================

#' Resolve knockoff persistence flags without touching the filesystem
#' @keywords internal
.resolve_knockoff_persistence <- function(pipeline_stage, save_knockoff,
                                           multi_pheno) {
  if (!pipeline_stage %in% c("full", "stage1_knockoff", "stage2_analysis"))
    stop("Unknown pipeline stage: ", pipeline_stage)
  if (!is.logical(multi_pheno) || length(multi_pheno) != 1L || is.na(multi_pheno))
    stop("'multi_pheno' must be one non-missing logical value.")
  if (!is.null(save_knockoff) &&
      (!is.logical(save_knockoff) || length(save_knockoff) != 1L || is.na(save_knockoff)))
    stop("'save_knockoff' must be TRUE, FALSE, or NULL.")

  if (is.null(save_knockoff))
    save_knockoff <- identical(pipeline_stage, "stage1_knockoff")
  if (identical(pipeline_stage, "stage1_knockoff") &&
      !isTRUE(save_knockoff))
    stop("pipeline_stage = 'stage1_knockoff' requires save_knockoff = TRUE (or NULL).")

  write_first_pass <- identical(pipeline_stage, "stage1_knockoff") ||
    (identical(pipeline_stage, "full") &&
       (isTRUE(save_knockoff) || isTRUE(multi_pheno)))

  list(
    save_knockoff      = isTRUE(save_knockoff),
    write_first_pass   = write_first_pass,
    load_first_pass    = identical(pipeline_stage, "stage2_analysis"),
    reuse_later_passes = identical(pipeline_stage, "full") &&
      isTRUE(multi_pheno),
    cleanup_after_run  = identical(pipeline_stage, "full") &&
      isTRUE(multi_pheno) && !isTRUE(save_knockoff)
  )
}


#' Resolve save/load flags for one phenotype pass
#' @keywords internal
.knockoff_pass_flags <- function(plan, knockoffs_ready) {
  if (!is.logical(knockoffs_ready) || length(knockoffs_ready) != 1L ||
      is.na(knockoffs_ready))
    stop("'knockoffs_ready' must be one non-missing logical value.")
  list(
    save_knockoff = isTRUE(plan$write_first_pass) && !knockoffs_ready,
    load_knockoff = isTRUE(plan$load_first_pass) ||
      (isTRUE(plan$reuse_later_passes) && knockoffs_ready)
  )
}


#' Create or validate a knockoff directory and report run ownership
#'
#' @return \code{TRUE} only when this call created \code{path}.
#' @keywords internal
.prepare_knockoff_directory <- function(path, create = FALSE,
                                         require_existing = FALSE) {
  existed <- dir.exists(path)
  if (isTRUE(require_existing) && !existed)
    stop("Knockoff directory not found: ", path,
         "\nRun pipeline_stage = 'stage1_knockoff' first.")

  if (isTRUE(create) && !existed) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(ok) && !dir.exists(path))
      stop("Unable to create knockoff directory: ", path)
    return(TRUE)
  }
  FALSE
}


#' Remove a temporary knockoff directory only when this run owns it
#' @keywords internal
.remove_owned_knockoff_directory <- function(path, cleanup,
                                              created_by_run) {
  if (!isTRUE(cleanup) || !isTRUE(created_by_run) || !dir.exists(path))
    return(FALSE)
  unlink(path, recursive = TRUE, force = TRUE)
  !dir.exists(path)
}


# =============================================================================
# Internal: restart-checkpoint compatibility helpers
# =============================================================================

#' Compact file identity without hashing a potentially very large binary file
#' @keywords internal
.file_stat_identity <- function(path) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) stop("Required input file not found: ", path)
  info <- file.info(path)
  list(
    path = normalizePath(path, mustWork = TRUE),
    size = as.double(info$size[[1L]]),
    mtime = as.double(info$mtime[[1L]])
  )
}


#' Identify the three files forming a PLINK bed/bim/fam dataset
#' @keywords internal
.plink_dataset_identity <- function(geno_file) {
  bed <- paste0(geno_file, ".bed")
  bim <- paste0(geno_file, ".bim")
  fam <- paste0(geno_file, ".fam")
  list(
    # The BED may be hundreds of gigabytes, so use path/size/mtime rather than
    # hashing the entire file on every restart.
    bed = .file_stat_identity(bed),
    bim = .reference_file_id(bim),
    fam = .reference_file_id(fam)
  )
}


#' Identity for an optional sparse-GRM/data-ID pair
#' @keywords internal
.optional_file_pair_identity <- function(first, second) {
  if (is.null(first) && is.null(second)) return(NULL)
  list(first = .file_stat_identity(first), second = .reference_file_id(second))
}


#' Build-specific reference identity used by an analysis checkpoint
#' @keywords internal
.analysis_reference_identity <- function(test_type, genome_build,
                                         chromosomes, ld_block_file = NULL) {
  if (identical(test_type, "Single_Window")) {
    return(.reference_file_id(
      .resolve_ld_block_file(genome_build, ld_block_file)
    ))
  }

  enhancer_files <- unlist(lapply(chromosomes, function(chr) {
    unname(.enhancer_reference_paths(genome_build, chr))
  }), use.names = FALSE)
  list(
    gene_annotation = .reference_file_id(
      .gene_annotation_path(genome_build)
    ),
    # Enhancer maps affect gene-level scores but can be moderately large. File
    # path/size/mtime detects normal replacement without costly repeated hashes.
    enhancer_maps = lapply(enhancer_files, .file_stat_identity)
  )
}


#' Create a normalized context for intermediate-score restart files
#' @keywords internal
.make_analysis_checkpoint_context <- function(
  test_type, phenotype, sample_ids, pheno_id, covar_cols, cat_covar_cols,
  out_type, M, seed, genome_build, sliding_window_length,
  geno_missing_imputation, sample_uncorrelated, relatedness_cutoff,
  n_markers_grm, fdr, sources
) {
  list(
    # Version 2 invalidates checkpoints written before the corrected Single
    # MAC alignment and gene knockoff skip-index logic.
    schema_version = 2L,
    test_type = as.character(test_type),
    phenotype = as.character(phenotype),
    sample_ids = as.character(sample_ids),
    pheno_id = if (is.null(pheno_id)) NA_character_ else as.character(pheno_id),
    covar_cols = sort(as.character(covar_cols)),
    cat_covar_cols = sort(as.character(cat_covar_cols)),
    out_type = as.character(out_type),
    M = as.integer(M),
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
    genome_build = as.character(genome_build),
    sliding_window_length = as.double(sliding_window_length),
    geno_missing_imputation = as.character(geno_missing_imputation),
    sample_uncorrelated = isTRUE(sample_uncorrelated),
    relatedness_cutoff = as.double(relatedness_cutoff),
    n_markers_grm = as.integer(n_markers_grm),
    fdr = as.double(fdr),
    sources = sources
  )
}


#' Create, validate, or deliberately replace an analysis checkpoint directory
#'
#' Existing intermediates without a manifest are rejected when resuming. Setting
#' read_mid_exist = FALSE means a true fresh run: only the dedicated mid
#' directory is removed and recreated before any results are written.
#' @keywords internal
.prepare_analysis_checkpoint <- function(mid_dir, context,
                                         read_mid_exist) {
  manifest <- file.path(mid_dir, "checkpoint_context.rds")

  if (!isTRUE(read_mid_exist)) {
    if (dir.exists(mid_dir)) unlink(mid_dir, recursive = TRUE, force = TRUE)
    if (!dir.create(mid_dir, recursive = TRUE, showWarnings = FALSE) &&
        !dir.exists(mid_dir))
      stop("Unable to create fresh intermediate directory: ", mid_dir)
    saveRDS(context, manifest)
    return(invisible("fresh"))
  }

  if (!dir.exists(mid_dir)) {
    if (!dir.create(mid_dir, recursive = TRUE, showWarnings = FALSE) &&
        !dir.exists(mid_dir))
      stop("Unable to create intermediate directory: ", mid_dir)
  }

  entries <- setdiff(list.files(mid_dir, all.files = TRUE), c(".", ".."))
  payload <- setdiff(entries, basename(manifest))
  if (!file.exists(manifest)) {
    if (length(payload) > 0L) {
      stop(
        "Existing intermediate files lack compatibility metadata in: ",
        mid_dir, ". Use a new output directory, or set read_mid_exist = FALSE ",
        "to replace this dedicated mid directory with a fresh run."
      )
    }
    saveRDS(context, manifest)
    return(invisible("new"))
  }

  saved <- tryCatch(
    readRDS(manifest),
    error = function(e) stop("Cannot read checkpoint manifest: ", manifest,
                             ". ", conditionMessage(e))
  )
  if (!identical(saved, context)) {
    fields <- union(names(saved), names(context))
    mismatched <- fields[!vapply(fields, function(field) {
      identical(saved[[field]], context[[field]])
    }, logical(1))]
    stop(
      "Intermediate checkpoint is incompatible with the current run ",
      "(mismatch: ", paste(mismatched, collapse = ", "), "): ", manifest,
      ". Use a new output directory, or set read_mid_exist = FALSE for a fresh run."
    )
  }
  invisible("resume")
}


# =============================================================================
# Internal: Stage-1 knockoff generation only (no association testing)
# Dispatches to Single_Window or Gene_Centric generation helpers.
# =============================================================================

.ordered_chr_genes <- function(genes_info, chromosome) {
  chr_genes <- genes_info[which(genes_info$chr == chromosome), ]
  chr_genes[order(chr_genes$start), ]
}


.run_knockoff_generation <- function(
  test_type, geno_file, Gsub.id, knockoff_dir, chr_vector,
  M, seed, genome_build, sliding_window_length, geno_missing_imputation,
  plink_path, batch_size, sample_uncorrelated, user_cores, read_mid_exist,
  plink_keep_file, ld_block_file, export_switch = NULL, plink_threads = 1L,
  temp_dir = NULL
) {
  # A minimal "null object" is not needed here: run_single_block /
  # run_batch_gene accept save_knockoff = TRUE without running tests.
  # We pass nullobj = NULL and the analysis branches will short-circuit.

  if (test_type == "Single_Window") {
    block_file <- .resolve_ld_block_file(genome_build, ld_block_file)
    reference_id <- .reference_file_id(block_file)

    blocks     <- data.table::fread(block_file)
    unique_chr <- intersect(sort(unique(blocks$chr)), chr_vector)
    if (length(unique_chr) == 0L) stop("No chromosomes after intersecting block file with requested chromosomes.")

    for (c in unique_chr) {
      message("--- chr ", c, " ---")
      bim_chr <- .read_plink_bim_chr(geno_file, c, plink_path)
      chr_ko_dir <- file.path(knockoff_dir, paste0("chr", c))
      .prepare_knockoff_directory(chr_ko_dir, create = TRUE)
      block_chr  <- blocks[blocks$chr == c]
      n_blocks   <- nrow(block_chr)
      message("  chr ", c, ": ", n_blocks, " blocks")

      generated <- parallel::mclapply(seq_len(n_blocks), function(kk) {
        message("  Block ", kk, " / ", n_blocks)
        ko_file <- .ko_file_single(chr_ko_dir, kk)
        reuse_existing <- isTRUE(read_mid_exist) && file.exists(ko_file)
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
          bim_metadata            = bim_chr,
          genome_build            = genome_build,
          reference_id            = reference_id,
          plink_keep_file         = plink_keep_file,
          export_switch           = export_switch,
          plink_threads           = plink_threads,
          temp_dir                = temp_dir,
          knockoff_seed           = .derive_unit_seed(
            seed, "Single_Window", c, kk
          ),
          save_knockoff           = !reuse_existing,
          load_knockoff           = reuse_existing,
          knockoff_file           = ko_file,
          knockoff_sample_ids     = Gsub.id,
          stage1_only             = TRUE    # skip association test
        )
      }, mc.cores = user_cores)
      failed <- vapply(generated, inherits, logical(1), what = "try-error")
      if (any(failed))
        stop("Stage-1 Single_Window knockoff generation failed on chr", c,
             ": ", paste(as.character(generated[failed]), collapse = "; "))
    }

  } else {
    # Gene_Centric
    gene_file <- .gene_annotation_path(genome_build)
    if (!file.exists(gene_file)) stop("Gene annotation file not found: ", gene_file)
    reference_id <- .reference_file_id(gene_file)

    genes_info <- data.table::fread(gene_file)
    genes_info$chr <- as.numeric(gsub("[^0-9]", "", genes_info$chr))
    genes_info      <- genes_info[!is.na(chr)]
    unique_chr      <- intersect(sort(unique(genes_info$chr)), chr_vector)
    if (length(unique_chr) == 0L) stop("No chromosomes remain.")

    for (c in unique_chr) {
      message("--- chr ", c, " ---")
      bim_chr <- .read_plink_bim_chr(geno_file, c, plink_path)
      chr_ko_dir <- file.path(knockoff_dir, paste0("chr", c))
      .prepare_knockoff_directory(chr_ko_dir, create = TRUE)

      chr_genes <- .ordered_chr_genes(genes_info, c)
      enhancer_files <- .enhancer_reference_paths(genome_build, c)
      abc_df    <- data.table::fread(enhancer_files$abc)
      gh_df     <- data.table::fread(enhancer_files$gh)

      batch_index <- split(seq_len(nrow(chr_genes)), ceiling(seq_len(nrow(chr_genes)) / batch_size))
      for (b in seq_along(batch_index)) {
        message("  Batch ", b, " / ", length(batch_index))
        run_batch_gene(
          genes               = chr_genes,
          b                   = b,
          batch_index         = batch_index,
          geno.file           = geno_file,
          obj_nullmodel       = NULL,
          window_length       = sliding_window_length,
          plink_prefix        = plink_path,
          M                   = M,
          genome_build        = genome_build,
          Gsub.id             = Gsub.id,
          bim_metadata        = bim_chr,
          plink_keep_file     = plink_keep_file,
          reference_id        = reference_id,
          export_switch       = export_switch,
          plink_threads       = plink_threads,
          temp_dir            = temp_dir,
          seed                = seed,
          use_glmm            = !sample_uncorrelated,
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
# Helpers for incremental intermediate result saving
# =============================================================================

#' Ensure the per-block result directory exists for a chromosome
#' @keywords internal
.ensure_block_dir <- function(mid_dir, chr) {
  d <- file.path(mid_dir, "blocks", paste0("chr", chr))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

#' Ensure the per-batch result directory exists for a chromosome
#' @keywords internal
.ensure_batch_dir <- function(mid_dir, chr) {
  d <- file.path(mid_dir, "batches", paste0("chr", chr))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

#' Path to the Gene_Centric progress file for a chromosome
#' @keywords internal
.progress_file <- function(mid_dir, chr) {
  file.path(mid_dir, paste0("progress_chr", chr, ".txt"))
}

#' Get block indices that already have per-block result files
#' @param mid_dir Character. Mid-results directory.
#' @param chr     Integer. Chromosome number.
#' @param prefix  Character. "Single" or "Window".
#' @return Integer vector of completed block indices (possibly empty).
#' @keywords internal
.get_completed_blocks <- function(mid_dir, chr) {
  d <- file.path(mid_dir, "blocks", paste0("chr", chr))
  if (!dir.exists(d)) return(integer(0))
  files <- list.files(d, pattern = "^done_(\\d+)\\.txt$")
  if (length(files) == 0) return(integer(0))
  as.integer(gsub("^done_|\\.txt$", "", files))
}

#' Write per-block results to disk (safe for use inside mclapply workers)
#' @keywords internal
.write_block_result <- function(mid_dir, chr, kk, single_df, window_df) {
  d <- .ensure_block_dir(mid_dir, chr)
  has_single <- !is.null(single_df) && nrow(single_df) > 0L
  has_window <- !is.null(window_df) && nrow(window_df) > 0L
  if (has_single) {
    data.table::fwrite(single_df,
      file.path(d, sprintf("Single_block_%04d.txt", kk)), sep = "\t")
  }
  if (has_window) {
    data.table::fwrite(window_df,
      file.path(d, sprintf("Window_block_%04d.txt", kk)), sep = "\t")
  }
  # Always write progress marker (even for empty blocks)
  .mark_block_done(mid_dir, chr, kk)
  invisible(structure(
    list(
      chr = as.integer(chr), block = as.integer(kk),
      empty = !has_single && !has_window
    ),
    class = "knockoff_pipeline_block_status"
  ))
}

#' Write a progress marker for a completed block
#' @keywords internal
.mark_block_done <- function(mid_dir, chr, kk) {
  d <- .ensure_block_dir(mid_dir, chr)
  writeLines(as.character(kk), file.path(d, sprintf("done_%04d.txt", kk)))
}

#' Merge per-block files into chromosome-level files (backward compat)
#' @keywords internal
.merge_per_block_files <- function(mid_dir, chr) {
  d <- file.path(mid_dir, "blocks", paste0("chr", chr))
  single_files <- list.files(d, pattern = "^Single_block_\\d+\\.txt$",
                             full.names = TRUE)
  window_files <- list.files(d, pattern = "^Window_block_\\d+\\.txt$",
                             full.names = TRUE)
  single_chr <- NULL
  window_chr <- NULL

  if (length(single_files) > 0L) {
    single_chr <- data.table::rbindlist(
      lapply(single_files, data.table::fread), fill = TRUE)
    data.table::fwrite(single_chr,
      file.path(mid_dir, paste0("Single_mid_results_chr", chr, ".txt")),
      sep = "\t")
    message("  Merged ", length(single_files),
            " per-block Single files for chr ", chr)
  }
  if (length(window_files) > 0L) {
    window_chr <- data.table::rbindlist(
      lapply(window_files, data.table::fread), fill = TRUE)
    data.table::fwrite(window_chr,
      file.path(mid_dir, paste0("Window_mid_results_chr", chr, ".txt")),
      sep = "\t")
    message("  Merged ", length(window_files),
            " per-block Window files for chr ", chr)
  }
  invisible(list(single = single_chr, window = window_chr))
}

#' Read completed gene IDs from the progress file
#' @return Character vector of gene IDs (empty if no progress file).
#' @keywords internal
.read_progress_genes <- function(mid_dir, chr) {
  f <- .progress_file(mid_dir, chr)
  if (!file.exists(f)) return(character(0))
  g <- readLines(f, warn = FALSE)
  g[nzchar(g)]
}

#' Append gene IDs to the progress file (one per line)
#' @keywords internal
.write_progress_genes <- function(mid_dir, chr, gene_ids) {
  f <- .progress_file(mid_dir, chr)
  cat(paste0(gene_ids, collapse = "\n"), "\n",
      file = f, append = TRUE, sep = "")
}

#' Write per-batch result file for Gene_Centric
#' @keywords internal
.write_batch_result <- function(mid_dir, chr, run_tag, b, result_df) {
  d <- .ensure_batch_dir(mid_dir, chr)
  out_file <- file.path(d,
    sprintf("GeneCentric_batch_%s_b%04d.txt", run_tag, b))
  data.table::fwrite(result_df, out_file, sep = "\t")
}

#' Merge all per-batch files into the chromosome-level file
#'
#' De-duplicates by gene_id to handle batch_size changes across runs.
#' @keywords internal
.merge_batch_files <- function(mid_dir, chr) {
  d <- file.path(mid_dir, "batches", paste0("chr", chr))
  if (!dir.exists(d)) {
    warning("No batch directory found for chr ", chr, " -- cannot merge.")
    return(invisible(NULL))
  }
  batch_files <- list.files(d, pattern = "^GeneCentric_batch_.*_b\\d+\\.txt$",
                            full.names = TRUE)
  if (length(batch_files) == 0L) {
    warning("No batch files found for chr ", chr)
    return(invisible(NULL))
  }

  result_chr <- data.table::rbindlist(
    lapply(batch_files, data.table::fread), fill = TRUE)
  if ("gene_id" %in% names(result_chr) && any(duplicated(result_chr$gene_id))) {
    n_dup <- sum(duplicated(result_chr$gene_id))
    result_chr <- result_chr[!duplicated(result_chr$gene_id), ]
    message("  Removed ", n_dup, " duplicate gene entries during batch merge")
  }
  data.table::fwrite(result_chr,
    file.path(mid_dir, paste0("GeneCentric_mid_results_chr", chr, ".txt")),
    sep = "\t")
  message("  Merged ", length(batch_files), " per-batch files for chr ", chr)
  invisible(result_chr)
}


# =============================================================================
# Internal: Single_Window analysis for one phenotype
# =============================================================================

.run_single_window <- function(
  outdir, mid_dir, geno_file, nullobj, Gsub.id, chr_vector, genome_build,
  ld_block_file, plink_keep_file,
  sliding_window_length, geno_missing_imputation, plink_path, M, seed,
  user_cores, read_mid_exist, fdr,
  save_knockoff, load_knockoff, knockoff_dir, knockoff_sample_ids,
  export_switch = NULL, plink_threads = 1L, temp_dir = NULL
) {
  block_file <- .resolve_ld_block_file(genome_build, ld_block_file)
  reference_id <- .reference_file_id(block_file)

  blocks     <- data.table::fread(block_file)
  unique_chr <- intersect(sort(unique(blocks$chr)), chr_vector)
  if (length(unique_chr) == 0L) stop("No chromosomes remain after intersecting block file.")
  last_chr <- utils::tail(unique_chr, 1L)
  last_merge <- list(single = NULL, window = NULL)

  for (c in unique_chr) {
    message("--- chr ", c, " (Single_Window) ---")
    bim_chr <- .read_plink_bim_chr(geno_file, c, plink_path)
    single_mid_file <- file.path(mid_dir, paste0("Single_mid_results_chr", c, ".txt"))
    window_mid_file <- file.path(mid_dir, paste0("Window_mid_results_chr", c, ".txt"))

    # Backward compat: skip if chromosome-level files already exist
    if (read_mid_exist && file.exists(single_mid_file) && file.exists(window_mid_file)) {
      message("  Existing intermediate files found -- skipping chr ", c); next
    }

    block_chr  <- blocks[blocks$chr == c]
    n_blocks   <- nrow(block_chr)
    chr_ko_dir <- file.path(knockoff_dir, paste0("chr", c))
    .prepare_knockoff_directory(
      chr_ko_dir,
      create           = save_knockoff,
      require_existing = load_knockoff
    )
    message("  chr ", c, ": ", n_blocks, " blocks")

    # ---- Determine pending blocks from per-block progress -----------------
    if (read_mid_exist) {
      completed_blocks <- .get_completed_blocks(mid_dir, c)
    } else {
      completed_blocks <- integer(0)
    }
    pending_blocks <- setdiff(seq_len(n_blocks), completed_blocks)

    if (length(pending_blocks) == 0L) {
      message("  All blocks already completed for chr ", c)
      if (identical(c, last_chr)) {
        last_merge <- .merge_per_block_files(mid_dir, c)
      } else {
        .merge_per_block_files(mid_dir, c)
      }
      next
    }

    message("  Processing ", length(pending_blocks), " pending block(s) of ",
            n_blocks, " total")

    # Ensure block directory exists BEFORE parallel fork (avoids race)
    .ensure_block_dir(mid_dir, c)

    out <- parallel::mclapply(pending_blocks, function(kk) {
      message("  Block ", kk, " / ", n_blocks)
      tryCatch({
        res <- run_single_block(
          blocks                  = block_chr,
          kk                      = kk,
          geno.file               = geno_file,
          obj_nullmodel           = nullobj,
          window_length           = sliding_window_length,
          plink_prefix            = plink_path,
          impute.method           = geno_missing_imputation,
          M                       = M,
          Gsub.id                 = Gsub.id,
          bim_metadata            = bim_chr,
          genome_build            = genome_build,
          reference_id            = reference_id,
          plink_keep_file         = plink_keep_file,
          export_switch           = export_switch,
          plink_threads           = plink_threads,
          temp_dir                = temp_dir,
          knockoff_seed           = .derive_unit_seed(
            seed, "Single_Window", c, kk
          ),
          save_knockoff           = save_knockoff,
          load_knockoff           = load_knockoff,
          knockoff_file           = .ko_file_single(chr_ko_dir, kk),
          knockoff_sample_ids     = knockoff_sample_ids,
          stage1_only             = FALSE
        )
        # Write per-block results IMMEDIATELY (different file per block = safe).
        # The helper also records a done marker for a valid empty block.
        status <- .write_block_result(
          mid_dir, c, kk,
          if (is.null(res)) NULL else
            data.table::as.data.table(res$result.single),
          if (is.null(res)) NULL else
            data.table::as.data.table(res$result.window)
        )
        rm(res)
        status
      }, error = function(e) {
        structure(
          list(chr = c, block = kk, message = conditionMessage(e)),
          class = "knockoff_pipeline_block_error"
        )
      })
    }, mc.cores = user_cores)

    failed <- vapply(out, function(x) {
      inherits(x, "knockoff_pipeline_block_error") ||
        inherits(x, "try-error")
    }, logical(1))
    if (any(failed)) {
      messages <- vapply(out[failed], function(x) {
        if (inherits(x, "knockoff_pipeline_block_error")) {
          paste0("chr", x$chr, " block ", x$block, ": ", x$message)
        } else {
          as.character(x)
        }
      }, character(1))
      stop("Single_Window block failure(s): ", paste(messages, collapse = "; "))
    }
    valid_status <- vapply(
      out, inherits, logical(1), what = "knockoff_pipeline_block_status"
    )
    if (!all(valid_status))
      stop("Single_Window worker returned an invalid completion status for chr ", c, ".")
    n_empty <- sum(vapply(out, `[[`, logical(1), "empty"))
    if (n_empty > 0L)
      message("  ", n_empty, " completed block(s) contained no testable variants")

    # Merge per-block files into chromosome-level files (backward compat)
    if (identical(c, last_chr)) {
      last_merge <- .merge_per_block_files(mid_dir, c)
    } else {
      .merge_per_block_files(mid_dir, c)
    }
    rm(out); gc()
  }

  # ---- Merge and summarise -------------------------------------------------
  message("Merging intermediate results ...")

  read_mid <- function(chr, prefix) {
    cache_key <- if (identical(prefix, "Single")) "single" else "window"
    if (identical(chr, last_chr) && !is.null(last_merge[[cache_key]]))
      return(last_merge[[cache_key]])
    f <- file.path(mid_dir, paste0(prefix, "_mid_results_chr", chr, ".txt"))
    if (!file.exists(f)) { warning("Missing intermediate file: ", f); return(NULL) }
    data.table::fread(f)
  }

  collect_mid <- function(prefix) {
    parts <- Filter(
      Negate(is.null), lapply(unique_chr, read_mid, prefix = prefix)
    )
    if (length(parts) == 0L) return(data.table::data.table())
    if (length(parts) == 1L) return(parts[[1L]])
    data.table::rbindlist(parts, fill = TRUE)
  }

  result.single.all <- collect_mid("Single")
  # Once the Single table has been collected, do not retain its cached merge
  # while collecting the Window table.
  last_merge$single <- NULL
  result.window.all <- collect_mid("Window")
  rm(last_merge)

  if (nrow(result.single.all) == 0L || nrow(result.window.all) == 0L)
    stop("No results across all chromosomes. Check intermediate files in: ", mid_dir)

  # FIX: was `summary <- ...` (name clash + undefined summary_res below)
  summary_res <- KS_summary(
    as.matrix(result.window.all),
    as.matrix(result.single.all),
    M, fdr = fdr
  )

  keep_cols  <- c("chr", "start", "end", "W", "Qvalue", "W_Threshold", "indicator")
  result_all <- summary_res[, keep_cols]
  colnames(result_all) <- c("chr", "start", "end", "W statistics", "q-value", "threshold", "indicator")

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
  plink_keep_file,
  sliding_window_length, plink_path, M, seed, user_cores, read_mid_exist, fdr,
  batch_size, sample_uncorrelated,
  save_knockoff, load_knockoff, knockoff_dir, knockoff_sample_ids,
  sparseSigma, ratio, export_switch = NULL, plink_threads = 1L,
  temp_dir = NULL
) {
  gene_file <- .gene_annotation_path(genome_build)
  if (!file.exists(gene_file)) stop("Gene annotation file not found: ", gene_file)
  reference_id <- .reference_file_id(gene_file)

  genes_info <- data.table::fread(gene_file)
  genes_info$chr <- as.numeric(gsub("[^0-9]", "", genes_info$chr))
  genes_info      <- genes_info[!is.na(chr)]
  unique_chr      <- intersect(sort(unique(genes_info$chr)), chr_vector)
  if (length(unique_chr) == 0L) stop("No chromosomes remain.")
  last_chr <- utils::tail(unique_chr, 1L)
  last_merge <- NULL
  glmm_precomputed <- NULL

  for (c in unique_chr) {
    message("--- chr ", c, " (Gene_Centric) ---")
    bim_chr <- .read_plink_bim_chr(geno_file, c, plink_path)
    mid_file_chr <- file.path(mid_dir, paste0("GeneCentric_mid_results_chr", c, ".txt"))

    # Backward compat: skip if chromosome-level file already exists
    if (read_mid_exist && file.exists(mid_file_chr)) {
      message("  Existing intermediate file found -- skipping chr ", c); next
    }

    chr_genes <- .ordered_chr_genes(genes_info, c)
    chr_ko_dir <- file.path(knockoff_dir, paste0("chr", c))
    .prepare_knockoff_directory(
      chr_ko_dir,
      create           = save_knockoff,
      require_existing = load_knockoff
    )
    enhancer_files <- .enhancer_reference_paths(genome_build, c)
    abc_df     <- data.table::fread(enhancer_files$abc)
    gh_df      <- data.table::fread(enhancer_files$gh)

    # ---- Read progress file: filter out already-completed genes ----------
    if (read_mid_exist) {
      completed_genes <- .read_progress_genes(mid_dir, c)
      if (length(completed_genes) > 0L) {
        n_before  <- nrow(chr_genes)
        chr_genes <- chr_genes[!chr_genes$id %in% completed_genes, ]
        message("  ", n_before - nrow(chr_genes), " gene(s) already completed, ",
                nrow(chr_genes), " remaining")
      }
    }

    # All genes done -- just ensure chromosome-level file is in place
    if (nrow(chr_genes) == 0L) {
      message("  All genes already completed for chr ", c)
      if (identical(c, last_chr)) {
        last_merge <- .merge_batch_files(mid_dir, c)
      } else {
        .merge_batch_files(mid_dir, c)
      }
      next
    }

    # These matrices depend on the phenotype/null model, not on chromosome,
    # batch, or gene.  Build them lazily so a fully resumed run does no work.
    if (!sample_uncorrelated && is.null(glmm_precomputed)) {
      glmm_precomputed <- .bigknock_glmm_precompute(
        nullobj$result.null.model.GLMM, sparseSigma
      )
    }

    # Re-batch remaining genes
    batch_index     <- split(seq_len(nrow(chr_genes)),
                             ceiling(seq_len(nrow(chr_genes)) / batch_size))
    run_tag         <- format(Sys.time(), "%Y%m%d_%H%M%S")

    # Ensure batch directory exists
    .ensure_batch_dir(mid_dir, c)

    for (b in seq_along(batch_index)) {

      batch_res <- run_batch_gene(
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
        bim_metadata        = bim_chr,
        plink_keep_file     = plink_keep_file,
        reference_id        = reference_id,
        export_switch       = export_switch,
        plink_threads       = plink_threads,
        temp_dir            = temp_dir,
        seed                = seed,
        use_glmm            = !sample_uncorrelated,
        abc_df              = abc_df,
        gh_df               = gh_df,
        sparseSigma         = sparseSigma,
        ratio               = ratio,
        glmm_precomputed    = glmm_precomputed,
        user_cores          = user_cores,
        save_knockoff       = save_knockoff,
        load_knockoff       = load_knockoff,
        knockoff_dir        = chr_ko_dir,
        knockoff_sample_ids = knockoff_sample_ids,
        stage1_only         = FALSE,
        read_mid_exist      = read_mid_exist
      )

      # ---- Incremental save: result file + progress record -------------
      if (!is.null(batch_res) && nrow(batch_res) > 0L) {
        .write_batch_result(mid_dir, c, run_tag, b, batch_res)
      }
      # Always record gene-level progress (prevents retry of empty genes)
      batch_gene_ids <- as.character(chr_genes[batch_index[[b]], ]$id)
      .write_progress_genes(mid_dir, c, batch_gene_ids)

      gc()
    }

    # Merge per-batch files into chromosome-level file (backward compat)
    if (identical(c, last_chr)) {
      last_merge <- .merge_batch_files(mid_dir, c)
    } else {
      .merge_batch_files(mid_dir, c)
    }
  }

  # ---- Merge and summarise -------------------------------------------------
  gene_parts <- Filter(
    Negate(is.null), lapply(unique_chr, function(c) {
      if (identical(c, last_chr) && !is.null(last_merge)) return(last_merge)
      f <- file.path(mid_dir, paste0("GeneCentric_mid_results_chr", c, ".txt"))
      if (!file.exists(f)) { warning("Missing intermediate file: ", f); return(NULL) }
      data.table::fread(f)
    })
  )
  if (length(gene_parts) == 0L) {
    result.all <- data.table::data.table()
  } else if (length(gene_parts) == 1L) {
    result.all <- gene_parts[[1L]]
  } else {
    result.all <- data.table::rbindlist(gene_parts, fill = TRUE)
  }
  rm(gene_parts, last_merge)
  if (nrow(result.all) == 0L)
    stop("No gene-centric results. Check intermediate files in: ", mid_dir)

  summary_res <- GeneScan3DKnock_Summary(result.all, M = M, fdr = fdr)

  keep_cols  <- c("chr", "gene_id", "gene_start", "gene_end", "W", "Qvalue", "W_Threshold", "indicator")
  result_all <- summary_res[, keep_cols]
  data.table::setnames(result_all, old = c("gene_start","gene_end"), new = c("start","end"))
  colnames(result_all) <- c("chr", "gene_id", "start", "end", "W statistics", "q-value", "threshold", "indicator")

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

#' Resolve build-specific gene and enhancer resources
#' @keywords internal
.gene_annotation_path <- function(genome_build) {
  filename <- switch(
    genome_build,
    hg19 = "coding.genes.TSS.hg19.tsv",
    hg38 = "coding.genes_TSS.hg38.tsv",
    stop("Unsupported genome build: ", genome_build)
  )
  .extdata_path(genome_build, filename)
}

#' @keywords internal
.enhancer_reference_paths <- function(genome_build, chr) {
  filenames <- switch(
    genome_build,
    hg19 = list(
      abc = paste0("ABC_combined_chr", chr, ".csv"),
      gh = paste0("GH.data_chr", chr, ".csv")
    ),
    hg38 = list(
      abc = paste0("ABC_combined.hg38_chr", chr, ".csv"),
      gh = paste0("GH.data.hg38_chr", chr, ".csv")
    ),
    stop("Unsupported genome build: ", genome_build)
  )
  lapply(filenames, function(filename) .extdata_path(genome_build, filename))
}

#' Resolve the bundled or user-supplied LD-block resource
#' @keywords internal
.resolve_ld_block_file <- function(genome_build, ld_block_file = NULL) {
  if (!is.null(ld_block_file)) {
    block_file <- normalizePath(ld_block_file, mustWork = TRUE)
  } else {
    block_filename <- if (identical(genome_build, "hg19"))
      "LAVA_s2500_m25_f1_w200.blocks" else "deCODE_EUR_LD_blocks.bed"
    block_file <- file.path(
      system.file("extdata", package = "KnockoffPipeline"), block_filename
    )
  }
  if (!file.exists(block_file)) stop("LD block reference not found: ", block_file)

  header <- names(data.table::fread(block_file, nrows = 0L, showProgress = FALSE))
  required <- c("chr", "start", "stop")
  if (!all(required %in% header))
    stop("LD block file must contain columns: ", paste(required, collapse = ", "))
  block_file
}

#' Stable identifier for a reference file used in saved-knockoff validation
#' @keywords internal
.reference_file_id <- function(path) {
  paste0(basename(path), ":md5:", unname(tools::md5sum(path)))
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
  row_map   <- match(as.character(target_ids), as.character(saved_ids))

  if (any(is.na(row_map)))
    stop(sum(is.na(row_map)), " target sample(s) not found in saved knockoff.")

  # Preserve file-backed matrices without copying when rows already match.
  rows_already_match <- identical(as.character(target_ids),
                                  as.character(saved_ids))

  # Single_Window: G_k is a list of M matrices (n x p).
  if (!is.null(ko_obj$G_k) && !rows_already_match) {
    ko_obj$G_k <- lapply(ko_obj$G_k, function(m) m[row_map, , drop = FALSE])
  }

  # Gene_Centric: knockoffs are M x n x p arrays
  if (!is.null(ko_obj$G_gene_buffer_knockoff)) {
    ko_obj$G_gene_buffer_knockoff <- ko_obj$G_gene_buffer_knockoff[, row_map, , drop = FALSE]
  }
  if (!is.null(ko_obj$G_EnhancerAll_knockoff) && length(dim(ko_obj$G_EnhancerAll_knockoff)) == 3L) {
    ko_obj$G_EnhancerAll_knockoff <- ko_obj$G_EnhancerAll_knockoff[, row_map, , drop = FALSE]
  }

  ko_obj$sample_ids <- target_ids
  ko_obj
}

#' Check whether saved knockoffs can be reused or must be regenerated
#'
#' Model-X knockoffs are generated and preprocessed for a fixed sample set.
#' Reuse therefore requires the same character-IID set; a different order is
#' harmless because rows are explicitly aligned.
#'
#' @param current_ids Numeric/character vector of sample IDs in the current run.
#' @param saved_ids   Numeric/character vector of sample IDs stored in the RDS.
#' @return \code{TRUE} when regeneration is required, \code{FALSE} when only
#'   row reordering may be needed.
#' @keywords internal
.need_regenerate_samples <- function(current_ids, saved_ids) {
  if (is.null(current_ids) || is.null(saved_ids)) return(FALSE)
  current_ids <- as.character(current_ids)
  saved_ids <- as.character(saved_ids)
  length(current_ids) != length(saved_ids) || !setequal(current_ids, saved_ids)
}
