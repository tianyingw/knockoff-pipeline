#' Run Knockoff Genome-wide Pipeline
#'
#' Main entry function for knockoffPipeline.
#'
#' @param outdir Output directory.
#' @param test_type Either "Single_Window" or "Gene_Centric".
#' @param pheno_file Path to phenotype file.
#' @param geno_file PLINK genotype prefix (without extension).
#' @param phenotype Column name of phenotype.
#' @param pheno_id Optional sample ID column.
#' @param covariates Character vector of covariate column names.
#' @param user_cores Number of cores for parallelization within chromosome.
#' @param sliding_window_length Sliding window size.
#' @param geno_missing_imputation Imputation method.
#' @param plink_path Path to plink executable.
#' @param M Number of knockoffs.
#' @param genome_build "hg19" or "hg38".
#' @param sample_uncorrelated Logical. TRUE = standard null model.
#' @param grm_file Sparse GRM file for GLMM.
#' @param grm_id_file Sparse GRM ID file.
#' @param fdr Target FDR level.
#' #' @param chromosomes Chromosomes to analyze. NULL or "all" = all autosomes in reference.
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
#' @export
run_pipeline <- function(
  outdir,
  test_type,
  pheno_file,
  geno_file,
  phenotype,
  pheno_id = NULL,
  covariates = NULL,
  user_cores = 1,
  sliding_window_length = NULL,
  geno_missing_imputation = "fixed",
  plink_path = "plink",
  M = 5,
  genome_build = "hg19",
  sample_uncorrelated = TRUE,
  grm_file = NULL,
  grm_id_file = NULL,
  fdr = 0.1,
  chromosomes = 1:22
) {

  ## -----------------------------
  ## 1. Validation
  ## -----------------------------

  if (!test_type %in% c("Single_Window","Gene_Centric"))
    stop("test_type must be Single_Window or Gene_Centric")

  if (!genome_build %in% c("hg19","hg38"))
    stop("genome_build must be hg19 or hg38")

  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)

  if (user_cores > 1)
    Sys.setenv(MKL_NUM_THREADS = 1)

  ## -----------------------------
  ## 2. Fit null model 
  ## -----------------------------

  pheno <- data.table::fread(pheno_file)

  pheno.value <- unique(pheno[[phenotype]])
  is.binary <- length(pheno.value)==2 &&
               all(sort(pheno.value)==c(0,1))

  out_type <- if(is.binary) "D" else "C"

  if (sample_uncorrelated) {

    if (!is.null(pheno_id)) {
      id_vec <- pheno[[pheno_id]]
      nullobj <- Fit_null_model(
        Y = pheno[[phenotype]],
        X = as.matrix(pheno[, covariates, with=FALSE]),
        id = as.numeric(id_vec),
        out_type = out_type
      )
    } else {
      nullobj <- Fit_null_model(
        Y = pheno[[phenotype]],
        X = as.matrix(pheno[, covariates, with=FALSE]),
        out_type = out_type
      )
    }

  } else {

    nullobj_results <- Fit_null_model_GLMM(
      geno_file,
      pheno_file,
      phenotype,
      plink_path,
      outcome_type = out_type,
      sample_id_col = pheno_id,
      covar_cols = covariates,
      sparse_grm_file = grm_file,
      sparse_grm_id_file = grm_id_file
    )
  }

  rm(pheno); gc()

  ## genotype id
  if (is.null(pheno_id)) {
    Gsub.id <- NULL
  } else {
    Gsub.id <- data.table::fread(paste0(geno_file,".fam"),header=FALSE)[[2]]
  }

  chr_vector <- intersect(as.numeric(chromosomes), 1:22)

  mid_dir <- file.path(outdir,"mid")
  if (!dir.exists(mid_dir))
    dir.create(mid_dir)

  ############################################################
  ################  Single Window ############################
  ############################################################

  if (test_type == "Single_Window") {

    block_file <- file.path(
      system.file("extdata", package="KnockoffPipeline"),
      ifelse(genome_build=="hg19",
             "LAVA_s2500_m25_f1_w200.blocks",
             "deCODE_EUR_LD_blocks.bed")
    )

    blocks <- data.table::fread(block_file)
    unique_chr <- intersect(sort(unique(blocks$chr)), chr_vector)

    for (c in unique_chr) {

      message("Processing chromosome ", c)

      block_chr <- blocks[blocks$chr==c]
      sub_seq_id <- seq_len(nrow(block_chr))

      safe_fun <- function(kk) {
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
          error=function(e) NULL
        )
      }

      out <- mclapply(sub_seq_id, safe_fun, mc.cores=user_cores)
      out <- Filter(Negate(is.null), out)

      if (length(out)==0) next

      single_chr <- data.table::rbindlist(
        lapply(out, `[[`, "result.single"), fill=TRUE)

      window_chr <- data.table::rbindlist(
        lapply(out, `[[`, "result.window"), fill=TRUE)

      data.table::fwrite(single_chr,
                         file.path(mid_dir,
                           paste0("Single_mid_results_chr",c,".txt")),
                         sep="\t")

      data.table::fwrite(window_chr,
                         file.path(mid_dir,
                           paste0("Window_mid_results_chr",c,".txt")),
                         sep="\t")

      rm(out,single_chr,window_chr); gc()
    }

    ## merge
    result.single.all <- data.table::rbindlist(
      lapply(chr_vector,function(c)
        data.table::fread(file.path(mid_dir,
          paste0("Single_mid_results_chr",c,".txt")))),
      fill=TRUE)

    result.window.all <- data.table::rbindlist(
      lapply(chr_vector,function(c)
        data.table::fread(file.path(mid_dir,
          paste0("Window_mid_results_chr",c,".txt")))),
      fill=TRUE)

    summary <- KS_summary(
      as.matrix(result.window.all),
      as.matrix(result.single.all),
      M,
      fdr=fdr
    )

    result_all <- summary[,c("chr","start","end",
                             "Qvalue","W_KS",
                             "W_Threshold","detect")]

    data.table::fwrite(result_all,
      file.path(outdir,"Single_Window_results.csv"))

    plot_manhattan(result_all,outdir,
                   "manhattan_plot_single.png")
  }

  ############################################################
  ################  Gene Centric #############################
  ############################################################

  if (test_type == "Gene_Centric") {

    gene_file <- file.path(
      system.file("extdata",package="KnockoffPipeline"),
      genome_build,
      paste0("coding.genes.TSS.",genome_build,".tsv")
    )

    genes_info <- data.table::fread(gene_file)
    genes_info$chr <- as.numeric(gsub("[^0-9]","",genes_info$chr))
    unique_chr <- intersect(sort(unique(genes_info$chr)), chr_vector)

    for (c in unique_chr) {

      message("Processing chromosome ",c)

      chr_genes <- genes_info[chr==c]
      sub_seq_id <- seq_len(nrow(chr_genes))

      safe_fun <- function(kk){

        tryCatch({

          if (sample_uncorrelated) {

            run_single_gene(
              genes = chr_genes,
              kk = kk,
              geno.file = geno_file,
              obj_nullmodel = nullobj,
              window_length = sliding_window_length,
              plink_prefix = plink_path,
              M = M,
              genome_build = genome_build,
              Gsub.id = Gsub.id
            )

          } else {

            run_single_gene(
              genes = chr_genes,
              kk = kk,
              geno.file = geno_file,
              obj_nullmodel =
                nullobj_results$result.null.model.GLMM,
              window_length = sliding_window_length,
              plink_prefix = plink_path,
              M = M,
              genome_build = genome_build,
              Gsub.id = Gsub.id,
              sparseSigma = nullobj_results$sparseSigma,
              ratio = nullobj_results$ratio
            )
          }

        }, error=function(e) NULL)
      }

      out <- mclapply(sub_seq_id,
                      safe_fun,
                      mc.cores=user_cores)

      out <- Filter(Negate(is.null), out)
      if (length(out)==0) next

      result_chr <- data.table::rbindlist(
        lapply(out, `[[`, "result"), fill=TRUE)

      data.table::fwrite(result_chr,
        file.path(mid_dir,
          paste0("GeneCentric_mid_results_chr",c,".txt")),
        sep="\t")

      rm(out,result_chr); gc()
    }

    result.all <- data.table::rbindlist(
      lapply(unique_chr,function(c)
        data.table::fread(file.path(mid_dir,
          paste0("GeneCentric_mid_results_chr",c,".txt")))),
      fill=TRUE)

    summary <- GeneScan3DKnock_Summary(
      result.all,
      M=M,
      fdr=fdr
    )

    result_all <- summary[,c("chr","gene_id",
                             "gene_start","gene_end",
                             "Qvalue","W",
                             "W_Threshold","detect")]

    colnames(result_all) <-
      c("chr","gene_id","start","end",
        "Qvalue","W_KS","W_Threshold","detect")

    data.table::fwrite(result_all,
      file.path(outdir,"GeneCentric_results.csv"))

    plot_manhattan(result_all,
                   outdir,
                   "manhattan_plot_gene.png")
  }

  invisible(TRUE)
}
