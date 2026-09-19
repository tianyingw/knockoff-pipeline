.build_glm_covariates <- function(data, covar_cols = NULL,
                                  cat_covar_cols = NULL) {
   all_cols <- c(covar_cols, cat_covar_cols)
   if (length(all_cols) == 0L) return(NULL)
   if (anyDuplicated(all_cols))
      stop("A covariate cannot appear in both covar_cols and cat_covar_cols.")

   data <- as.data.frame(data, stringsAsFactors = FALSE)
   frame <- data[, all_cols, drop = FALSE]
   for (nm in covar_cols) {
      if (!is.numeric(frame[[nm]]))
         stop("Continuous covariate '", nm, "' must be numeric.")
      if (any(!is.finite(frame[[nm]])))
         stop("Continuous covariate '", nm, "' contains non-finite values.")
   }
   for (nm in cat_covar_cols) {
      frame[[nm]] <- factor(as.character(frame[[nm]]))
      if (nlevels(frame[[nm]]) < 2L)
         stop("Categorical covariate '", nm,
              "' has fewer than two observed levels.")
   }

   design <- stats::model.matrix(~ ., data = frame)
   design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
   if (ncol(design) == 0L) NULL else design
}

Fit_null_model<-function(Y, X=NULL, id=NULL, out_type="C", resampling=FALSE,B=1000){
   
   Y<-as.matrix(Y);n<-nrow(Y)
   
   if(length(X)!=0){
      X <- as.matrix(X)
      if (nrow(X) != n) stop("Covariate rows do not match the phenotype length.")
      if (any(!is.finite(X))) stop("The covariate design contains non-finite values.")
   }else{X<-NULL}
   # Use an orthonormal basis for the complete design (including the
   # intercept), retaining only nonzero singular directions. This preserves
   # the fitted covariate space while avoiding singular cross-products from
   # redundant dummy variables or constant columns.
   design <- cbind(`(Intercept)` = rep(1, n), X)
   sx <- svd(design)
   tol <- max(dim(design)) * max(sx$d, 0) * .Machine$double.eps
   keep <- which(sx$d > tol)
   X0 <- sx$u[, keep, drop = FALSE]
   
   if(out_type=="C"){nullglm<-glm(Y~0+X0,family=gaussian)}
   if(out_type=="D"){nullglm<-glm(Y~0+X0,family=binomial)}
   
   if (length(id)==0){id<-1:n}
   
   mu<-nullglm$fitted.values;Y.res<-Y-mu;
   #permute the residuals for B times when sample size is small
   re.Y.res=NULL
   if(resampling==TRUE){
      index<-sapply(1:B,function(x)sample(1:length(Y)));temp.Y.res<-Y.res[as.vector(index)]
      re.Y.res<-matrix(temp.Y.res,length(Y),B)
   }
   
   #prepare invserse matrix for covariates
   if(out_type=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),length(Y))}
   inv.X0<-solve(t(X0)%*%(v*X0))
   inv.vX0<-inv.X0
   
   #prepare the preliminary features
   result.null.model<-list(Y=Y,id=id,n=n,mu=mu,res=Y.res,v=v,
                           X0=X0,nullglm=nullglm,out_type=out_type,
                           re.Y.res=re.Y.res,inv.X0=inv.X0,inv.vX0=inv.vX0)
   return(result.null.model)
}

.align_sparse_grm <- function(grm, grm_ids, model_ids) {
  grm_ids <- .as_sample_id(grm_ids, "sparse GRM sample IDs")
  model_ids <- .as_sample_id(model_ids, "SAIGE null-model sample IDs")
  if (nrow(grm) != length(grm_ids) || ncol(grm) != length(grm_ids))
    stop("Sparse GRM dimensions do not match its sample-ID file.")

  grm_order <- match(model_ids, grm_ids)
  if (anyNA(grm_order)) {
    missing_ids <- model_ids[is.na(grm_order)]
    stop(
      length(missing_ids),
      " SAIGE model sample(s) are absent from the sparse GRM ID file. Examples: ",
      paste(utils::head(missing_ids, 5L), collapse = ", ")
    )
  }
  if (identical(grm_order, seq_along(model_ids)) &&
      nrow(grm) == length(model_ids)) return(grm)
  grm[grm_order, grm_order, drop = FALSE]
}

.saige_covariate_args <- function(covar_cols = NULL, cat_covar_cols = NULL) {
  out <- list()
  all_covar_cols <- unique(c(covar_cols, cat_covar_cols))
  if (length(all_covar_cols) > 0L)
    out$covarColList <- as.character(all_covar_cols)
  if (length(cat_covar_cols) > 0L)
    out$qCovarCol <- as.character(cat_covar_cols)
  out
}

Fit_null_model_GLMM <- function(plink_file,
                                pheno_file,
                                pheno_col,
                                plink_prefix,
                                outcome_type = "C",
                                sample_id_col = NULL,
                                covar_cols = NULL,
                                cat_covar_cols = NULL,
                                output_prefix = "saige_output",
                                n_threads = 4,
                                sparse_grm_file = NULL,
                                sparse_grm_id_file = NULL,
                                thin_target_markers = 5000L,
                                num_random_marker_for_sparse_kin = 1000L,
                                min_maf_for_grm = 0.01,
                                max_missing_rate_for_grm = 0.15,
                                relatedness_cutoff = 0.125,
                                random_seed = NULL) {
  # if ("package:SAIGE" %in% search()) {
  #   try(closeGenoFile_plink(), silent = TRUE)
  # }
  
  # 验证必需参数
  if (missing(plink_file) || missing(pheno_file) || missing(pheno_col)) {
    stop("plink_file, pheno_file, and pheno_col are required parameters")
  }
  if (!file.exists(paste0(plink_file, ".bed"))) {
    stop("PLINK file not found: ", plink_file)
  }
  if (!file.exists(pheno_file)) {
    stop("Phenotype file not found: ", pheno_file)
  }
  if (!outcome_type %in% c("D", "C")) {
    stop("outcome must be 'D' or 'C'")
  }
  if (xor(is.null(sparse_grm_file), is.null(sparse_grm_id_file))) {
    stop("sparse_grm_file and sparse_grm_id_file must be supplied together")
  }
  if (!is.numeric(thin_target_markers) || length(thin_target_markers) != 1L || thin_target_markers < 1) {
    stop("'thin_target_markers' must be a positive integer")
  }
  if (!is.numeric(num_random_marker_for_sparse_kin) || length(num_random_marker_for_sparse_kin) != 1L || num_random_marker_for_sparse_kin < 1) {
    stop("'num_random_marker_for_sparse_kin' must be a positive integer")
  }
  if (!is.null(random_seed)) .derive_unit_seed(random_seed, "SAIGE-validation")
  trait_type <- ifelse(outcome_type == "D", 'binary', 'quantitative')
  output_prefix <- normalizePath(output_prefix, winslash = "/", mustWork = FALSE)
  grm_prefix <- file.path(output_prefix, "GRM")
  thin_path <- file.path(output_prefix, "thinned")
  total_markers <- nrow(data.table::fread(paste0(plink_file, ".bim"), header = FALSE, select = 1L, showProgress = FALSE))
  analysis_prefix <- plink_file

  # 准备SAIGE参数列表
  saige_args <- list(
    plinkFile = analysis_prefix,
    phenoFile = pheno_file,
    phenoCol = pheno_col,
    traitType = trait_type,
    sampleIDColinphenoFile = sample_id_col,
    outputPrefix = output_prefix,
    nThreads = n_threads,
    useSparseGRMtoFitNULL = TRUE,
    usePCGwithSparseGRM = TRUE,       # iterative PCG solver (avoids SuperLU OOM)
    skipVarianceRatioEstimation = FALSE,
    IsOverwriteVarianceRatioFile = TRUE
  )
  
  # 添加协变量
  # SAIGE requires every categorical covariate in qCovarCol to also be
  # present in covarColList.  Keep the public distinction between continuous
  # and categorical columns, but pass their union to SAIGE's covariate list.
  saige_args <- c(
    saige_args,
    .saige_covariate_args(covar_cols, cat_covar_cols)
  )
  # 创建输出目录
  # If user supplies a GRM, do NOT delete the output directory — the GRM may
  # be inside it.  If no GRM is supplied, clean up any stale SAIGE output so
  # that createSparseGRM starts from a clean state.
  if (is.null(sparse_grm_file)) {
    if (dir.exists(output_prefix)) {
      unlink(output_prefix, recursive = TRUE, force = TRUE)
    }
  }
  dir.create(output_prefix, recursive = TRUE, showWarnings = FALSE)

  if (total_markers > as.integer(thin_target_markers)) {
    thin_fraction <- as.integer(thin_target_markers) / total_markers
    message(sprintf(
      "Thinning PLINK markers from %d to about %d (fraction %.6f).",
      total_markers, as.integer(thin_target_markers), thin_fraction
    ))
    seed_arg <- if (is.null(random_seed)) "" else
      paste("--seed", as.integer(random_seed))
    thin_status <- system(sprintf(
      "%s --bfile %s --thin %s %s --make-bed --out %s --silent",
      shQuote(plink_prefix),
      shQuote(plink_file),
      format(thin_fraction, scientific = FALSE, trim = TRUE),
      seed_arg,
      shQuote(thin_path)
    ))
    if (!identical(thin_status, 0L) ||
        !all(file.exists(paste0(thin_path, c(".bed", ".bim", ".fam"))))) {
      stop("PLINK failed while creating the marker-thinned dataset for SAIGE.")
    }
    analysis_prefix <- thin_path
  } else {
    message(sprintf(
      "PLINK file has %d markers only; skipping thinning and using the original dataset.",
      total_markers
    ))
  }
  # 添加稀疏GRM参数
  if (!is.null(sparse_grm_file)) {
      saige_args$sparseGRMFile <- sparse_grm_file
      saige_args$sparseGRMSampleIDFile <- sparse_grm_id_file
  }else{
    create_sparse_grm_once <- function(prefix) {
      createSparseGRM(
        bedFile = paste0(prefix, ".bed"),
        bimFile = paste0(prefix, ".bim"),
        famFile = paste0(prefix, ".fam"),
        outputPrefix = grm_prefix,
        numRandomMarkerforSparseKin = as.integer(num_random_marker_for_sparse_kin),
        relatednessCutoff = relatedness_cutoff,
        nThreads = n_threads,
        minMAFforGRM = min_maf_for_grm,
        maxMissingRateforGRM = max_missing_rate_for_grm
      )
    }

    tryCatch(
      create_sparse_grm_once(analysis_prefix),
      error = function(e) {
        if (analysis_prefix == plink_file) {
          stop(e)
        }
        warning(
          "Sparse GRM creation failed on the thinned dataset; retrying with the original PLINK dataset. Original error: ",
          conditionMessage(e)
        )
        analysis_prefix <<- plink_file
        create_sparse_grm_once(analysis_prefix)
      }
    )
    saige_args$sparseGRMFile <- paste0(grm_prefix, "_relatednessCutoff_", relatedness_cutoff, "_", as.integer(num_random_marker_for_sparse_kin), "_randomMarkersUsed.sparseGRM.mtx")
    saige_args$sparseGRMSampleIDFile <- paste0(grm_prefix, "_relatednessCutoff_", relatedness_cutoff, "_", as.integer(num_random_marker_for_sparse_kin), "_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
  }

  # Validate that GRM exists (either user-supplied or auto-generated)
  if (!is.null(saige_args$sparseGRMFile)) {
    if (!file.exists(saige_args$sparseGRMFile)) {
      stop("Sparse GRM file not found: ", saige_args$sparseGRMFile)
    }
    message("Using sparse GRM: ", saige_args$sparseGRMFile)
  }

  if (!file.exists(saige_args$sparseGRMSampleIDFile)) {
    stop("Sparse GRM sample-ID file not found: ",
         saige_args$sparseGRMSampleIDFile)
  }
  sparse_grm_check <- Matrix::readMM(saige_args$sparseGRMFile)
  grm_id_table <- data.table::fread(
    saige_args$sparseGRMSampleIDFile, header = FALSE,
    keepLeadingZeros = TRUE, showProgress = FALSE
  )
  if (ncol(grm_id_table) < 1L)
    stop("Sparse GRM sample-ID file is empty: ",
         saige_args$sparseGRMSampleIDFile)
  # SAIGE-generated files contain one sample-ID column. Accept a FID/IID-style
  # two-column file as well, using the final (IID) column.
  grm_ids <- .as_sample_id(
    grm_id_table[[ncol(grm_id_table)]], "sparse GRM sample IDs"
  )
  if (nrow(sparse_grm_check) != length(grm_ids) ||
      ncol(sparse_grm_check) != length(grm_ids)) {
    stop("Sparse GRM dimensions do not match its sample-ID file.")
  }
  # Bundled SAIGE expects one IID per line. Accept a conventional two-column
  # FID/IID file at the public interface, but normalize it before delegation.
  normalized_grm_ids <- tempfile(
    "KnockoffPipeline_sparseGRM_IID_", fileext = ".txt"
  )
  writeLines(grm_ids, normalized_grm_ids)
  on.exit(unlink(normalized_grm_ids, force = TRUE), add = TRUE)
  saige_args$sparseGRMSampleIDFile <- normalized_grm_ids
  # 执行SAIGE null model拟合
  saige_args$plinkFile <- analysis_prefix
  rda_file <- paste0(output_prefix, ".rda")
  do.call(fitNULLGLMM, saige_args)
  load(rda_file)

  ratio <- as.matrix(read.table(paste0(output_prefix,".varianceRatio.txt")))[1,1]

  # fitNULLGLMM determines the actual analyzed sample order. Never overwrite
  # that order with the raw phenotype-file order, which may contain excluded or
  # differently ordered rows.
  if (is.null(modglmm$sampleID) || length(modglmm$sampleID) == 0L) {
    stop("SAIGE null model did not return analyzed sample IDs; safe genotype/GRM alignment is impossible.")
  }
  model_ids <- .as_sample_id(modglmm$sampleID, "SAIGE null-model sample IDs")
  if (!is.null(sample_id_col)) {
    pheno_ids <- data.table::fread(
      pheno_file, select = sample_id_col, keepLeadingZeros = TRUE,
      showProgress = FALSE
    )[[1L]]
    pheno_ids <- .as_sample_id(pheno_ids, "phenotype sample IDs")
    if (any(!model_ids %in% pheno_ids))
      stop("SAIGE returned sample IDs absent from the phenotype file.")
  }

  sparse_grm_check <- .align_sparse_grm(
    sparse_grm_check, grm_ids = grm_ids, model_ids = model_ids
  )
  modglmm$sampleID <- model_ids

  diag_entries <- sum(Matrix::diag(sparse_grm_check) != 0)
  off_diag_nnz <- Matrix::nnzero(sparse_grm_check) - diag_entries
  if (off_diag_nnz <= 0) {
    stop(
      "No related sample pairs remain in the analyzed-sample sparse GRM at relatedness_cutoff = ",
      relatedness_cutoff,
      ". Use sample_uncorrelated = TRUE to select the standard GLM path."
    )
  }

  modglmm$traitType <- ifelse(modglmm$traitType == "binary", 'D', 'C')
  results <- list(
      result.null.model.GLMM = modglmm,
      sparseSigma = sparse_grm_check,   # reuse the GRM already loaded at validation step
      ratio = as.numeric(ratio))
  return(results)
}
