Fit_null_model<-function(Y, X=NULL, id=NULL, out_type="C", resampling=FALSE,B=1000){
   
   Y<-as.matrix(Y);n<-nrow(Y)
   
   if(length(X)!=0){X0<-svd(as.matrix(X))$u}else{X0<-NULL}
   X0<-cbind(rep(1,n),X0)
   
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
   inv.vX0<-solve(t(X0)%*%(v*X0))
   
   #prepare the preliminary features
   result.null.model<-list(Y=Y,id=id,n=n,mu=mu,res=Y.res,v=v,
                           X0=X0,nullglm=nullglm,out_type=out_type,
                           re.Y.res=re.Y.res,inv.X0=inv.X0,inv.vX0=inv.vX0)
   return(result.null.model)
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
                                relatedness_cutoff = 0.125) {
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
  if (!is.numeric(thin_target_markers) || length(thin_target_markers) != 1L || thin_target_markers < 1) {
    stop("'thin_target_markers' must be a positive integer")
  }
  if (!is.numeric(num_random_marker_for_sparse_kin) || length(num_random_marker_for_sparse_kin) != 1L || num_random_marker_for_sparse_kin < 1) {
    stop("'num_random_marker_for_sparse_kin' must be a positive integer")
  }
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
    skipVarianceRatioEstimation = FALSE,
    IsOverwriteVarianceRatioFile = TRUE
  )
  
  # 添加协变量
  if (!is.null(covar_cols)) {
    saige_args$covarColList <- as.character(covar_cols)
  }
  if (!is.null(cat_covar_cols)) {
    saige_args$qCovarCol <- as.character(cat_covar_cols)
  }
  # 创建输出目录
  if (dir.exists(output_prefix)) {
    unlink(output_prefix, recursive = TRUE, force = TRUE)
  }
  dir.create(output_prefix, recursive = TRUE, showWarnings = FALSE)

  if (total_markers > as.integer(thin_target_markers)) {
    thin_fraction <- as.integer(thin_target_markers) / total_markers
    message(sprintf(
      "Thinning PLINK markers from %d to about %d (fraction %.6f).",
      total_markers, as.integer(thin_target_markers), thin_fraction
    ))
    system(sprintf(
      "%s --bfile %s --thin %s --make-bed --out %s --silent",
      shQuote(plink_prefix),
      shQuote(plink_file),
      format(thin_fraction, scientific = FALSE, trim = TRUE),
      shQuote(thin_path)
    ))
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

  sparse_grm_check <- Matrix::readMM(saige_args$sparseGRMFile)
  diag_entries <- sum(Matrix::diag(sparse_grm_check) != 0)
  off_diag_nnz <- Matrix::nnzero(sparse_grm_check) - diag_entries
  if (off_diag_nnz <= 0) {
    stop(
      "No related sample pairs were detected in the sparse GRM at relatedness_cutoff = ",
      relatedness_cutoff,
      ". Samples appear unrelated; use sample_uncorrelated = TRUE instead of SAIGE/GLMM."
    )
  }
  # 执行SAIGE null model拟合
  saige_args$plinkFile <- analysis_prefix
  rda_file <- paste0(output_prefix, ".rda")
  do.call(fitNULLGLMM, saige_args)
  load(rda_file)
 
  ratio <- as.matrix(read.table(paste0(output_prefix,".varianceRatio.txt")))[1,1]
  if (length(sample_id_col) == 0){
    modglmm$sampleID = 1:length(modglmm$sampleID)
  }else{
    modglmm$sampleID = fread(pheno_file)[[sample_id_col]]
  }
  
  modglmm$traitType <- ifelse(modglmm$traitType == "binary", 'D', 'C')
  results <- list(
      result.null.model.GLMM = modglmm,
      sparseSigma = readMM(paste0(grm_prefix, "_relatednessCutoff_", relatedness_cutoff, "_", as.integer(num_random_marker_for_sparse_kin), "_randomMarkersUsed.sparseGRM.mtx")),
      ratio = as.numeric(ratio))
  return(results)
}
