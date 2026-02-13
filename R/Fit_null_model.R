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
                                sparse_grm_id_file = NULL) {
  
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
  trait_type <- ifelse(outcome_type == "D", 'binary', 'quantitative')

  # 准备SAIGE参数列表
  saige_args <- list(
    plinkFile = plink_file,
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
  if (!dir.exists(output_prefix)) {
    dir.create(output_prefix, recursive = TRUE)
  }
  # 添加稀疏GRM参数
  if (!is.null(sparse_grm_file)) {
      saige_args$sparseGRMFile <- sparse_grm_file
      saige_args$sparseGRMSampleIDFile <- sparse_grm_id_file
  }else{
    thin.path <- file.path(output_prefix,"thinned")
    system(sprintf("%s --bfile %s --thin 0.001 --make-bed --out %s --silent",
                  plink_prefix, plink_file, thin.path))
    createSparseGRM(bedFile = paste0(thin.path, ".bed"),
                bimFile = paste0(thin.path, ".bim"),
                famFile = paste0(thin.path, ".fam"),
                outputPrefix= file.path(output_prefix,"GRM"),
                nThreads = n_threads)
    saige_args$sparseGRMFile <- file.path(output_prefix,"GRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
    saige_args$sparseGRMSampleIDFile <- file.path(output_prefix,"GRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
  }
  # 执行SAIGE null model拟合
  saige_args$plinkFile <- file.path(output_prefix,"thinned")
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
      sparseSigma = readMM(file.path(output_prefix,"GRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")),
      ratio = as.numeric(ratio))
  return(results)
}