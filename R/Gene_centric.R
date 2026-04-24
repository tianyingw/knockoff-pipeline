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
  abc_df,
  gh_df,
  sparseSigma=NULL,
  ratio=NULL,
  user_cores=1
){
  kk_vec <- batch_index[[b]]
  tmpdir <- tempdir()
  chr <- as.numeric(gsub("chr", "", genes[kk_vec[1], chr]))

  gene_buffer_extension <- 5000 + 50000

  start_all <- min(genes[kk_vec, start]) - gene_buffer_extension
  end_all   <- max(genes[kk_vec, end]) + gene_buffer_extension

  batch_prefix <- file.path(tmpdir,
                     sprintf("temp_chr%d_batch_%d_%d",
                             chr,
                             min(kk_vec),
                             max(kk_vec)))

  ## ===== plink IO =====
  system(sprintf(
    "%s --bfile %s --chr %s --from-bp %d --to-bp %d --recode A --out %s --silent",
    plink_prefix, geno.file, chr,
    start_all, end_all,
    batch_prefix), ignore.stdout = TRUE, ignore.stderr = TRUE)

  raw_file <- paste0(batch_prefix, ".raw")
  if (!file.exists(raw_file)) 
    return(NULL)
  
  raw <- data.table::fread(raw_file, data.table = FALSE)

  unlink(paste0(batch_prefix, c(".raw",".log",".nosex")),force=TRUE)

  if (ncol(raw) <= 6) return(NULL)

  message("  Batch ", b, " / ", length(batch_index),
          " (snp ", start_all, "-",
          end_all, ")")

  G_batch <- as.matrix(raw[, -(1:6), drop = FALSE])
  variants_batch <- extract_position_universal(colnames(G_batch))
  rm(raw); gc()

  ## ===== batch 内并行 =====
  safe_fun <- function(kk){

    tryCatch({

      gene_start <- genes[kk, start]
      gene_end   <- genes[kk, end]
      gene_id    <- genes[kk, id]

      ## -------- gene buffer SNP --------
      idx_gene_buffer <- which(variants_batch >= gene_start-5000 & variants_batch <= gene_end+5000)
      idx_gene_surround <- which(variants_batch >= gene_start-gene_buffer_extension & variants_batch <= gene_end+gene_buffer_extension)
      if(length(idx_gene_buffer)<=1) return(NULL)
      print(paste0("Gene ", gene_id, ": ", length(idx_gene_buffer), " SNPs in buffer region, ", length(idx_gene_surround), " SNPs in surrounding region."))

      G_gene <- G_batch[, idx_gene_surround, drop=FALSE]
      gene_buffer.pos <- c(min(variants_batch[idx_gene_buffer]),max(variants_batch[idx_gene_buffer]))

      ## -------- enhancer --------
      abc_enhancers <- abc_df[TargetGene == gene_id,.(start,end)]
      gh_enhancers <- gh_df[gene == gene_id,.(start=GH_start,end=GH_end)]

      enhancers <- unique(rbind(abc_enhancers, gh_enhancers))

      G_EnhancerAll_surround <- NULL
      variants_EnhancerAll_surround <- NULL
      Enhancer.pos <- NULL
      p_EnhancerAll_surround <- NULL
      p_EnhancerAll <- NULL
      R <- 0

      if(nrow(enhancers)>0){

        for(r in seq_len(nrow(enhancers))){

          e_start <- enhancers$start[r]
          e_end   <- enhancers$end[r]
          idx_e_surrond <- which(variants_batch >= e_start-5000 & variants_batch <= e_end+5000)
          idx_e <- which(variants_batch >= e_start & variants_batch <= e_end)

          if(length(idx_e)>5){
            temp.G <- G_batch[, idx_e_surrond, drop=FALSE]
            G_EnhancerAll_surround <-cbind(G_EnhancerAll_surround, temp.G)
            Enhancer.pos <-rbind(Enhancer.pos,c(e_start, e_end))
            variants_EnhancerAll_surround <-c(variants_EnhancerAll_surround,variants_batch[idx_e_surrond])
            p_EnhancerAll_surround <- c(p_EnhancerAll_surround, length(idx_e_surrond))
            p_EnhancerAll <- c(p_EnhancerAll,length(idx_e))

            R <- R + 1
          }
        }
      }

      ## -------- GeneScan --------
      if (is.null(sparseSigma)) {
        print("GeneScan3DKnock")
        print(R)
        full_results <- GeneScan3D.KnockoffGeneration(
          G_gene_buffer_surround=G_gene,
          variants_gene_buffer_surround=variants_batch[idx_gene_surround],
          gene_buffer.pos=gene_buffer.pos,
          R=R,
          G_EnhancerAll_surround=G_EnhancerAll_surround,
          variants_EnhancerAll_surround=variants_EnhancerAll_surround,
          p_EnhancerAll_surround=p_EnhancerAll_surround,
          Enhancer.pos=Enhancer.pos,
          p.EnhancerAll=p_EnhancerAll,
          Z=NULL,Z.promoter=NULL,Z.EnhancerAll=NULL,promoter.pos=NULL,
          window.size=window_length,
          result.null.model=obj_nullmodel,
          M=M,
          Gsub.id=Gsub.id)
        print(str(full_results))

      } else {
        # print("BIGKnock")
        full_results <- GeneScan3D.UKB.GLMM.KnockoffGeneration(
          G_gene_buffer_surround=G_gene,
          variants_gene_buffer_surround=variants_batch[idx_gene_surround],
          gene_buffer.pos=gene_buffer.pos,
          R=R,
          G_EnhancerAll_surround=G_EnhancerAll_surround,
          variants_EnhancerAll_surround=variants_EnhancerAll_surround,
          p_EnhancerAll_surround=p_EnhancerAll_surround,
          Enhancer.pos=Enhancer.pos,
          p.EnhancerAll=p_EnhancerAll,
          Z=NULL,Z.promoter=NULL,Z.EnhancerAll=NULL,promoter.pos=NULL,
          window.size=window_length,
          result.null.model=obj_nullmodel,
          M=M,
          Gsub.id=Gsub.id,
          sparseSigma=sparseSigma,
          ratio=ratio)
      }
      
      results <- data.frame(
        chr = chr,
        gene_id = gene_id,
        gene_start = gene_start,
        gene_end = gene_end,
        GeneScan3D.Cauchy = full_results$GeneScan3D.Cauchy[1],
        t(full_results$GeneScan3D.Cauchy_knockoff[,1,drop=FALSE]),
        stringsAsFactors = FALSE
      )

      colnames(results)[6:ncol(results)] <- paste0("GeneScan3D.Cauchy_knockoff_",1:M)
      # print(str(results))
      return(results)

    }, error=function(e) NULL)
  }

  out <- mclapply(kk_vec,
                  safe_fun,
                  mc.cores=user_cores)

  out <- Filter(Negate(is.null), out)
  # print(str(out))
  rm(G_batch); gc()

  if(length(out)==0) return(NULL)

  return(data.table::rbindlist(out, fill=TRUE))
}

GeneScan3D.UKB.GLMM.KnockoffGeneration <- function(G_gene_buffer_surround=G_gene_buffer_surround,
                                                   variants_gene_buffer_surround=variants_gene_buffer_surround,
                                                   gene_buffer.pos=gene_buffer.pos,R=R,
                                                   G_EnhancerAll_surround=G_EnhancerAll_surround,
                                                   variants_EnhancerAll_surround=variants_EnhancerAll_surround,
                                                   p_EnhancerAll_surround=p_EnhancerAll_surround,
                                                   Enhancer.pos=Enhancer.pos,p.EnhancerAll=p_EnhancerAll,
                                                   window.size=window_length,result.null.model=obj_nullmodel,M=M,
                                                   sparseSigma=sparseSigma,ratio=ratio,
                                                   MAC.threshold=10,MAF.threshold=0.01,Gsub.id=NULL){
   mu<-as.vector(result.null.model$fitted.values)
   Y.res<-as.vector(result.null.model$residuals)
   # print("null model prepared")
   impute.method='fixed'
   ## Prelimanry checking and filtering the variants
   #match phenotype id and genotype id
   if(length(Gsub.id)==0){match.index<-match(result.null.model$sampleID,1:nrow(G_gene_buffer_surround))}else{
      match.index<-match(result.null.model$sampleID,Gsub.id)
   }

   if(mean(is.na(match.index))>0){
      msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
      warning(msg,call.=F)
   }
   #individuals ids are matched with genotype
   G_gene_buffer_surround=Matrix(G_gene_buffer_surround[match.index,])
   #missing genotype imputation
   G_gene_buffer_surround[G_gene_buffer_surround==-9 | G_gene_buffer_surround==9]=NA
   N_MISS=sum(is.na(G_gene_buffer_surround))
   MISS.freq=apply(is.na(G_gene_buffer_surround),2,mean)
   if(N_MISS>0){
      msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G_gene_buffer_surround)/ncol(G_gene_buffer_surround))
      warning(msg,call.=F)
      G_gene_buffer_surround=Impute(G_gene_buffer_surround,impute.method)
   }
   
   #MAF filtering
   MAF<-apply(G_gene_buffer_surround,2,mean)/2 #MAF of nonfiltered variants
   G_gene_buffer_surround[,MAF>0.5 & !is.na(MAF)]<-2-G_gene_buffer_surround[,MAF>0.5 & !is.na(MAF)]
   MAF<-apply(G_gene_buffer_surround,2,mean)/2
   MAC<-apply(G_gene_buffer_surround,2,sum) #minor allele count
   s<-apply(G_gene_buffer_surround,2,sd)
   SNP.index<-which(MAF>0 & s!=0 & !is.na(MAF) & MISS.freq<0.1) 
   
   check.index<-which(MAF>0 & s!=0 & !is.na(MAF)  & MISS.freq<0.1)
   if(length(check.index)<=1){
      warning('Number of variants with missing rate <=10% in the gene is <=1')
   }
   
   G_gene_buffer_surround<-Matrix(G_gene_buffer_surround[,SNP.index])
   variants_gene_buffer_surround_filter=variants_gene_buffer_surround[SNP.index]
   # print("geno prepared")
   ###Generate multiple knockoffs
   n=length(mu)
   colnames(G_gene_buffer_surround)<-extract_position_universal(colnames(G_gene_buffer_surround))
   G_gene_buffer_knockoff<-NULL
   invisible(capture.output(G_gene_buffer_knockoff<-Knockoffgeneration.gene.buffer(G_gene_buffer_surround=G_gene_buffer_surround,
                                                         gene_buffer_start=gene_buffer.pos[1],gene_buffer_end=gene_buffer.pos[2],M=M)))
   # print("knockoff success")
   ##obtain knockoff genotypes for gene buffer region
   positions_gene_buffer=variants_gene_buffer_surround_filter[variants_gene_buffer_surround_filter<=gene_buffer.pos[2]&variants_gene_buffer_surround_filter>=gene_buffer.pos[1]]
   G_gene_buffer=G_gene_buffer_surround[,variants_gene_buffer_surround_filter%in%positions_gene_buffer]

   ## R enhancers ##
   G_EnhancerAll=c()
   p_EnhancerAll=c()
   G_EnhancerAll_knockoff=c()

   if (R!=0){
      for (r in 1:R){
         ##genotype Enhancer_surround_region
         if (r==1){
            G_Enhancer_surround=G_EnhancerAll_surround[,1:cumsum(p_EnhancerAll_surround)[r]]
            positions_Enhancer_surround=variants_EnhancerAll_surround[1:cumsum(p_EnhancerAll_surround)[r]]
         }else{
            G_Enhancer_surround=G_EnhancerAll_surround[,(cumsum(p_EnhancerAll_surround)[r-1]+1):cumsum(p_EnhancerAll_surround)[r]]
            positions_Enhancer_surround=variants_EnhancerAll_surround[(cumsum(p_EnhancerAll_surround)[r-1]+1):cumsum(p_EnhancerAll_surround)[r]]
         }
         
         #individuals ids are matched with genotype
         G_Enhancer_surround=Matrix(G_Enhancer_surround[match.index,])
         #missing genotype imputation
         G_Enhancer_surround[G_Enhancer_surround==-9 | G_Enhancer_surround==9]=NA
         N_MISS=sum(is.na(G_Enhancer_surround))
         MISS.freq=apply(is.na(G_Enhancer_surround),2,mean)
         if(N_MISS>0){
            msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G_Enhancer_surround)/ncol(G_Enhancer_surround))
            warning(msg,call.=F)
            G_Enhancer_surround=Impute(G_Enhancer_surround,impute.method)
         }
         
         #MAF filtering
         MAF<-apply(G_Enhancer_surround,2,mean)/2 #MAF of nonfiltered variants
         G_Enhancer_surround[,MAF>0.5 & !is.na(MAF)]<-2-G_Enhancer_surround[,MAF>0.5 & !is.na(MAF)]
         MAF<-apply(G_Enhancer_surround,2,mean)/2
         MAC<-apply(G_Enhancer_surround,2,sum) #minor allele count
         s<-apply(G_Enhancer_surround,2,sd)
         SNP.index<-which(MAF>0 & s!=0 & !is.na(MAF) & MISS.freq<0.1) 
         
         check.index<-which(MAF>0 & s!=0 & !is.na(MAF)  & MISS.freq<0.1)
         if(length(check.index)<=1){
            warning('Number of variants with missing rate <=10% in the gene is <=1')
         }
         
         G_Enhancer_surround<-Matrix(G_Enhancer_surround[,SNP.index])
         positions_Enhancer_surround_filter=positions_Enhancer_surround[SNP.index]
         colnames(G_Enhancer_surround)<-extract_position_universal(colnames(G_Enhancer_surround))
         G_Enhancer_knockoff<-NULL
         invisible(capture.output(G_Enhancer_knockoff<-Knockoffgeneration.enhancer(G_enhancer_surround=G_Enhancer_surround,
                                                         enhancer_start=as.numeric(Enhancer.pos[r,1]),enhancer_end=as.numeric(Enhancer.pos[r,2]),M=5)))
         
         positions_enhancer=positions_Enhancer_surround_filter[positions_Enhancer_surround_filter<=Enhancer.pos[r,2]&positions_Enhancer_surround_filter>=Enhancer.pos[r,1]]
         G_enhancer=Matrix(G_Enhancer_surround[,positions_Enhancer_surround_filter%in%positions_enhancer])
         G_EnhancerAll=cbind(G_EnhancerAll,G_enhancer)
         
         # p_Enhancer=length(positions_enhancer)
         p_Enhancer = dim(G_Enhancer_knockoff)[3]
         p_EnhancerAll=c(p_EnhancerAll,p_Enhancer)
         G_EnhancerAll_knockoff=abind::abind(G_EnhancerAll_knockoff,G_Enhancer_knockoff)
      }
   }
   print("enhancer prepared")
    ####GeneScan3D.UKB.GLMM: conduct gene-based test on the gene buffer region, adding R enhancers ################
   ##original p-values
   tmp<-NULL
   # invisible(capture.output(tmp<-GeneScan3D.UKB.GLMM(G=G_gene_buffer,G.EnhancerAll=G_EnhancerAll,R=R,
   #                                     p_Enhancer=p_EnhancerAll,window.size=window.size,pos=positions_gene_buffer,
   #                                     MAC.threshold=MAC.threshold,MAF.threshold=MAF.threshold,Gsub.id=row.names(G_gene_buffer),
   #                                     result.null.model.GLMM=result.null.model,outcome=result.null.model$traitType,
   #                                     sparseSigma=sparseSigma,ratio=ratio)$GeneScan3D.Cauchy.pvalue))
   X<-result.null.model$X #covariates include intercept
   tmp<-GeneScan3D.UKB.GLMM(G=G_gene_buffer,Z=NULL,G.promoter=NULL,Z.promoter=NULL,
                              G.EnhancerAll=G_EnhancerAll,Z.EnhancerAll=NULL,R=R,
                              p_Enhancer=p_EnhancerAll,window.size=window.size,pos=positions_gene_buffer,
                              MAC.threshold=MAC.threshold,MAF.threshold=MAF.threshold,Gsub.id=Gsub.id[match.index],
                              result.null.model.GLMM=result.null.model,outcome=result.null.model$traitType,
                              sparseSigma=sparseSigma,ratio=ratio)$GeneScan3D.Cauchy.pvalue
   GeneScan3D.Cauchy = tmp
   print("GeneScan3D success")
   #M knockoff p-values 
   GeneScan3D.Cauchy_knockoff=matrix(NA,nrow=M,ncol=3)
   for (k in 1:M){
      G_gene_buffer_knockoff_k=G_gene_buffer_knockoff[k,,]
      tmp<-NULL
      invisible(capture.output(tmp<-GeneScan3D.UKB.GLMM(G=G_gene_buffer_knockoff_k,G.EnhancerAll=G_EnhancerAll_knockoff[k,,], R=R,
                                                                  p_Enhancer=p_EnhancerAll,window.size=window.size,pos=positions_gene_buffer,
                                                                  MAC.threshold=MAC.threshold,MAF.threshold=MAF.threshold,Gsub.id=Gsub.id[match.index],
                                                                  result.null.model.GLMM=result.null.model,outcome=result.null.model$traitType,
                                                                  sparseSigma=sparseSigma,ratio=ratio)$GeneScan3D.Cauchy.pvalue))
      GeneScan3D.Cauchy_knockoff[k,] = tmp
   }
   print("GeneScan3D knockoff success")
   return(list(GeneScan3D.Cauchy=GeneScan3D.Cauchy,GeneScan3D.Cauchy_knockoff=GeneScan3D.Cauchy_knockoff))
}

extract_position_universal <- function(col_names) {
  positions <- sapply(col_names, function(col) {
    # 提取所有数字序列
    numbers <- regmatches(col, gregexpr("\\d+", col))[[1]]
    
    if (length(numbers) == 0) {
      return(NA)  # 没有数字
    } else if (length(numbers) == 1) {
      # 只有一个数字，直接使用
      return(as.numeric(numbers))
    } else {
      # 有多个数字，使用启发式规则
      # 规则1: 优先选择看起来像基因组位置的数字（3-9位）
      candidate <- numbers[nchar(numbers) >= 3 & nchar(numbers) <= 9]
      if (length(candidate) > 0) {
        return(as.numeric(tail(candidate, 1)))  # 取最后一个符合条件的
      }
      
      # 规则2: 选择最大的数字（通常是位置）
      return(as.numeric(max(numbers)))
    }
  })
  
  return(as.numeric(positions))
}

 
preprocess_for_GeneScan3DKnock <- function(p0, p_ko, M) {
  p0 <- as.numeric(p0)
  p_ko <- as.matrix(p_ko)
  
#   # 处理极端 p 值
#   p0[p0 == 0] <- min(p0[p0 > 0]) / 10
#   p0[p0 == 1] <- 1 - 1e-10
  
#   # 对 knockoff p 值做同样处理
#   for (j in 1:ncol(p_ko)) {
#     col_vals <- p_ko[, j]
#     col_vals[col_vals == 0] <- min(col_vals[col_vals > 0], na.rm = TRUE) / 10
#     col_vals[col_vals == 1] <- 1 - 1e-10
#     p_ko[, j] <- col_vals
#   }
  
  # NA
  p0[is.na(p0)] <- 0.5
  p_ko[is.na(p_ko)] <- 0.5
  
  return(list(p0 = p0, p_ko = p_ko))
}

GeneScan3DKnock_Summary <- function(result, M, fdr = 0.1) {
  result <- as.data.frame(result)
  result <- result[order(result[,4]),]
  result <- result[order(result[,3]),]
  
  p0 <- as.numeric(result[,5])
  pk <- as.matrix(result[,6:(5+M), drop=F])
  
  preprocessed <- preprocess_for_GeneScan3DKnock(p0, pk, M)
  p0_clean <- preprocessed$p0
  pk_clean <- preprocessed$p_ko
   # p0_clean <- p0
   # pk_clean <- pk
  result.GeneScan3DKnock <- GeneScan3DKnock(
    M = M, 
    p0 = p0_clean, 
    p_ko = pk_clean, 
    gene_id = result[,2], 
    fdr = fdr
  )
  
  result$W <- result.GeneScan3DKnock$W
  result$W_Threshold <- rep(result.GeneScan3DKnock$W.threshold, nrow(result))
  result$Qvalue <- result.GeneScan3DKnock$Qvalue
  result$detect <- result.GeneScan3DKnock$Qvalue <= fdr
  
  return(result)
}