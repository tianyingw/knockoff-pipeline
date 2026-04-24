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