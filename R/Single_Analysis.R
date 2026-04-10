run_single_block <- function(blocks, kk, geno.file, obj_nullmodel, window_length, plink_prefix,
                              impute.method, M, Gsub.id) {
  chr   <- blocks[kk, chr]
  start <- blocks[kk, start]
  stop  <- blocks[kk, stop]
  # ---- extract block SNP ----
  tmpdir <- tempdir()
  block_prefix <- file.path(tmpdir, sprintf("temp_chr%d_block%d", chr, kk))

  cmd <- sprintf(
    "%s --bfile %s --chr %s --from-bp %d --to-bp %d --recode A --out %s --silent",
    plink_prefix, geno.file, chr, start, stop, block_prefix
  )

  system(paste(cmd, "2>/dev/null"), ignore.stdout = TRUE)

  # ---- geno df ----
  raw_file <- paste0(block_prefix, ".raw")
  if (!file.exists(raw_file)) 
    return(NULL)
  raw <- data.table::fread(raw_file, data.table = FALSE)
  unlink(paste0(block_prefix, c(".raw", ".log", ".nosex")),force = TRUE)
  
  if (ncol(raw) <= 6)
    return(NULL)
  print(kk)
  df <- as.matrix(raw[, -(1:6), drop = FALSE])
  rm(raw)
  cat(sprintf("chr: %s, start: %d, end: %d, snp count: %d\n", chr, start, stop, ncol(df)))
  
  # ----  KS analysis ----
  invisible(capture.output(
    results <- Single_Window_Analysis(obj_nullmodel, df, chr, 
                                    window_length = window_length, 
                                    M = M, impute.method = impute.method, Gsub.id=Gsub.id)
  ))
  rm(df)
  gc()
  return(results)
}


Single_Window_Analysis<-function(nullobj,geno,chr,window_length=NULL,M=5,thres.single=0.01,thres.ultrarare=25,thres.missing=0.10,midout.dir=NULL,jobtitle=NULL,impute.method='fixed',Gsub.id=NULL,bigmemory=T,leveraging=T,LD.filter=NULL){
  preprocess <- Preprocess(geno = geno, chr = chr, window = window_length,
                             impute.method = impute.method)
  G <- preprocess$G
  pos <- preprocess$pos
  window.bed <- preprocess$window.bed
  G_k <- create.KS(X = G, pos = pos, M = M)
  fit <- KS.chr(result.prelim = nullobj,
                input.X = G,
                window.bed = window.bed,
                input.G_k = G_k,
                M = M,
                thres.single = thres.single,
                Gsub.id = Gsub.id)
  rm(geno,preprocess,G_k,G,pos,window.bed)
  gc()
  return(fit)
}

## preprocess genoinfo: maf>0, missing, cluster, order pos, window bed
Preprocess<-function(geno, chr, window=NULL, thres.maf=0, thres.missing=0.1, impute.method='fixed'){
    # preprocess geno matrix
    G = as.matrix(geno)[,-c(1:6)]
    G = 2-G
    m <- nrow(G)
    if(length(G)==0){
      msg<-'Number of variants in the specified range is 0'
      warning(msg,call.=F)
      next
    }else{
      if(ncol(G)==1){
        msg<-'Number of variants in the specified range is 1'
        warning(msg,call.=F)
        next
      }
    }
    G<-G[,match(unique(colnames(G)),colnames(G)),drop=F]

    # missing rate filtering
    # MISS<-colMeans(G<0 | G>2)
    # G<-G[,MISS<=thres.missing,drop=F]
    # missing genotype imputation
    G[G<0 | G>2]<-NA
    G<-Impute(G,"fixed")

    # filter out constant variants
    s<-apply(G,2,sd)
    G<-G[,s!=0]

    #get positions and reorder G
    pos <- extract_position_universal(colnames(G))
    G<-G[,order(pos),drop=F]
    MAF<-colMeans(G)/2
    G<-as.matrix(G)
    G[,MAF>0.5 & !is.na(MAF)]<-2-G[,MAF>0.5 & !is.na(MAF)]
    MAF<-colMeans(G)/2;MAC<-colSums(G)
    G<-Matrix(G,sparse=T)

    pos <- extract_position_universal(colnames(G))    
    colnames(G) <- pos
    start<-min(pos);end<-max(pos)

    # clustering
    SNP.set <- G
    cor.X<-sparse.cor(Matrix(SNP.set))$cor
    Sigma.distance = as.dist(1 - abs(cor.X))
    fit = hclust(Sigma.distance, method="single")
    corr_max = 0.75
    clusters = cutree(fit, h=1-corr_max)
    cluster.index<-match(unique(clusters),clusters)
    SNP.set<-SNP.set[,cluster.index];MAF<-MAF[cluster.index];MAC<-MAC[cluster.index];pos<-pos[cluster.index]
    unique.index <- match(unique(pos), pos)
    SNP.set<-SNP.set[,unique.index,drop=F];MAF<-MAF[unique.index];MAC<-MAC[unique.index];pos<-pos[unique.index]
    G <- SNP.set
    N.SNP<-ncol(G)

    # Define the window.bed file
    if(length(window)!=0){
      pos.min<-start;pos.max<-end
      window.size <- window
      window.bed<-c();
      for(size in window.size){
          pos.tag<-seq(pos.min,pos.max,by=size*1/2)
          window.bed<-rbind(window.bed,cbind(chr,pos.tag,pos.tag+size))
      }
      window.bed<-window.bed[order(as.numeric(window.bed[,2])),]

      return(list(G=G,chr=chr,pos=pos,window.bed=window.bed))
    }else{
      return(list(G=G,chr=chr,pos=pos,window.bed=NULL))
    }
}

extract_position_universal <- function(col_names) {
  positions <- sapply(col_names, function(col) {
    numbers <- regmatches(col, gregexpr("\\d+", col))[[1]]
    
    if (length(numbers) == 0) {
      return(NA)  
    } else if (length(numbers) == 1) {
      return(as.numeric(numbers))
    } else {
      candidate <- numbers[nchar(numbers) >= 3 & nchar(numbers) <= 9]
      if (length(candidate) > 0) {
        return(as.numeric(tail(candidate, 1)))  
      }
      return(as.numeric(max(numbers)))
    }
  })
  
  return(as.numeric(positions))
}
