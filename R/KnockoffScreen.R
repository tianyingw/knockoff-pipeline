.get_p_in_chunks <- function(X, result.prelim, column_index = NULL,
                             chunk_size = 1000L) {
  if (is.null(column_index)) column_index <- seq_len(ncol(X))
  column_index <- as.integer(column_index)
  if (length(column_index) == 0L) return(matrix(numeric(0), ncol = 1L))
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      is.na(chunk_size) || chunk_size < 1L)
    stop("'chunk_size' must be one positive integer.")

  groups <- split(
    seq_along(column_index),
    ceiling(seq_along(column_index) / as.integer(chunk_size))
  )
  answer <- NULL
  for (target in groups) {
    source <- column_index[target]
    chunk <- as.matrix(Get.p(X[, source, drop = FALSE], result.prelim))
    if (nrow(chunk) != length(target))
      stop("Get.p returned an unexpected number of rows.")
    if (is.null(answer)) {
      answer <- matrix(
        NA_real_, nrow = length(column_index), ncol = ncol(chunk)
      )
      if (!is.null(colnames(chunk))) colnames(answer) <- colnames(chunk)
    } else if (ncol(chunk) != ncol(answer)) {
      stop("Get.p returned an inconsistent number of columns across chunks.")
    }
    answer[target, ] <- chunk
  }
  answer
}


KS.chr<-function(result.prelim,input.X,window.bed,beta=NULL,input.G_k=NULL,region.pos=NULL,tested.pos=NULL,excluded.pos=NULL,M=5,thres.single=0.01,thres.ultrarare=25,thres.missing=0.10,midout.dir=NULL,temp.dir=NULL,jobtitle=NULL,Gsub.id=NULL,impute.method='fixed',bigmemory=T,leveraging=T,LD.filter=NULL,prevalidated=FALSE){
  #region.step=0.1*10^5
  chr<-window.bed[1,1]
  if(length(region.pos)==0){
    region.pos=c(min(window.bed[,2:3]),max(window.bed[,2:3]))
  }

  result.summary.single<-c();result.summary.window<-c();variant.info<-c()
  start.index<-1
  while (start.index<length(region.pos)){
    start<-region.pos[start.index]
    end<-region.pos[start.index+1]-1

    start.index<-start.index+1

    # print('processing genotype data')

    if(length(Gsub.id)==0){match.index<-match(result.prelim$id,rownames(input.X))}else{
      match.index<-match(result.prelim$id,Gsub.id)
    }
    if(mean(is.na(match.index))>0){
      msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
      warning(msg,call.=F)
      match.index <- match.index[!is.na(match.index)]
    }
    if (identical(match.index, seq_len(nrow(input.X)))) {
      G <- input.X
    } else {
      G <- input.X[match.index,,drop=FALSE]
    }
    #sparse matrix operation
    MAF<-.kp_col_means(G)/2;MAC<-.kp_col_sums(G)
    flip_to_minor <- MAF > 0.5 & !is.na(MAF)
    MAC[flip_to_minor] <- nrow(G)*2-MAC[flip_to_minor]
    MAF[flip_to_minor] <- 1-MAF[flip_to_minor]
    if (isTRUE(prevalidated)) {
      SNP.index <- seq_len(ncol(G))
    } else {
      s <- .kp_col_means(G^2) - .kp_col_means(G)^2
      SNP.index<-which(MAF>0 & MAC>=thres.ultrarare & s!=0 & !is.na(MAF))# & MAC>10
    }

    #check.index<-which(MAF>0 & MAC>=thres.ultrarare & s!=0 & !is.na(MAF)  & MISS.freq<0.1)
    if(length(SNP.index)<=1 ){
      msg<-'Number of variants with missing rate <=10% in the specified range is <=1'
      warning(msg,call.=F)
      next
    }
    if (!identical(SNP.index, seq_len(ncol(G)))) {
      G<-G[,SNP.index,drop=FALSE]#;beta<-beta[SNP.index,]
    }
    MAF <- MAF[SNP.index]
    MAC <- MAC[SNP.index]
    # pos<-as.numeric(gsub("^X8\\.","",colnames(G)))
    pos <- as.numeric(colnames(G))
    if(length(beta)==0){beta<-rep(0,ncol(G))}

    ##single variant test for all variants
    p.single <- .get_p_in_chunks(G, result.prelim)

    #output info for tested variants
    temp.variant.info<-cbind(pos,p.single,MAF,MAC)
    colnames(temp.variant.info)<-c('pos','pvalue','MAF','MAC')

    variant.info<-rbind(variant.info,temp.variant.info)
    
    # print(sum(MAF>0.01))
    if (!inherits(G, "Matrix")) G<-Matrix::Matrix(G,sparse=TRUE)

    #get positions
    # pos<-as.numeric(gsub("^X8\\.","",colnames(G)))
    pos <- as.numeric(colnames(G))

    # tic('generating knockoffs')
    #generate knockoffs
    if(is.null(input.G_k)){
      gc()
      if(leveraging==T){
        G_k<-create.KS(G,pos,M=M,corr_max=0.75,maxN.neighbor=Inf,maxBP.neighbor=0.1*10^6,thres.ultrarare=thres.ultrarare,bigmemory=bigmemory,R2.thres=0.75)
      }else{
        G_k<-create.KS(G,pos,M=M,corr_max=0.75,maxN.neighbor=Inf,maxBP.neighbor=0.1*10^6,thres.ultrarare=thres.ultrarare,bigmemory=bigmemory,n.AL=nrow(G),R2.thres=0.75)
      }
      G_k_column_index <- seq_len(ncol(G))
    }else{
      G_k <- input.G_k
      if (length(G_k) != M)
        stop("The number of supplied knockoff matrices does not match M.")
      if (any(vapply(G_k, nrow, numeric(1)) != nrow(input.X)))
        stop("A supplied knockoff matrix has the wrong number of rows.")
      if (any(vapply(G_k, ncol, numeric(1)) < max(SNP.index)))
        stop("A supplied knockoff matrix has too few columns for the QC index.")
      G_k_column_index <- SNP.index
    }
    # toc()
    
    # tic('knockoff analysis')
    ##single variant test for all variants
    p.single_k <- matrix(NA_real_, nrow = ncol(G), ncol = M)
    for(k in 1:M){
      temp.p <- .get_p_in_chunks(
        G_k[[k]], result.prelim, column_index = G_k_column_index
      )
      if (ncol(temp.p) != 1L)
        stop("Single-variant Get.p must return one p-value column.")
      p.single_k[, k] <- temp.p[, 1L]
    }
    # toc()

    # tic('common variants')
    ##common variants
    common.index<-which(MAF>=thres.single)
    # print(length(common.index))
    if(length(common.index)>=1){
      p.common<-p.single[common.index,,drop=F]
      p.common_k<-p.single_k[common.index,,drop=F]

      MK.stat<-MK.statistic(-log10(p.common),-log10(p.common_k),method='median')
      W<-MK.stat[,'tau']*(MK.stat[,'kappa']==0)
      W[!is.finite(W)]<-0
      temp.summary.single<-cbind(chr,pos[common.index],pos[common.index],
                                 pos[common.index],pos[common.index],
                                 (beta!=0)[common.index],
                                 MK.stat,W,p.common,
                                 p.common_k,MAF[common.index])

      colnames(temp.summary.single)<-c('chr','start','end','actual_start','actual_end',
                                       'signal',
                                       'kappa','tau',
                                       'W','P_KS',paste0('P_KS_k',1:M),'MAF')
      if(length(midout.dir)!=0){
        write.table(temp.summary.single,paste0(midout.dir,jobtitle,'_single_',chr,':',start,'-',end,'.txt'),sep='\t',row.names=F,col.names=T,quote=F)
      }
      result.summary.single<-rbind(result.summary.single,temp.summary.single)
    }
    # toc()

    # tic('rare variants')
    ##rare variants
    rare.index<-which(MAF<thres.single & MAC>=thres.ultrarare)
    # print(length(rare.index))
    if(length(rare.index)>=1){
      MAF<-MAF[rare.index];MAC<-MAC[rare.index]
      # pos<-as.numeric(gsub("^X8\\.","",colnames(G)[rare.index]))
      pos <- as.numeric(colnames(G)[rare.index])
      p.rare<-p.single[rare.index]
      p.rare_k<-p.single_k[rare.index,,drop=F]
      beta.rare<-beta[rare.index]
      # print(length(beta))
      #set beta weights and windows
      weight.beta<-dbeta(MAF,1,25)
      weight.matrix<-as.matrix(weight.beta)
      colnames(weight.matrix)<-c(paste0('MAF<',thres.single,'&MAC>',thres.ultrarare,'&Beta'))
      window.temp<-window.bed[window.bed[,2]<max(pos) & window.bed[,3]>min(pos),]
      if(length(nrow(window.temp))==0){window.temp<-matrix(window.temp,1,ncol(window.bed))}
      window.matrix0<-matrix(apply(window.temp,1,function(x)as.numeric(pos>=x[2] & pos<x[3])),length(pos),nrow(window.temp))
      window.string<-apply(window.matrix0,2,function(x)paste(as.character(x),collapse = ""))

      window.MAC<-apply(MAC*window.matrix0,2,sum)
      window.index<-intersect(match(unique(window.string),window.string),which(apply(window.matrix0,2,sum)>1 & window.MAC>=10))
      # window.index<-intersect(match(unique(window.string),window.string),which(apply(window.matrix0,2,sum)>1))
      if(length(window.index)==0){
        start<-end
        next
      }
      window.matrix0<-as.matrix(window.matrix0[,window.index])
      window.matrix<-Matrix(window.matrix0)
      window.summary<-cbind(window.temp[window.index,2],window.temp[window.index,3],t(apply(window.matrix,2,function(x)c(min(pos[which(x==1)]),max(pos[which(x==1)])))))
      p.KS<-c();p.KS_k<-c();p.individual<-c()
      for(i in 1:ceiling(ncol(window.matrix)/100)){
        temp.window.matrix<-window.matrix[,(1+(i-1)*100):min(ncol(window.matrix),i*100),drop=F]
        temp.index<-which(rowSums(temp.window.matrix)!=0)
        temp.weight.matrix<-weight.matrix[temp.index,,drop=F]
        temp.window.matrix<-temp.window.matrix[temp.index,,drop=F]

        temp.G<-G[,rare.index[temp.index],drop=F]
        temp.G_k<-list()
        for(k in 1:M){
          temp.G_k[[k]]<-Matrix(
            G_k[[k]][,G_k_column_index[rare.index[temp.index]],drop=F],
            sparse=T
          )
        }
        KS.fit<-KS.test(temp.G,temp.G_k,p.rare[temp.index],p.rare_k[temp.index,,drop=F],result.prelim,window.matrix=temp.window.matrix,weight.matrix=temp.weight.matrix)
        p.KS<-rbind(p.KS,KS.fit$p.KS)
        p.KS_k<-rbind(p.KS_k,KS.fit$p.KS_k)
        p.individual<-rbind(p.individual,KS.fit$p.individual)
      }
      #######
      p.A<-p.KS;p.A_k<-p.KS_k
      #Knockoff statistics
      MK.stat<-MK.statistic(-log10(p.A),-log10(p.A_k),method='median')
      W<-MK.stat[,'tau']*(MK.stat[,'kappa']==0)
      W[!is.finite(W)]<-0

      temp.summary.window<-cbind(chr,window.summary,
                                 c(t(beta.rare!=0)%*%window.matrix0!=0),
                                 MK.stat,
                                 W,p.A,p.A_k,
                                 p.individual)
      colnames(temp.summary.window)[1:(grep('burden_MAF',colnames(temp.summary.window))-1)]<-c('chr','start','end','actual_start','actual_end',
                                                                                               'signal',
                                                                                               'kappa','tau',
                                                                                               'W','P_KS',paste0('P_KS_k',1:M))
      if(length(midout.dir)!=0){
        write.table(temp.summary.window,paste0(midout.dir,jobtitle,'_window_',chr,':',start,'-',end,'.txt'),sep='\t',row.names=F,col.names=T,quote=F)
      }

      result.summary.window<-rbind(result.summary.window,temp.summary.window)
    #   toc()
    }

    #start<-end
  }
  # rownames(result.summary.window) <- paste(result.summary.window[, "actual_start"], result.summary.window[, "actual_end"], sep = "-")
  return(list(result.single=result.summary.single,result.window=result.summary.window,variant.info=variant.info))
}

create.KS<- function(X,pos,M=5,corr_max=0.75,maxN.neighbor=Inf,maxBP.neighbor=100000,n.AL=floor(10*nrow(X)^(1/3)*log(nrow(X))),thres.ultrarare=25,R2.thres=1,method='shrinkage',bigmemory=T,backing_path=NULL,backing_prefix=NULL,cor.X.precomputed=NULL,preclustered=FALSE) {

  if(class(X)[1]!='dgCMatrix'){X<-Matrix(X,sparse=T)} #convert it to sparse matrix format

  if (is.null(cor.X.precomputed)) {
    cor.X <- .kp_sparse_cov_cor(
      X, need_cov = FALSE, need_cor = TRUE
    )$cor
  } else {
    cor.X <- as.matrix(cor.X.precomputed)
    if (!identical(dim(cor.X), c(ncol(X), ncol(X))) || any(!is.finite(cor.X)))
      stop("The precomputed correlation matrix does not match X.")
  }

  #svd to get leverage score, can be optimized;update: tried fast leveraging, but the R matrix is singular possibly because X is sparse.
  if(method=='svd.irlba'){
    svd.X.u<-irlba(X,nv=floor(sqrt(ncol(X)*log(ncol(X)))))$u
    h<-rowSums(svd.X.u^2)
    prob<-h/sum(h)
  }
  if(method=='svd.full'){
    svd.X.u<-svd(X)$u
    h<-rowSums(svd.X.u^2)
    prob<-h/sum(h)
  }
  if(method=='uniform'){
    h<-rep(1,nrow(X))
    prob<-h/sum(h)
    svd.X.u<-c()
  }
  if(method=='shrinkage'){
    svd.X.u<-irlba(X,nv=floor(sqrt(ncol(X)*log(ncol(X)))))$u
    h1<-rowSums(svd.X.u^2)
    h2<-rep(1,nrow(X))
    prob1<-h1/sum(h1)
    prob2<-h2/sum(h2)
    prob<-0.5*prob1+0.5*prob2
  }
  index.AL<-sample(1:nrow(X),min(n.AL,nrow(X)),replace = FALSE,prob=prob)
  w<-1/sqrt(n.AL*prob[index.AL])
  rm(svd.X.u) #remove temp file

  X.AL<-w*X[index.AL,]
  cov.X.AL <- .kp_sparse_cov_cor(
    X.AL, need_cov = TRUE, need_cor = FALSE
  )$cov
  skip.index <- which(colSums(X.AL != 0) <= thres.ultrarare)

  if (isTRUE(preclustered)) {
    clusters <- seq_len(ncol(X))
  } else if(ncol(X)>1){
    Sigma.distance = as.dist(1 - abs(cor.X))
    fit = hclust(Sigma.distance, method="single")
    clusters = cutree(fit, h=1-corr_max)
  }else{clusters<-1}

  X_k<-list()
  for(k in 1:M){
    if(bigmemory==T && !is.null(backing_path)){
      if (is.null(backing_prefix) || !nzchar(backing_prefix))
        stop("backing_prefix is required when backing_path is supplied.")
      dir.create(backing_path, recursive = TRUE, showWarnings = FALSE)
      backing_file <- paste0(backing_prefix, "_", k, ".bin")
      descriptor_file <- paste0(backing_prefix, "_", k, ".desc")
      if (any(file.exists(file.path(
        backing_path, c(backing_file, descriptor_file)
      )))) stop("Refusing to overwrite existing file-backed knockoff files.")
      X_k[[k]] <- bigmemory::filebacked.big.matrix(
        nrow = nrow(X), ncol = ncol(X), init = 0,
        backingfile = backing_file, descriptorfile = descriptor_file,
        backingpath = backing_path
      )
    }else if(bigmemory==T){X_k[[k]]<-big.matrix(nrow=nrow(X),ncol=ncol(X),init=0)}else{
      X_k[[k]]<-matrix(0,nrow=nrow(X),ncol=ncol(X))
    }
  }

  index.exist<-c()
  for (k in unique(clusters)){
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)){
      #print(i)
      rate<-1;R2<-1;temp.maxN.neighbor<-maxN.neighbor

      while(R2>=R2.thres){

        temp.maxN.neighbor<-floor(temp.maxN.neighbor/rate)

        index.pos<-which(pos>=max(pos[i]-maxBP.neighbor,pos[1]) & pos<=min(pos[i]+maxBP.neighbor,pos[length(pos)]))
        temp<-abs(cor.X[i,]);temp[which(clusters==k)]<-0;temp[-index.pos]<-0
        temp[which(temp<=0.05)]<-0

        index<-order(temp,decreasing=T)
        if(sum(temp!=0,na.rm=T)==0 | temp.maxN.neighbor==0){index<-NULL}else{
          index<-setdiff(index[1:min(length(index),floor((nrow(X))^(1/3)),temp.maxN.neighbor,sum(temp!=0,na.rm=T))],i)
        }
        #index<-setdiff(index[1:min(length(index),sum(temp>corr_base),floor((nrow(X.AL))^(1/3)),maxN.neighbor)],i)

        y<-X[,i]
        if(length(index)==0){fitted.values<-0}
        if(i %in% skip.index){fitted.values<-0}
        if(!(i %in% skip.index |length(index)==0)){

          x.AL<-X.AL[,index,drop=F];
          n.exist<-length(intersect(index,index.exist))
          x.exist.AL<-matrix(0,nrow=nrow(X.AL),ncol=n.exist*M)
          if(length(intersect(index,index.exist))!=0){
            for(j in 1:M){ # this is the most time-consuming part
              x.exist.AL[,((j-1)*n.exist+1):(j*n.exist)]<-w*X_k[[j]][index.AL,intersect(index,index.exist),drop=F]
              #x.exist[,((j-1)*n.exist+1):(j*n.exist)]<-X_k[j,,intersect(index,index.exist),drop=F]
            }
          }
          y.AL<-w*X[index.AL,i];#x.exist.AL<-w*x.exist[index.AL,,drop=F];x.AL<-w*x[index.AL,,drop=F] #create re-scaled data

          temp.xy<-rbind(mean(y.AL),crossprod(x.AL,y.AL)/length(y.AL)-colMeans(x.AL)*mean(y.AL))
          temp.xy<-rbind(temp.xy,crossprod(x.exist.AL,y.AL)/length(y.AL)-colMeans(x.exist.AL)*mean(y.AL))

          temp.cov.cross <- .kp_sparse_cross_cov(x.AL, x.exist.AL)
          temp.cov <- .kp_sparse_cov_cor(
            x.exist.AL, need_cov = TRUE, need_cor = FALSE
          )$cov
          temp.xx<-cov.X.AL[index,index]
          temp.xx<-rbind(cbind(temp.xx,temp.cov.cross),cbind(t(temp.cov.cross),temp.cov))

          temp.xx<-cbind(0,temp.xx)
          temp.xx<-rbind(c(1,rep(0,ncol(temp.xx)-1)),temp.xx)

          svd.fit<-svd(temp.xx)
          v<-svd.fit$v
          cump<-cumsum(svd.fit$d)/sum(svd.fit$d)
          n.svd<-which(cump>=0.999)[1]
          svd.index<-intersect(1:n.svd,which(svd.fit$d!=0))
          temp.inv<-v[,svd.index,drop=F]%*%(svd.fit$d[svd.index]^(-1)*t(v[,svd.index,drop=F]))
          temp.beta<-temp.inv%*%temp.xy
          x<-X[,index,drop=F]
          temp.j<-1
          fitted.values<-temp.beta[1]+x%*%temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F]-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])

          if(length(intersect(index,index.exist))!=0){
            temp.j<-temp.j+ncol(x)
            for(j in 1:M){
              temp.x<-X_k[[j]][,intersect(index,index.exist),drop=F]
              if(ncol(temp.x)>=1){
                fitted.values<-fitted.values+temp.x%*%temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F]-sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])
              }
              temp.j<-temp.j+ncol(temp.x)
            }
          }

        }
        residuals<-as.numeric(y-fitted.values)
        #overfitted model
        R2<-1-var(residuals,na.rm=T)/var(y,na.rm=T)
        rate<-rate*2;temp.maxN.neighbor<-length(index)
      }

      #print(R2)

      #if(var(residuals,na.rm=T)/var(y,na.rm=T)<=0.05){fitted.values<-y;residuals<-y-fitted.values}
      #if(var(residuals,na.rm=T)/var(y,na.rm=T)>=0.95){fitted.values<-0;residuals<-y-fitted.values}
      #print(1-var(residuals,na.rm=T)/var(y,na.rm=T))

      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)

      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[[j]][,which(clusters==k)]<-round(cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F],digits=1)
      #X_k[j,,which(clusters==k)]<-cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F]
    }
  }
  return(X_k)
}

KS_summary<-function(result.window,result.single,M,fdr=0.1){

  temp<-result.single[,match(colnames(result.window),colnames(result.single)),drop=FALSE]
  colnames(temp)<-colnames(result.window)

  result<-rbind(result.window,temp)
  result<-result[order(result[,2]),,drop=FALSE]
  result<-result[order(result[,1]),,drop=FALSE]

  q<-MK.q.byStat(result[,'kappa'],result[,'tau'],M=M)
  threshold<-MK.threshold.byStat(result[,'kappa'],result[,'tau'],M=M,fdr=fdr,Rej.Bound=10000)
  result.summary <- data.frame(
    result[, 1:5, drop = FALSE],  # drop=FALSE 保持数据框结构
    Qvalue = q,
    W_Threshold = threshold,
    result[, -(1:5), drop = FALSE],
    indicator = q <= fdr,
    stringsAsFactors = FALSE
  )
  return(result.summary)
}



KS.test<-function(temp.G,temp.G_k,p.rare,p.rare_k,result.prelim,window.matrix,weight.matrix){
  mu<-result.prelim$nullglm$fitted.values;
  Y.res<-result.prelim$Y-mu;re.Y.res<-result.prelim$re.Y.res
  X0<-result.prelim$X0;outcome<-result.prelim$out_type
  M<-length(temp.G_k)

  #Burden test
  p.burden<-matrix(NA,ncol(window.matrix),ncol(weight.matrix))
  p.burden_k<-array(NA,dim=c(length(temp.G_k),ncol(window.matrix),ncol(weight.matrix)))
  for (k in 1:ncol(weight.matrix)){
    temp.window.matrix<-weight.matrix[,k]*window.matrix
    p.burden[,k]<-Get.p(temp.G%*%temp.window.matrix,result.prelim)
    temp <- matrix(NA_real_, nrow = ncol(window.matrix), ncol = M)
    for(i in 1:M){
      temp[, i] <- as.numeric(
        Get.p(temp.G_k[[i]][]%*%temp.window.matrix,result.prelim)
      )
    }
    p.burden_k[,,k]<-t(temp)
  }

  #Dispersion test
  p.dispersion<-matrix(NA,ncol(window.matrix),ncol(weight.matrix))
  p.dispersion_k<-array(NA,dim=c(length(temp.G_k),ncol(window.matrix),ncol(weight.matrix)))
  
  if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(temp.G))}
  A<-crossprod(temp.G,v*temp.G)#t(temp.G)%*%(v*temp.G)
  B<-crossprod(temp.G,v*X0)#t(temp.G)%*%(v*X0)
  C<-result.prelim$inv.vX0#solve(t(X0)%*%(v*X0))
  K<-A-B%*%C%*%t(B) #here we use the same K for original and knockoffs, due to the exchangeability

  score<-t(temp.G)%*%Y.res#;re.score<-t(t(temp.G)%*%re.Y.res)
  score_k <- matrix(NA_real_, nrow = ncol(temp.G), ncol = M)
  for(i in 1:M){
    score_k[, i] <- as.numeric(crossprod(temp.G_k[[i]][,],Y.res))
  }
  for (k in 1:ncol(weight.matrix)){
    # print(k)
    skat_prepared <- .prepare_skat_ks(
      K, window.matrix, weight = weight.matrix[, k]
    )
    p.dispersion[,k] <- .get_p_skat_ks_prepared(
      score, window.matrix, weight.matrix[, k], skat_prepared
    )
    p.dispersion_k[,,k] <- t(vapply(
      seq_len(M),
      function(s) as.numeric(.get_p_skat_ks_prepared(
        score_k[, s], window.matrix, weight.matrix[, k], skat_prepared
      )),
      numeric(ncol(window.matrix))
    ))
  }
  p.V1<-Get.cauchy.scan(p.rare,window.matrix)
  p.V1_k<-apply(p.rare_k,2,Get.cauchy.scan,window.matrix=window.matrix)
  if(ncol(window.matrix)==1){p.V1_k<-matrix(p.V1_k,1,length(temp.G_k))}

  p.individual<-cbind(p.burden,p.dispersion,p.V1);
  colnames(p.individual)<-c(paste0('burden_',colnames(weight.matrix)),
                            paste0('dispersion_',colnames(weight.matrix)),
                            paste0('singleCauchy'))

  p.KS<-as.matrix(apply(p.individual,1,Get.cauchy))
  p.KS_k<-matrix(sapply(1:length(temp.G_k),function(s){apply(cbind(matrix(p.burden_k[s,,],dim(p.burden_k)[2],dim(p.burden_k)[3]),
                                                                   matrix(p.dispersion_k[s,,],dim(p.dispersion_k)[2],dim(p.dispersion_k)[3]),
                                                                   p.V1_k[,s]),1,Get.cauchy)}),dim(p.burden_k)[2],dim(p.burden_k)[1])

  return(list(p.KS=p.KS,p.KS_k=p.KS_k,p.individual=p.individual))
}
Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, log.p=FALSE){
  param<-Get_Liu_Params_Mod_Lambda(lambda)
  Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
  Q.Norm1<-Q.Norm * param$sigmaX + param$muX
  p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
  return(p.value)
}

Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation

  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }

  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

.prepare_skat_ks <- function(K, window.matrix, weight) {
  K.temp <- weight * t(weight * K)
  lapply(seq_len(ncol(window.matrix)), function(i) {
    member <- window.matrix[, i] != 0
    temp <- K.temp[member, member, drop = FALSE]
    if (sum(temp^2) == 0) return(NULL)
    lambda <- eigen(temp, symmetric = TRUE, only.values = TRUE)$values
    if (anyNA(lambda)) return(NULL)
    lambda
  })
}

.get_p_skat_ks_prepared <- function(score, window.matrix, weight,
                                    prepared) {
  Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2)
  p<-rep(NA,length(Q))
  for(i in 1:length(Q)){
    lambda <- prepared[[i]]
    if (is.null(lambda)) {p[i] <- NA; next}

    #temp.p<-SKAT_davies(Q[i],lambda,acc=10^(-6))$Qq
    temp.p<-davies(Q[i],lambda,acc=10^(-6))$Qq

    if(temp.p > 1 || temp.p <= 0 ){
      temp.p<-Get_Liu_PVal.MOD.Lambda(Q[i],lambda)
    }
    p[i]<-temp.p
  }

  return(as.matrix(p))
}

Get.p.SKAT.KS<-function(score,K,window.matrix,weight,result.prelim){
  prepared <- .prepare_skat_ks(K, window.matrix, weight)
  .get_p_skat_ks_prepared(score, window.matrix, weight, prepared)
}



max_nth<-function(x,n){return(sort(x,partial=length(x)-(n-1))[length(x)-(n-1)])}

Get.p.base<-function(X,result.prelim){
  #X<-Matrix(X)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  outcome<-result.prelim$out_type
  if(outcome=='D'){
    v<-mu*(1-mu)
    A<-(t(X)%*%Y.res)^2
    B<-colSums(v*X^2)
    C<-t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0
    D<-t(t(result.prelim$X0)%*%as.matrix(v*X))
    p<-pchisq(as.numeric(A/(B-rowSums(C*D))),df=1,lower.tail=F)
  }else{
    p <- .kp_continuous_score_p(X, result.prelim)
  }
  #p<-pchisq(as.numeric((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum))),df=1,lower.tail=F)
  #p[is.na(p)]<-NA
  return(as.matrix(p))
}

Get.p<-function(X,result.prelim){
  #X<-as.matrix(X)
  outcome<-result.prelim$out_type
  if(outcome=='D'){
    invisible(capture.output(
      p <- WGScan::ScoreTest_SPA(
        t(X), result.prelim$Y, result.prelim$X0,
        method = c("fastSPA"), minmac = -Inf
      )$p.value
    ))
  }else{
    p <- .kp_continuous_score_p(X, result.prelim)
    #p<-pchisq(as.numeric((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum))),df=1,lower.tail=F)
  }
  return(as.matrix(p))
}

Get.Z<-function(X,result.prelim){
  #X<-Matrix(X)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  sd.X<-apply(X,2,sd)
  Z<-t(apply(X,2,scale))%*%as.matrix(scale(Y.res))/sqrt(length(Y.res))
  Z[sd.X==0,]<-0
  return(as.matrix(Z))
}

MK.statistic<-function (T_0,T_k,method='median'){
  if (length(method) != 1L || is.na(method) ||
      !method %in% c("median", "max"))
    stop("method must be 'median' or 'max'.", call. = FALSE)
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  T.temp<-cbind(T_0,T_k)
  invalid <- apply(is.na(T.temp) | is.nan(T.temp) |
                     (is.infinite(T.temp) & T.temp < 0), 1, any)
  # Association-test p-values can underflow to zero.  Retain their ordering
  # without allowing an infinite W to collide with Inf's use as the
  # no-rejection threshold sentinel.
  T.temp[is.infinite(T.temp) & T.temp > 0] <-
    -log10(.Machine$double.xmin)
  if (any(invalid)) T.temp[invalid, ] <- 0

  which.max.alt<-function(x){
    temp.index<-which(x==max(x))
    if(length(temp.index)!=1){return(temp.index[2])}else{return(temp.index[1])}
  }
  kappa<-apply(T.temp,1,which.max.alt)-1

  if(method=='max'){tau<-apply(T.temp,1,max)-apply(T.temp,1,max_nth,n=2)}
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  kappa[invalid] <- NA_integer_
  tau[invalid] <- 0
  return(cbind(kappa,tau))
}

.MK.fdp.path <- function(kappa, tau, M, Rej.Bound = 10000) {
  if (length(kappa) != length(tau))
    stop("kappa and tau must have the same length.", call. = FALSE)
  if (length(M) != 1L || is.na(M) || !is.finite(M) ||
      M < 1 || M > .Machine$integer.max || M != floor(M))
    stop("M must be a positive integer.", call. = FALSE)
  if (length(Rej.Bound) != 1L || is.na(Rej.Bound) || Rej.Bound <= 0 ||
      (is.finite(Rej.Bound) &&
       (Rej.Bound > .Machine$integer.max || Rej.Bound != floor(Rej.Bound))))
    stop("Rej.Bound must be a positive integer or Inf.", call. = FALSE)
  M <- as.integer(M)
  observed_kappa <- kappa[!is.na(kappa)]
  if (any(!is.finite(observed_kappa)) ||
      any(observed_kappa < 0 | observed_kappa > M |
          observed_kappa != floor(observed_kappa)))
    stop("Non-missing kappa values must be integers between 0 and M.",
         call. = FALSE)
  if (any(is.infinite(tau) & tau < 0, na.rm = TRUE))
    stop("Non-missing tau values cannot be negative infinity.",
         call. = FALSE)
  # Positive infinity can occur in older intermediate files when a p-value
  # underflowed to zero.  Use the same finite cap as MK.statistic so that Inf
  # remains reserved for the no-rejection threshold sentinel.
  tau[is.infinite(tau) & tau > 0] <- -log10(.Machine$double.xmin)
  eligible <- which(!is.na(kappa) & !is.na(tau) & tau > 0)
  if (length(eligible) == 0L)
    return(list(index = integer(), group = integer(), threshold = numeric(),
                fdp = numeric()))

  b <- eligible[order(tau[eligible], decreasing = TRUE)]
  if (is.finite(Rej.Bound) && length(b) > Rej.Bound) {
    boundary <- max(1L, as.integer(Rej.Bound))
    boundary_tau <- tau[b[boundary]]
    boundary <- max(which(tau[b] == boundary_tau))
    b <- b[seq_len(boundary)]
  }

  original <- kappa[b] == 0
  n_original <- cumsum(original)
  n_knockoff <- seq_along(b) - n_original
  group_end <- which(!duplicated(tau[b], fromLast = TRUE))
  fdp <- (1 / M + n_knockoff[group_end] / M) /
    pmax(1, n_original[group_end])
  list(
    index = b,
    group = match(tau[b], tau[b][group_end]),
    threshold = tau[b][group_end],
    fdp = fdp
  )
}


MK.threshold.byStat<-function (kappa,tau,M,fdr = 0.1,Rej.Bound=10000){
  if (length(fdr) != 1L || is.na(fdr) || !is.finite(fdr) ||
      fdr < 0 || fdr > 1)
    stop("fdr must be a finite number between 0 and 1.", call. = FALSE)
  path <- .MK.fdp.path(kappa, tau, M, Rej.Bound)
  ok <- which(path$fdp <= fdr)
  if (length(ok) > 0L) path$threshold[ok[length(ok)]] else Inf
}

MK.threshold<-function (T_0,T_k, fdr = 0.1,method='median',Rej.Bound=10000){
  stat<-MK.statistic(T_0,T_k,method=method)
  kappa<-stat[,1];tau<-stat[,2]
  t<-MK.threshold.byStat(kappa,tau,M=ncol(T_k),fdr=fdr,Rej.Bound=Rej.Bound)
  return(t)
}


MK.q.byStat<-function (kappa,tau,M,Rej.Bound=10000){
  path <- .MK.fdp.path(kappa, tau, M, Rej.Bound)
  q <- rep(1, length(tau))
  if (length(path$index) == 0L) return(q)

  group_q <- rev(cummin(rev(path$fdp)))
  original <- kappa[path$index] == 0
  q[path$index[original]] <- pmin(
    1, group_q[path$group[original]]
  )
  q
}



Get.cauchy<-function(p){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))

  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}

Get.cauchy.scan<-function(p,window.matrix){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16) & !is.na(p)
  temp<-rep(0,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[!is.small]<-as.numeric(tan((0.5-p[!is.small])*pi))
  #window.matrix.MAC10<-(MAC>=10)*window.matrix0

  cct.stat<-as.numeric(t(temp)%*%window.matrix/apply(window.matrix,2,sum))
  #cct.stat<-as.numeric(t(temp)%*%window.matrix.MAC10/apply(window.matrix.MAC10,2,sum))
  is.large<-cct.stat>1e+15 & !is.na(cct.stat)
  is.regular<-cct.stat<=1e+15 & !is.na(cct.stat)
  pval<-rep(NA,length(cct.stat))
  pval[is.large]<-(1/cct.stat[is.large])/pi
  pval[is.regular]<-1-pcauchy(cct.stat[is.regular])
  return(pval)
}

Get.p.moment<-function(Q,re.Q){ #Q a A*q matrix of test statistics, re.Q a B*q matrix of resampled test statistics
  re.mean<-apply(re.Q,2,mean)
  re.variance<-apply(re.Q,2,var)
  re.kurtosis<-apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
  re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
  re.p<-t(pchisq((t(Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df,lower.tail=F))
  #re.p[re.p==1]<-0.99
  return(re.p)
}
