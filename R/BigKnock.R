utils::globalVariables(c('create.MK.AL_gene_buffer','G_gene_buffer_surround','LD.filter',
                         'surround.region','G_gene_buffer','G_EnhancerAll','p_EnhancerAll',
                         'pos_gene_buffer','G_Enhancer','n','G_enhancer_surround','pos_enhancer'))

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
                                                         enhancer_start=as.numeric(Enhancer.pos[r,1]),enhancer_end=as.numeric(Enhancer.pos[r,2]),M=M)))
         
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

Knockoffgeneration.gene.buffer=function(G_gene_buffer_surround=G_gene_buffer_surround,
                                        gene_buffer_start=gene_buffer_start,
                                        gene_buffer_end=gene_buffer_end,
                                        M=5,surround.region=100000,LD.filter=0.75){

  #missing genotype imputation
  G_gene_buffer_surround[G_gene_buffer_surround<0 | G_gene_buffer_surround>2]<-NA
  N_MISS<-sum(is.na(G_gene_buffer_surround))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G_gene_buffer_surround)/ncol(G_gene_buffer_surround))
    #print(msg,call.=F)
    colmean<-colMeans(x = G_gene_buffer_surround, na.rm = T)
    index <- which(is.na(G_gene_buffer_surround), arr.ind=TRUE)
    G_gene_buffer_surround[index] <- colmean[index[,2]]
  }

  #sparse matrix operation
  MAF<-colMeans(G_gene_buffer_surround)/2;MAC<-colSums(G_gene_buffer_surround)
  MAF[MAF>0.5]<-1-MAF[MAF>0.5]
  MAC[MAF>0.5]<-nrow(G_gene_buffer_surround)*2-MAC[MAF>0.5]
  s<-colMeans(G_gene_buffer_surround^2)-colMeans(G_gene_buffer_surround)^2
  SNP.index<-which(MAF>0 & MAC>=25 & s!=0 & !is.na(MAF))

  if(length(SNP.index)<=1 ){
    msg<-'Number of variants with missing rate <=10% in the specified range is <=1'
    #print(msg,call.=F)
    stop
  }
  G_gene_buffer_surround<-G_gene_buffer_surround[,SNP.index,drop=F]

  #get positions and reorder G_gene_buffer_surround
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_gene_buffer_surround)))
  G_gene_buffer_surround<-G_gene_buffer_surround[,order(pos),drop=F]

  MAF<-colMeans(G_gene_buffer_surround)/2
  G_gene_buffer_surround<-as.matrix(G_gene_buffer_surround)
  G_gene_buffer_surround[,MAF>0.5 & !is.na(MAF)]<-2-G_gene_buffer_surround[,MAF>0.5 & !is.na(MAF)]
  MAF<-colMeans(G_gene_buffer_surround)/2;MAC<-colSums(G_gene_buffer_surround)

  G_gene_buffer_surround<-Matrix(G_gene_buffer_surround,sparse=T)
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_gene_buffer_surround)))
  n=dim(G_gene_buffer_surround)[1]

  max.corr=1
  while(max.corr>=LD.filter){ #max corr < 0.75
    #clustering and filtering
    G_gene_buffer_surround=G_gene_buffer_surround
    sparse.fit<-sparse.cor(G_gene_buffer_surround)
    cor.X<-sparse.fit$cor;cov.X<-sparse.fit$cov
    range(c(cor.X)[round(c(cor.X),digits = 2)!=1.00])
    max.corr=max(abs(c(cor.X)[round(c(cor.X),digits = 2)!=1.00]))

    Sigma.distance = as.dist(1 - abs(cor.X))
    if(ncol(G_gene_buffer_surround)>1){
      fit = hclust(Sigma.distance, method="complete")
      corr_max = 0.75
      clusters = cutree(fit, h=1-corr_max)
    }else{clusters<-1}

    ##apply the LD filter before knockoff generation
    #One variant is randomly selected as the representative per cluster.
    #If a cluster is inside the gene-buffer region, we prioritize to keep one variant inside the gene buffer region instead of outsides
    gene_buffer_ind=(pos>=gene_buffer_start&pos<=gene_buffer_end)

    set.seed(12345)
    temp.index.gene_buffer<-sample(sum(gene_buffer_ind))
    temp.index.gene_buffer<-temp.index.gene_buffer[match(unique(clusters[gene_buffer_ind]),clusters[gene_buffer_ind][temp.index.gene_buffer])]
    if(length(temp.index.gene_buffer)<=1 ){
      msg<-'Number of variants after LD filtering in the gene buffer is <=1'
      warning(msg,call.=F)
      break
    }
    gene_buffer.index=which(gene_buffer_ind)[temp.index.gene_buffer]

    ##Then filter other variants in +-100kb surrounding region
    temp.index.surround<-sample(length(pos)-sum(gene_buffer_ind))
    temp.index.surround<-temp.index.surround[match(unique(clusters[!gene_buffer_ind]),clusters[!gene_buffer_ind][temp.index.surround])]
    surround.index=which(!gene_buffer_ind)[temp.index.surround]
    surround.index=surround.index[!clusters[which(!gene_buffer_ind)[temp.index.surround]]%in%unique(clusters[gene_buffer_ind])]

    temp.index=unique(c(gene_buffer.index,surround.index))

    G_gene_buffer_surround<-G_gene_buffer_surround[,temp.index,drop=F]
    pos=pos[temp.index]
  }

  #print('generating knockoffs of gene buffer region') #knockoff-AL for gene buffer, adapt the code of KnockoffScreen-AL
  set.seed(12345)
  G_gene_buffer_knockoff=create.MK.AL_gene_buffer(X=G_gene_buffer_surround,pos=pos,
                                                  gene_buffer_start=gene_buffer_start,gene_buffer_end=gene_buffer_end,M=M,
                                                  corr_max=LD.filter,maxN.neighbor=Inf,
                                                  maxBP.neighbor=surround.region,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                                  thres.ultrarare=25,R2.thres=LD.filter)

  return(G_gene_buffer_knockoff)
}

Knockoffgeneration.enhancer=function(G_enhancer_surround=G_enhancer_surround,
                                     enhancer_start=enhancer_start,
                                     enhancer_end=enhancer_start,
                                     M=5,surround.region=50000,LD.filter=0.75){

  #missing genotype imputation
  G_enhancer_surround[G_enhancer_surround<0 | G_enhancer_surround>2]<-NA
  N_MISS<-sum(is.na(G_enhancer_surround))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G_enhancer_surround)/ncol(G_enhancer_surround))
    print(msg,call.=F)
    colmean<-colMeans(x = G_enhancer_surround, na.rm = T)
    index <- which(is.na(G_enhancer_surround), arr.ind=TRUE)
    G_enhancer_surround[index] <- colmean[index[,2]]
  }

  #sparse matrix operation
  MAF<-colMeans(G_enhancer_surround)/2;MAC<-colSums(G_enhancer_surround)
  MAF[MAF>0.5]<-1-MAF[MAF>0.5]
  MAC[MAF>0.5]<-nrow(G_enhancer_surround)*2-MAC[MAF>0.5]
  s<-colMeans(G_enhancer_surround^2)-colMeans(G_enhancer_surround)^2
  SNP.index<-which(MAF>0 & MAC>=25 & s!=0 & !is.na(MAF))

  if(length(SNP.index)<=1 ){
    msg<-'Number of variants with missing rate <=10% in the specified range is <=1'
    print(msg,call.=F)
    stop
  }
  G_enhancer_surround<-G_enhancer_surround[,SNP.index,drop=F]

  #get positions and reorder G_enhancer_surround
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_enhancer_surround)))
  G_enhancer_surround<-G_enhancer_surround[,order(pos),drop=F]

  MAF<-colMeans(G_enhancer_surround)/2
  G_enhancer_surround<-as.matrix(G_enhancer_surround)
  G_enhancer_surround[,MAF>0.5 & !is.na(MAF)]<-2-G_enhancer_surround[,MAF>0.5 & !is.na(MAF)]
  MAF<-colMeans(G_enhancer_surround)/2;MAC<-colSums(G_enhancer_surround)

  G_enhancer_surround<-Matrix(G_enhancer_surround,sparse=T)
  pos<-as.numeric(gsub("^.*\\:","",colnames(G_enhancer_surround)))
  n=dim(G_enhancer_surround)[1]

  max.corr=1
  while(max.corr>=LD.filter){ #max corr < 0.75
    #clustering and filtering
    G_enhancer_surround=G_enhancer_surround
    sparse.fit<-sparse.cor(G_enhancer_surround)
    cor.X<-sparse.fit$cor;cov.X<-sparse.fit$cov
    range(c(cor.X)[round(c(cor.X),digits = 2)!=1.00])
    max.corr=max(abs(c(cor.X)[round(c(cor.X),digits = 2)!=1.00]))

    Sigma.distance = as.dist(1 - abs(cor.X))
    if(ncol(G_enhancer_surround)>1){
      fit = hclust(Sigma.distance, method="complete")
      corr_max = 0.75
      clusters = cutree(fit, h=1-corr_max)
    }else{clusters<-1}

    ##apply the LD filter before knockoff generation
    #One variant is randomly selected as the representative per cluster.
    #If a cluster is inside the gene-buffer region, we prioritize to keep one variant inside the gene buffer region instead of outsides
    enhancer_ind=(pos>=enhancer_start&pos<=enhancer_end)

    set.seed(12345)
    temp.index.enhancer<-sample(sum(enhancer_ind))
    temp.index.enhancer<-temp.index.enhancer[match(unique(clusters[enhancer_ind]),clusters[enhancer_ind][temp.index.enhancer])]
    if(length(temp.index.enhancer)<=1 ){
      msg<-'Number of variants after LD filtering in the gene buffer is <=1'
      warning(msg,call.=F)
      break
    }
    enhancer.index=which(enhancer_ind)[temp.index.enhancer]

    ##Then filter other variants in +-100kb surrounding region
    temp.index.surround<-sample(length(pos)-sum(enhancer_ind))
    temp.index.surround<-temp.index.surround[match(unique(clusters[!enhancer_ind]),clusters[!enhancer_ind][temp.index.surround])]
    surround.index=which(!enhancer_ind)[temp.index.surround]
    surround.index=surround.index[!clusters[which(!enhancer_ind)[temp.index.surround]]%in%unique(clusters[enhancer_ind])]

    temp.index=unique(c(enhancer.index,surround.index))

    G_enhancer_surround<-G_enhancer_surround[,temp.index,drop=F]
    pos=pos[temp.index]
  }

#   print('generating knockoffs of enhancer')
  set.seed(12345)
  G_enhancer_knockoff=create.MK.AL_enhancer(X=G_enhancer_surround,pos=pos,
                                            enhancer_start=enhancer_start,enhancer_end=enhancer_end,M=M,
                                            corr_max=LD.filter,maxN.neighbor=Inf,
                                            maxBP.neighbor=surround.region,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                            thres.ultrarare=25,R2.thres=LD.filter)

  return(G_enhancer_knockoff)
}


######### Other functions #########
#Optimize create.MK.AL function provided by Zihuai
#Knockoff generation for gene buffer regions
create.MK.AL_gene_buffer <- function(X=G_gene_buffer_surround,pos,gene_buffer_start,gene_buffer_end,M,corr_max=LD.filter,maxN.neighbor=Inf,
                                     maxBP.neighbor=surround.region,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                     thres.ultrarare=25,R2.thres=LD.filter) {

  method='shrinkage'
  sparse.fit<-sparse.cor(X)
  cor.X<-sparse.fit$cor;cov.X<-sparse.fit$cov  #correlation

  #svd to get leverage score, can be optimized;update: tried fast leveraging, but the R matrix is singular possibly because X is sparse.
  #Fast Truncated Singular Value Decomposition
  if(method=='shrinkage'){
    svd.X.u<-irlba(X,nv=floor(sqrt(ncol(X)*log(ncol(X)))))$u #U is the orthogonal singular vectors
    h1<-rowSums(svd.X.u^2)
    h2<-rep(1,nrow(X))
    prob1<-h1/sum(h1)
    prob2<-h2/sum(h2)
    prob<-0.5*prob1+0.5*prob2 #shrinkage leveraging estimator, probability weights for sampling
  }

  index.AL<-sample(1:nrow(X),min(n.AL,nrow(X)),replace = FALSE,prob=prob) #sampling r samples from n samples, using shrinkage leveraging estimator
  w<-1/sqrt(n.AL*prob[index.AL])
  rm(svd.X.u) #remove temp file

  X.AL<-w*X[index.AL,] #n.AL samples
  sum(is.na(X.AL)) #0

  sparse.fit<-sparse.cor(X.AL)
  cor.X.AL<-sparse.fit$cor;cov.X.AL<-sparse.fit$cov
  skip.index<-colSums(X.AL!=0)<=thres.ultrarare #skip features that are ultra sparse, permutation will be directly applied to generate knockoffs

  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single") #hierarchical clustering
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)  #variants from two different clusters do not have a correlation greater than 0.75.
  }else{clusters<-1}

  gc()
  X_k<-list()
  for(k in 1:M){
    X_k[[k]]<-matrix(0,nrow=nrow(X),ncol=ncol(X))
    #X_k[[k]]<-big.matrix(nrow=nrow(X),ncol=ncol(X),init=0,shared=FALSE)
  }

  ##only run snps within gene buffer
  snps_ind=which(pos<=gene_buffer_end&pos>=gene_buffer_start)

  index.exist<-c()
  for (k in unique(clusters[snps_ind])){
    #print(paste0('cluster',k))
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)[which(clusters==k)%in%snps_ind]){
      #print(i)
      rate<-1;R2<-1;temp.maxN.neighbor<-maxN.neighbor
      while(R2>=R2.thres){ #avoid over-fitting
        temp.maxN.neighbor<-floor(temp.maxN.neighbor/rate)
        snp.pos=as.numeric(gsub("^.*\\:","",names(clusters[i])))
        #+-100kb surrounding region
        index.pos<-which(pos>=max(snp.pos-maxBP.neighbor,pos[1]) & pos<=min(snp.pos+maxBP.neighbor,pos[length(pos)]))
        #correlation between this snp with other snps in +-100kb surrounding region
        temp<-abs(cor.X[i,])
        temp[which(clusters==k)]<-0 #exclude variants if they are in the same cluster as the target variant
        temp[-index.pos]<-0 #only focus on +-100kb surrounding region
        temp[which(temp<=corr_base)]<-0
        index<-order(temp,decreasing=T)
        if(sum(temp!=0,na.rm=T)==0 | temp.maxN.neighbor==0){index<-NULL}else{
          index<-setdiff(index[1:min(length(index),floor((nrow(X))^(1/3)),temp.maxN.neighbor,sum(temp!=0,na.rm=T))],i)
        } #top K snps up to K=n^1/3=75

        y<-X[,i] #n samples
        if(length(index)==0){fitted.values<-0}
        if(i %in% skip.index){fitted.values<-0}
        if(!(i %in% skip.index |length(index)==0)){
          x.AL<-X.AL[,index,drop=F]; #n.AL by K
          n.exist<-length(intersect(index,index.exist))
          x.exist.AL<-matrix(0,nrow=nrow(X.AL),ncol=n.exist*M)
          if(length(intersect(index,index.exist))!=0){
            for(j in 1:M){ # this is the most time-consuming part
              x.exist.AL[,((j-1)*n.exist+1):(j*n.exist)]<-w*X_k[[j]][index.AL,intersect(index,index.exist),drop=F]
            }
          }
          y.AL<-w*X[index.AL,i]; #n.AL

          temp.xy<-rbind(mean(y.AL),crossprod(x.AL,y.AL)/length(y.AL)-colMeans(x.AL)*mean(y.AL))
          temp.xy<-rbind(temp.xy,crossprod(x.exist.AL,y.AL)/length(y.AL)-colMeans(x.exist.AL)*mean(y.AL))
          temp.cov.cross<-sparse.cov.cross(x.AL,x.exist.AL)$cov
          temp.cov<-sparse.cor(x.exist.AL)$cov
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
          temp.beta<-temp.inv%*%temp.xy #least square estimate for regression coefficient, alpha and beta_k

          x<-X[,index,drop=F]
          temp.j<-1
          fitted.values<-temp.beta[1]+x%*%temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F]-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])
          length(fitted.values) #n samples

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
      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)
      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[[j]][,which(clusters==k)]<-round(cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F],digits=1)
    }
  }

  #save knockoffs of gene buffer region
#   print('saving knockoffs of gene buffer region')
  #G_gene_buffer=X[,snps_ind]
  G_gene_buffer_knockoff <- array(0, dim = c(M, nrow(X), length(snps_ind)))
  for (j in 1:M) {
    G_gene_buffer_knockoff[j, ,] <-X_k[[j]][,snps_ind]
  }
  rm(X_k)

  #G_gene_buffer_knockoff=list(G_gene_buffer=G_gene_buffer,G_gene_buffer_knockoff=G_gene_buffer_knockoff)
  return(G_gene_buffer_knockoff)
}

create.MK.AL_enhancer <- function(X=G_enhancer_surround,pos,enhancer_start,enhancer_end,M,corr_max=0.75,maxN.neighbor=Inf,
                                  maxBP.neighbor=50000,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                  thres.ultrarare=25,R2.thres=0.75) {

  method='shrinkage'
  sparse.fit<-sparse.cor(X)
  cor.X<-sparse.fit$cor;cov.X<-sparse.fit$cov

  #svd to get leverage score, can be optimized;update: tried fast leveraging, but the R matrix is singular possibly because X is sparse.
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
  sparse.fit<-sparse.cor(X.AL)
  cor.X.AL<-sparse.fit$cor;cov.X.AL<-sparse.fit$cov
  skip.index<-colSums(X.AL!=0)<=thres.ultrarare #skip features that are ultra sparse, permutation will be directly applied to generate knockoffs

  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single")
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)
  }else{clusters<-1}

  X_k<-list()
  ##only focus on snps within gene buffer
  for(k in 1:M){
    #X_k[[k]]<-big.matrix(nrow=nrow(X),ncol=ncol(X),init=0,shared=FALSE)
    X_k[[k]]<-matrix(0,nrow=nrow(X),ncol=ncol(X))
  }

  snps_ind=which(pos<=enhancer_end&pos>=enhancer_start)

  index.exist<-c()
  for (k in unique(clusters[snps_ind])){
    #print(paste0('cluster',k))
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    for(i in which(clusters==k)[which(clusters==k)%in%snps_ind]){
      #print(i)
      rate<-1;R2<-1;temp.maxN.neighbor<-maxN.neighbor

      while(R2>=R2.thres){

        temp.maxN.neighbor<-floor(temp.maxN.neighbor/rate)
        snp.pos=as.numeric(gsub("^.*\\:","",names(clusters[i])))
        index.pos<-which(pos>=max(snp.pos-maxBP.neighbor,pos[1]) & pos<=min(snp.pos+maxBP.neighbor,pos[length(pos)]))

        temp<-abs(cor.X[i,]);temp[which(clusters==k)]<-0;temp[-index.pos]<-0
        temp[which(temp<=corr_base)]<-0

        index<-order(temp,decreasing=T)
        if(sum(temp!=0,na.rm=T)==0 | temp.maxN.neighbor==0){index<-NULL}else{
          index<-setdiff(index[1:min(length(index),floor((nrow(X))^(1/3)),temp.maxN.neighbor,sum(temp!=0,na.rm=T))],i)
        }

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
            }
          }
          y.AL<-w*X[index.AL,i];

          temp.xy<-rbind(mean(y.AL),crossprod(x.AL,y.AL)/length(y.AL)-colMeans(x.AL)*mean(y.AL))
          temp.xy<-rbind(temp.xy,crossprod(x.exist.AL,y.AL)/length(y.AL)-colMeans(x.exist.AL)*mean(y.AL))
          temp.cov.cross<-sparse.cov.cross(x.AL,x.exist.AL)$cov
          temp.cov<-sparse.cor(x.exist.AL)$cov
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
      cluster.fitted[,match(i,which(clusters==k))]<-as.vector(fitted.values)
      cluster.residuals[,match(i,which(clusters==k))]<-as.vector(residuals)
      index.exist<-c(index.exist,i)
    }
    #sample mutiple knockoffs
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[[j]][,which(clusters==k)]<-round(cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F],digits=1)
    }
  }

  #save knockoffs of enhancer
#   print('saving knockoffs of enhancer')
  G_enhancer_knockoff <- array(0, dim = c(M, nrow(X), length(snps_ind)))
  for (j in 1:M) {
    G_enhancer_knockoff[j, ,] <-X_k[[j]][,snps_ind]
  }
  rm(X_k)

  return(G_enhancer_knockoff)
}



GeneScan3D.UKB.GLMM<-function(G=G_gene_buffer,G.EnhancerAll=G_EnhancerAll,R=length(p_EnhancerAll),
                              p_Enhancer=p_EnhancerAll,window.size=c(1000,5000,10000),pos=pos_gene_buffer,
                              MAC.threshold=10,MAF.threshold=0.01,Gsub.id=Gsub.id,
                              result.null.model.GLMM=result.null.model.GLMM,outcome='C',
                              sparseSigma=sparseSigma,ratio=ratio){
  #load preliminary features
  mu<-as.vector(result.null.model.GLMM$fitted.values)
  Y.res<-as.vector(result.null.model.GLMM$residuals)
  X<-result.null.model.GLMM$X #covariates include intercept
  
  invSigma_X<-solve(sparseSigma, X, sparse=T)
  C<-solve(t(X)%*%invSigma_X)
  #genotype filtering/checking/missing values imputation
  G_filter=Genotype_filter(G,pos,impute.method='fixed')
  G=G_filter$G
  pos=G_filter$pos

  #match phenotype id (phecode) and genotype id
  if(length(Gsub.id)==0){match.index<-match(as.numeric(result.null.model.GLMM$sampleID),1:nrow(G))}else{
    match.index<-match(result.null.model.GLMM$sampleID,Gsub.id)
  }
  if(mean(is.na(match.index))>0){
    msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
    warning(msg,call.=F)
  }
  
  #individuals ids are matched with genotype
  G=Matrix(G[match.index,])
#   print("match")
  #generate window matrix to specify the variants in each window
  window.matrix0_gene_buffer<-c()
  for(size in window.size){
    if (size==1){next}
    pos.tag<-seq(min(pos),max(pos),by=size*1/2)
    pos.tag<-sapply(pos.tag,function(x)pos[which.min(abs(x-pos))])
    window.matrix0_gene_buffer<-cbind(window.matrix0_gene_buffer,sapply(pos.tag,function(x)as.numeric(pos>=x & pos<x+size)))
  }

  window.string_gene_buffer<-apply(window.matrix0_gene_buffer,2,function(x)paste(as.character(x),collapse = ""))
  window.matrix_gene_buffer<-Matrix(window.matrix0_gene_buffer[,match(unique(window.string_gene_buffer),window.string_gene_buffer)])
  #Number of 1-D windows to scan the gene buffer region
  M_gene_buffer=dim(window.matrix_gene_buffer)[2]
#   print("window matrix")
  ##single variant score tests, related samples using SAIGE null GLMM
  if(outcome=='D'){v=as.numeric((mu*(1-mu)))}
  if(outcome=='C'){v=1/result.null.model.GLMM$theta[1]} #phi is residual variance

  #covariate adjusted genotypes
  G_tilde=G-X%*%solve(t(X)%*%(v*X))%*%(t(X)%*%(v*G))

  #variance-adjusted score statistics
  #as.vector(t(G_tilde)%*%Y.res)==as.vector(t(G)%*%Y.res)
  S=as.vector(t(G_tilde)%*%Y.res)/result.null.model.GLMM$theta[1]

  ##GLMM
  #adjusted score statistics, without SPA
  if(outcome=='C'){
    invSigma_G_tilde<-solve(sparseSigma, G_tilde, sparse=T)
    V=t(G_tilde)%*%invSigma_G_tilde
    p.single=pchisq(S^2/(ratio*diag(V)),df=1,lower.tail=F)
  }

  #with SPA
  if(outcome=='D'){
    qtilde =S/sqrt(ratio) +as.vector(t(G_tilde)%*%mu)
    #The term as.vector(t(G_tilde)%*%mu) would be removed in SPAtest:::Saddle_Prob
    #keep the ratio to estimate variance of scores in Saddle_Prob
    p.single=rep(NA,ncol(G))
    for (p in 1:ncol(G)){
      p.single[p]=Saddle_Prob(q=as.vector(qtilde)[p], mu = mu, g = G_tilde[,p])$p.value
    }
  }
#   print("GLMM")
  GeneScan1D.Cauchy.window=matrix(NA,nrow=M_gene_buffer,ncol=3)
  #Burden test: for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically; for binary traits, use SPA gene- or region-based score test
  #SKAT test: for continuous traits, compute p-value use Davies; for binary traits, use SPA gene- or region-based score test

  for (m in 1:M_gene_buffer){
    # print(paste0('1D-window',m))

    #Create index for each window
    index.window<-(window.matrix_gene_buffer[,m]==1)
    G.window=G[,index.window]
    G.window=Matrix(G.window)
    #if there is no variant in this window, then do not conduct combined test in this window, move to the next one
    if(dim(G.window)[2]<=1){
      next
    }

    MAF.window<-apply(G.window,2,mean)/2
    MAC.window<-apply(G.window,2,sum)
    weight.beta_125<-dbeta(MAF.window,1,25)
    weight.beta_1<-dbeta(MAF.window,1,1)

    weight.matrix<-cbind(MAC.window<MAC.threshold,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*weight.beta_125,(MAF.window>=MAF.threshold)*weight.beta_1)
    #ultra-rare variants, rare and common variants
    colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
    weight.matrix<-Matrix(weight.matrix)

    #Single variant score test for all variants in the window, SPA p-values for binary traits
    p.single.window<-p.single[index.window]

    #approximation the covariance matrix for GLMM: t(G) P_S G
    #G.window=G.window-X%*%solve(t(X)%*%(v*X))%*%(t(X)%*%(v*G.window))
    invSigma_G.window<-solve(sparseSigma, G.window, sparse=T)

    A<-t(G.window)%*%invSigma_G.window
    B<-t(X)%*%invSigma_G.window
    K_S=A-t(B)%*%C%*%B
    #adjusted covariance matrix
    K=K_S*ratio

    #SPA gene-based tests
    if(outcome=='D'){
      V=diag(K)
      #adjusted variance
      v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
      #adjusted covariance matrix
      K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
    }

    #Burden test: for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically
    #for binary traits, calculate the SPA gene-based p-value of Burden
    p.burden<-matrix(NA,1,ncol(weight.matrix))
    for (j in 1:ncol(weight.matrix)){
      if (sum(weight.matrix[,j]!=0)>1){
        #only conduct Burden test for at least 1 variants
        temp.window.matrix<-weight.matrix[,j]
        G.window2<-as.matrix(G.window%*%temp.window.matrix)
        weights=as.vector(weight.matrix[,j])
        if(outcome=='D'){ #SPA-adjusted
          p.burden[,j]<-pchisq(as.numeric((t(G.window2)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) ;
        }else{
          #continuous
          p.burden[,j]<-pchisq(as.numeric((t(G.window2)%*%Y.res/result.null.model.GLMM$theta[1])^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) ;
        }
      }
    }

    score<-as.vector(S)[index.window]
    p.dispersion<-matrix(NA,1,ncol(weight.matrix))
    #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling based moment matching
    weight.matrix0=(MAC.window>=MAC.threshold)*weight.matrix
    for (j in 2:ncol(weight.matrix)){
      if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
        if(outcome=='D'){
          #binary
          p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K_tilde,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j])
        }else{
          #continuous
          p.dispersion[,j]<-Get.p.SKAT_noMA(score,K=K,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j])
        }
      }
    }

    p.individual1<-Get.cauchy.scan(p.single.window,as.matrix((MAC.window>=MAC.threshold & MAF.window<MAF.threshold))) #rare variants
    p.individual2<-Get.cauchy.scan(p.single.window,as.matrix((MAF.window>=MAF.threshold))) #common and low frequency variants
    p.individual<-cbind(p.burden,p.dispersion,p.individual1,p.individual2);
    colnames(p.individual)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')

    p.Cauchy<-as.matrix(apply(p.individual,1,Get.cauchy))
    #aggregated Cauchy association test
    test.common<-grep('MAF>=MAF.threshold',colnames(p.individual))
    p.Cauchy.common<-as.matrix(apply(p.individual[,test.common,drop=FALSE],1,Get.cauchy))
    p.Cauchy.rare<-as.matrix(apply(p.individual[,-test.common,drop=FALSE],1,Get.cauchy))
    GeneScan1D.Cauchy.window[m,]=c(p.Cauchy,p.Cauchy.common,p.Cauchy.rare)
  }

  GeneScan1D.Cauchy=c(Get.cauchy(GeneScan1D.Cauchy.window[,1]),Get.cauchy(GeneScan1D.Cauchy.window[,2]),Get.cauchy(GeneScan1D.Cauchy.window[,3]))
#   print("1d scan")
  ###Obtain p-values for R enhancers
  GeneScan3D.Cauchy.EnhancerAll=c()
  if(R!=0){
    for (r in 1:R){ #Loop for each enhancer
    #   print(paste0('Enhancer',r))
      if (r==1){
        G.Enhancer=as.matrix(G.EnhancerAll[,1:cumsum(p_Enhancer)[r]])
      }else{
        G.Enhancer=as.matrix(G.EnhancerAll[,(cumsum(p_Enhancer)[r-1]+1):cumsum(p_Enhancer)[r]])
      }

      G.Enhancer=Genotype_filter_Enhancer(G.Enhancer=G.Enhancer,impute.method='fixed')

      #individuals ids are matched with genotype
      G.window.Enhancer=Matrix(G.Enhancer[match.index,])
      MAF.window.Enhancer<-apply(G.window.Enhancer,2,mean)/2
      MAC.window.Enhancer<-apply(G.window.Enhancer,2,sum)

      weight.beta_125<-dbeta(MAF.window.Enhancer,1,25)
      weight.beta_1<-dbeta(MAF.window.Enhancer,1,1)
      weight.matrix<-cbind(MAC.window.Enhancer<MAC.threshold,(MAF.window.Enhancer<MAF.threshold&MAC.window.Enhancer>=MAC.threshold)*weight.beta_125,(MAF.window.Enhancer>=MAF.threshold)*weight.beta_1)
      colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
      weight.matrix<-Matrix(weight.matrix)

      #Single variant score test for all variants in the enhancer
      G_tilde.Enhancer=G.window.Enhancer-X%*%solve(t(X)%*%(v*X))%*%(t(X)%*%(v*G.window.Enhancer))

      S.Enhancer=as.vector(t(G_tilde.Enhancer)%*%Y.res)/result.null.model.GLMM$theta[1]

      ##GLMM
      #adjusted score statistics, without SPA
      if(outcome=='C'){
        invSigma_G_tilde.Enhancer<-solve(sparseSigma, G_tilde.Enhancer, sparse=T)
        V.Enhancer=t(G_tilde.Enhancer)%*%invSigma_G_tilde.Enhancer
        p.single.Enhancer=pchisq(S.Enhancer^2/(ratio*diag(V.Enhancer)),df=1,lower.tail=F)
      }
      #with SPA
      if(outcome=='D'){
        #Observed test statistic
        qtilde.Enhancer =as.vector(S.Enhancer)/sqrt(ratio) +as.vector(t(G_tilde.Enhancer)%*%mu)
        #The term as.vector(t(G_tilde.Enhancer)%*%mu) would be removed in SPAtest:::Saddle_Prob
        #keep the ratio to estimate variance of scores in Saddle_Prob
        p.single.Enhancer=rep(NA,ncol(G_tilde.Enhancer))
        for (p in 1:ncol(G_tilde.Enhancer)){
          p.single.Enhancer[p]=Saddle_Prob(q=as.vector(qtilde.Enhancer)[p], mu = mu, g = G_tilde.Enhancer[,p])$p.value
        }
      }

      p.burden.Enhancer<-matrix(NA,1,ncol(weight.matrix))
      p.dispersion.Enhancer<-matrix(NA,1,ncol(weight.matrix))

      if(length(p.single.Enhancer)>1){
        #enhancer have more than 1 variant, then conduct SKAT and burden; otherwise only conduct single variant score test
        #approximation the covariance matrix for GLMM
        #t(G) P_S G
        invSigma_G.window.Enhancer<-solve(sparseSigma, G.window.Enhancer, sparse=T)
        A<-t(G.window.Enhancer)%*%invSigma_G.window.Enhancer
        B<-t(X)%*%invSigma_G.window.Enhancer
        K_S=A-t(B)%*%C%*%B
        #adjusted covariance matrix
        K=K_S*ratio

        #SPA gene-based tests
        if(outcome=='D'){
          V=diag(K)
          #adjusted variance
          v_tilde=as.vector(S.Enhancer^2)/qchisq(p.single.Enhancer,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
          #adjusted covariance matrix
          K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
        }
      }

      #Burden
      for (j in 1:ncol(weight.matrix)){
        if (sum(weight.matrix[,j]!=0)>1){
          #only conduct Burden test for at least 1 variants
          temp.window.matrix<-weight.matrix[,j]
          G.window.Enhancer2<-as.matrix(G.window.Enhancer%*%temp.window.matrix)
          weights=as.vector(weight.matrix[,j])
          if(outcome=='D'){ #SPA-adjusted
            p.burden.Enhancer[,j]<-pchisq(as.numeric((t(G.window.Enhancer2)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F)
          }else{
            #continuous
            p.burden.Enhancer[,j]<-pchisq(as.numeric((t(G.window.Enhancer2)%*%Y.res/result.null.model.GLMM$theta[1])^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F)
          }
        }
      }

      #SKAT
      #For extremely rare variants, do not conduct SKAT
      for (j in 2:ncol(weight.matrix)){
        if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
          if(outcome=='D'){
            #binary
            p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(S.Enhancer,K=K_tilde,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j])
          }else{
            #continuous
            p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(S.Enhancer,K=K,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j])
          }
        }
      }

      p.individual1.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAC.window.Enhancer>=MAC.threshold & MAF.window.Enhancer<MAF.threshold))) #rare variants
      p.individual2.Enhancer<-Get.cauchy.scan(p.single.Enhancer,as.matrix((MAF.window.Enhancer>=MAF.threshold))) #common and low frequency variants
      p.individual.Enhancer<-cbind(p.burden.Enhancer ,p.dispersion.Enhancer,p.individual1.Enhancer,p.individual2.Enhancer);
      colnames(p.individual.Enhancer)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')

      #aggregated Cauchy association test
      p.Cauchy.Enhancer<-as.matrix(apply(p.individual.Enhancer,1,Get.cauchy))
      test.common<-grep('MAF>=MAF.threshold',colnames(p.individual.Enhancer))
      p.Cauchy.common.Enhancer<-as.matrix(apply(p.individual.Enhancer[,test.common,drop=FALSE],1,Get.cauchy))
      p.Cauchy.rare.Enhancer<-as.matrix(apply(p.individual.Enhancer[,-test.common,drop=FALSE],1,Get.cauchy))
      GeneScan3D.Cauchy.Enhancer=c(p.Cauchy.Enhancer,p.Cauchy.common.Enhancer,p.Cauchy.rare.Enhancer)

      GeneScan3D.Cauchy.EnhancerAll=rbind(GeneScan3D.Cauchy.EnhancerAll,GeneScan3D.Cauchy.Enhancer)
    }  #end of the loop of R enhancers
  }
#   print("enhancer scan")
  ##Obtain 3D windows and p-values
  #do not add promoter
  #M 1D windows + Enhancer r, r=1, ..., R
  GeneScan3D.window.EnhancerAll=c()
  if(R!=0){
    for (r in 1:dim(GeneScan3D.Cauchy.EnhancerAll)[1]){
      GeneScan3D.window.enhancer=data.frame(apply(cbind(GeneScan1D.Cauchy.window[,1],GeneScan3D.Cauchy.EnhancerAll[r,1]),1,Get.cauchy),
                                            apply(cbind(GeneScan1D.Cauchy.window[,2],GeneScan3D.Cauchy.EnhancerAll[r,2]),1,Get.cauchy),
                                            apply(cbind(GeneScan1D.Cauchy.window[,3],GeneScan3D.Cauchy.EnhancerAll[r,3]),1,Get.cauchy))
      colnames(GeneScan3D.window.enhancer)=c('all','common','rare')
      GeneScan3D.window.EnhancerAll=rbind(GeneScan3D.window.EnhancerAll,GeneScan3D.window.enhancer)
    }
  }else{
    GeneScan3D.window.enhancer=data.frame(Get.cauchy(GeneScan1D.Cauchy.window[,1]),
                                          Get.cauchy(GeneScan1D.Cauchy.window[,2]),
                                          Get.cauchy(GeneScan1D.Cauchy.window[,3]))
    colnames(GeneScan3D.window.enhancer)=c('all','common','rare')
    GeneScan3D.window.EnhancerAll=rbind(GeneScan3D.window.EnhancerAll,GeneScan3D.window.enhancer)
  }

  GeneScan3D.Cauchy.RE=GeneScan3D.window.EnhancerAll
  GeneScan3D.Cauchy=c(Get.cauchy(GeneScan3D.Cauchy.RE[,1]), Get.cauchy(GeneScan3D.Cauchy.RE[,2]), Get.cauchy(GeneScan3D.Cauchy.RE[,3]))

  ###min-p and RE with min-p
  RE_minp.all=NA;RE_minp.common=NA;RE_minp.rare=NA
  if(R!=0){
    RE.indicator=c(rep(1:R,each=M_gene_buffer))

    if(!is.infinite(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))){
      RE_minp.all=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,1]==min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))])
    }

    if(!is.infinite(min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))){
      RE_minp.common=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,2]==min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))])
    }

    if(!is.infinite(min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))){
      RE_minp.rare=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,3]==min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))])
    }
  }

  #best enhancer
  return(list(GeneScan3D.Cauchy.pvalue=GeneScan3D.Cauchy,M=M_gene_buffer,R=R,
              minp=c(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE)),
              RE_minp=cbind(RE_minp.all,RE_minp.common,RE_minp.rare)))  #GeneScan1D.Cauchy.pvalue=GeneScan1D.Cauchy,
}


######### Other functions #########
Get.p.SKAT_noMA<-function(score,K,window.matrix,weight){

  Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2) #SKAT statistics
  K.temp<-weight*t(weight*K)

  temp<-K.temp[window.matrix[,1]!=0,window.matrix[,1]!=0]
  if(sum(temp^2)==0){p<-NA}else{
    lambda=eigen(temp,symmetric=T,only.values=T)$values #eigenvalues, mixture of chi-square
    temp.p<-SKAT_davies(Q,lambda,acc=10^(-6))$Qq

    if(length(temp.p)==0 || temp.p > 1 || temp.p <= 0){
      temp.p<-Get_Liu_PVal.MOD.Lambda(Q,lambda)
    }
    p<-temp.p
  }
  return(p)
}
SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  r <- length(lambda)
  if (length(h) != r) warning("lambda and h should have the same length!")
  if (length(delta) != r) warning("lambda and delta should have the same length!")
  #out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0),PACKAGE="SKAT")
  out=davies(q, lambda, h = rep(1, length(lambda)), delta = rep(0,length(lambda)), sigma = 0, lim = 10000, acc = 0.0001)
  out$res <- 1 - out$res
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
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

Get.cauchy.scan<-function(p,window.matrix){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16)
  temp<-rep(0,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[!is.small]<-as.numeric(tan((0.5-p[!is.small])*pi))

  cct.stat<-as.numeric(t(temp)%*%window.matrix/apply(window.matrix,2,sum))
  is.large<-cct.stat>1e+15 & !is.na(cct.stat)
  is.regular<-cct.stat<=1e+15 & !is.na(cct.stat)
  pval<-rep(NA,length(cct.stat))
  pval[is.large]<-(1/cct.stat[is.large])/pi
  pval[is.regular]<-1-pcauchy(cct.stat[is.regular])
  return(pval)
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
Impute<-function(Z, impute.method){
  p<-dim(Z)[2]
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  return(as.matrix(Z))
}

Genotype_filter=function(G,pos,impute.method='fixed'){

  if(ncol(G)==0|ncol(G)==1){
    stop('Number of variants in the gene buffer region is 0 or 1')
  }

  #missing genotype imputation
  G[G==-9 | G==9]=NA
  N_MISS=sum(is.na(G))
  MISS.freq=apply(is.na(G),2,mean)

  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G=Impute(G,impute.method)
  }

  #MAF filtering
  MAF<-apply(G,2,mean)/2 #MAF of nonfiltered variants
  G[,MAF>0.5 & !is.na(MAF)]<-2-G[,MAF>0.5 & !is.na(MAF)]
  MAF<-apply(G,2,mean)/2
  s<-apply(G,2,sd)
  SNP.index<-which(MAF>0 & s!=0 & !is.na(MAF))

  check.index<-which(MAF>0 & s!=0 & !is.na(MAF)  & MISS.freq<0.1)
  if(length(check.index)<=1 ){
    stop('Number of variants with missing rate <=10% in the gene plus buffer region is <=1')
  }

  G<-Matrix(G[,SNP.index])
  pos=pos[SNP.index]
  genotype_filter=list(G=G,pos=pos)
  return(genotype_filter)
}
Genotype_filter_Enhancer=function(G.Enhancer,impute.method='fixed'){

  #missing genotype imputation
  G.Enhancer[G.Enhancer==-9 | G.Enhancer==9]=NA
  N_MISS.Enhancer=sum(is.na(G.Enhancer))
  MISS.freq.Enhancer=apply(is.na(G.Enhancer),2,mean)
  if(N_MISS.Enhancer>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS.Enhancer/nrow(G.Enhancer)/ncol(G.Enhancer))
    warning(msg,call.=F)
    G.Enhancer=Impute(G.Enhancer,impute.method)
  }

  #MAF filtering
  MAF.Enhancer<-apply(G.Enhancer,2,mean)/2 #MAF of nonfiltered variants
  G.Enhancer[,MAF.Enhancer>0.5 & !is.na(MAF.Enhancer)]<-2-G.Enhancer[,MAF.Enhancer>0.5 & !is.na(MAF.Enhancer)]
  MAF.Enhancer<-apply(G.Enhancer,2,mean)/2
  s.Enhancer<-apply(G.Enhancer,2,sd)
  SNP.index.Enhancer<-which(MAF.Enhancer>0 & s.Enhancer!=0 & !is.na(MAF.Enhancer))

  G.Enhancer<-Matrix(G.Enhancer[,SNP.index.Enhancer])
  return(G.Enhancer)
}

#####knockoff AL functions
#percentage notation
percent <- function(x, digits = 3, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}
sparse.cov.cross <- function(x,y){
  n <- nrow(x)
  cMeans.x <- colMeans(x);cMeans.y <- colMeans(y)
  covmat <- (as.matrix(crossprod(x,y)) - n*tcrossprod(cMeans.x,cMeans.y))/(n-1)
  list(cov=covmat)
}
max_nth<-function(x,n){return(sort(x,partial=length(x)-(n-1))[length(x)-(n-1)])}

