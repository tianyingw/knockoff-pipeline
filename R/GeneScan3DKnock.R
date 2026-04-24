utils::globalVariables(c('G_Enhancer1_surround','G_Enhancer2_surround',
                         'variants_Enhancer1_surround','variants_Enhancer2_surround',
                         'Enhancer1.pos','Enhancer2.pos',
                         'create.MK.AL_gene_buffer','create.MK.AL_Enhancer',
                         'KnockoffGeneration.example','GeneScan3DKnock','GeneScan3DKnock.example',
                         'G_EnhancerAll','Z_EnhancerAll','p_EnhancerAll',
                         "G_gene_buffer", "Z_gene_buffer", 'pos_gene_buffer',
                         'n','G_promoter','Z_promoter',
                         'G_Enhancer1','Z_Enhancer1','G_Enhancer2','Z_Enhancer2',
                         'G_Enhancer_surround','G_gene_buffer_surround','qchisq'))

GeneScan1D<-function(G=G_gene_buffer,Z=NULL,window.size=c(1000,5000,10000), pos=pos_gene_buffer,
                     MAC.threshold=10,MAF.threshold=0.01,Gsub.id=NULL,resampling=FALSE,result.null.model=result.null.model){
   
   #load preliminary features
   mu<-result.null.model$nullglm$fitted.values;
   Y.res<-result.null.model$Y-mu
   re.Y.res<-result.null.model$re.Y.res 
   X0<-result.null.model$X0
   outcome<-result.null.model$out_type
   
   impute.method='fixed'
   #match phenotype id and genotype id
   if(length(Gsub.id)==0){match.index<-match(result.null.model$id,1:nrow(G))}else{
      match.index<-match(result.null.model$id,Gsub.id)
   }
   if(mean(is.na(match.index))>0){
      msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
      warning(msg,call.=F)
   }
   #individuals ids are matched with genotype
   G=Matrix(G[match.index,])
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
   if(!is.null(Z)){Z<-Matrix(Z[SNP.index,])}
   
   pos=pos[SNP.index]
   
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
   
   ##single variant score tests, using fastSPA in ScoreTest_SPA function for binary traits 
   p.single<-Get.p(G,result.null.model) 
   length(p.single) 
   #score statistics
   S=t(G)%*%Y.res
   
   GeneScan1D.Cauchy.window=matrix(NA,nrow=M_gene_buffer,ncol=3)
   for (m in 1:M_gene_buffer){
      #print(paste0('1D-window',m))
      
      #Create index for each window
      index.window<-(window.matrix_gene_buffer[,m]==1)
      G.window=G[,index.window]
      G.window=Matrix(G.window)
      if(!is.null(Z)){
         Z.window=Z[index.window,]
         Z.window=Matrix(Z.window)
      }else{
         Z.window=NULL
      }
      #if there is only 1 variant in this window, then do not conduct combined test in this window, move to the next one
      if(dim(G.window)[2]==1){
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
      
      #rare variants
      if (!is.null(Z.window)){
         colnames(Z.window)<-paste0('MAF<MAF.threshold&MAC>=MAC.threshold&',1:ncol(Z.window))
         weight.matrix<-cbind(weight.matrix,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*Z.window)
         weight.matrix<-Matrix(weight.matrix)
      } 
      
      #Single variant score test for all variants in the window, SPA p-values for binary traits
      p.single.window<-p.single[index.window]
      
      if(outcome=='D'){v<-result.null.model$v}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window))}
      A<-t(G.window)%*%(v*G.window)
      B<-t(G.window)%*%(v*X0)
      C<-solve(t(X0)%*%(v*X0))
      K<-A-B%*%C%*%t(B) #covariance matrix
      
      #apply SPA gene-based tests for binary trait, deal with imbalance case-control
      if(outcome=='D'){ 
         V=diag(K)
         #adjusted variance
         v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
         #adjusted covariance matrix
         K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
      }
      
      ##Burden test
      #for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically
      #for binary traits, calculate the SPA gene-based p-value of Burden
      p.burden<-matrix(NA,1,ncol(weight.matrix))
      if(resampling==TRUE){
         for (j in 1:ncol(weight.matrix)){
            temp.window.matrix<-weight.matrix[,j]
            X<-as.matrix(G.window%*%temp.window.matrix)
            p.burden[,j]<-Get.p.base(X,result.null.model)
         }
      }else{ #do not conduct resampling-based moment matching for large sample size
         for (j in 1:ncol(weight.matrix)){
            if (sum(weight.matrix[,j]!=0)>1){ 
               #only conduct Burden test for at least 1 variants
               temp.window.matrix<-weight.matrix[,j]
               X<-as.matrix(G.window%*%temp.window.matrix)
               weights=as.vector(weight.matrix[,j])
               if(outcome=='D'){ #SPA-adjusted
                  p.burden[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) 
               }else{ 
                  #continuous
                  p.burden[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) 
               }
            }
         } 
      }
      #SKAT test
      p.dispersion<-matrix(NA,1,ncol(weight.matrix))
      score<-as.vector(S)[index.window]
      if(resampling==TRUE){
         re.score<-t(t(G.window)%*%re.Y.res) #resampling for 1000 times
         for (j in 1:ncol(weight.matrix)){
            #For extremely rare variants, do not conduct SKAT
            p.dispersion[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j]) 
         }  
      }else{
         #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling-based moment matching
         weight.matrix0=(MAC.window>=MAC.threshold)*weight.matrix
         for (j in 1:ncol(weight.matrix)){
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
   return(list(GeneScan1D.Cauchy.pvalue=GeneScan1D.Cauchy,M=M_gene_buffer))
}

GeneScan3D<-function(G=G_gene_buffer,Z=Z_gene_buffer,G.promoter=G_promoter,Z.promoter=Z_promoter,G.EnhancerAll=G_EnhancerAll,Z.EnhancerAll=Z_EnhancerAll, R=length(p_EnhancerAll),
                     p_Enhancer=p_EnhancerAll,window.size=c(1000,5000,10000),pos=pos_gene_buffer,
                     MAC.threshold=10,MAF.threshold=0.01,Gsub.id=NULL,resampling=FALSE,result.null.model=result.null.model){
   #load preliminary features
   mu<-result.null.model$nullglm$fitted.values;
   Y.res<-result.null.model$Y-mu
   re.Y.res<-result.null.model$re.Y.res 
   X0<-result.null.model$X0
   outcome<-result.null.model$out_type
   # print(length(Y.res))
   
   impute.method='fixed'
   #match phenotype id and genotype id
   if(length(Gsub.id)==0){match.index<-match(result.null.model$id,1:nrow(G))}else{
      match.index<-match(result.null.model$id,Gsub.id)
   }
   if(mean(is.na(match.index))>0){
      msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
      warning(msg,call.=F)
   }
   #individuals ids are matched with genotype
   G=Matrix(G[match.index,])
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
   if(!is.null(Z)){Z<-Matrix(Z[SNP.index,])}
   pos=pos[SNP.index]

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
   
   ##single variant score tests, using fastSPA in ScoreTest_SPA function for binary traits 
   p.single<-Get.p(G,result.null.model) 
   #score statistics
   S=t(G)%*%Y.res
   # print("single variant tests done")
   
   GeneScan1D.Cauchy.window=matrix(NA,nrow=M_gene_buffer,ncol=3)
   #Burden test: for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically; for binary traits, use SPA gene- or region-based score test
   #SKAT test: for continuous traits, compute p-value use Davies and if Davies fail to converge, we use resampling moment-based adjustment (MA); for binary traits, use SPA gene- or region-based score test
   for (m in 1:M_gene_buffer){
      #print(paste0('1D-window',m))
      
      #Create index for each window
      index.window<-(window.matrix_gene_buffer[,m]==1)
      G.window=G[,index.window]
      G.window=Matrix(G.window)
      if(!is.null(Z)){
         Z.window=Z[index.window,]
         Z.window=Matrix(Z.window)
      }else{
         Z.window=NULL
      }
      #if there is only 1 variant in this window, then do not conduct combined test in this window, move to the next one
      if(dim(G.window)[2]==1){
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
      
      #rare variants
      if (!is.null(Z.window)){
         colnames(Z.window)<-paste0('MAF<MAF.threshold&MAC>=MAC.threshold&',1:ncol(Z.window))
         weight.matrix<-cbind(weight.matrix,(MAF.window<MAF.threshold&MAC.window>=MAC.threshold)*Z.window)
         weight.matrix<-Matrix(weight.matrix)
      } 
      
      #Single variant score test for all variants in the window, SPA p-values for binary traits
      p.single.window<-p.single[index.window]
      
      if(outcome=='D'){v<-result.null.model$v}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window))}
      A<-t(G.window)%*%(v*G.window)
      B<-t(G.window)%*%(v*X0)
      C<-solve(t(X0)%*%(v*X0))
      K<-A-B%*%C%*%t(B) #covariance matrix
      
      #apply SPA gene-based tests for binary trait, deal with imbalance case-control
      if(outcome=='D'){ 
         V=diag(K)
         #adjusted variance
         v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
         #adjusted covariance matrix
         K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
      }
      
      ##Burden test
      #for continuous traits, compute p-value of Q_Burden/Scale from chi-square 1 analytically
      #for binary traits, calculate the SPA gene-based p-value of Burden
      p.burden<-matrix(NA,1,ncol(weight.matrix))
      if(resampling==TRUE){
         for (j in 1:ncol(weight.matrix)){
            temp.window.matrix<-weight.matrix[,j]
            X<-as.matrix(G.window%*%temp.window.matrix)
            p.burden[,j]<-Get.p.base(X,result.null.model)
         }
      }else{ #do not conduct resampling-based moment matching for large sample size
         for (j in 1:ncol(weight.matrix)){
            if (sum(weight.matrix[,j]!=0)>1){ 
               #only conduct Burden test for at least 1 variants
               temp.window.matrix<-weight.matrix[,j]
               X<-as.matrix(G.window%*%temp.window.matrix)
               weights=as.vector(weight.matrix[,j])
               if(outcome=='D'){ #SPA-adjusted
                  p.burden[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) 
               }else{ 
                  #continuous
                  p.burden[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) 
               }
            }
         } 
      }
      
      
      #SKAT test
      p.dispersion<-matrix(NA,1,ncol(weight.matrix))
      score<-as.vector(S)[index.window]
      if(resampling==TRUE){
         re.score<-t(t(G.window)%*%re.Y.res) #resampling for 1000 times
         for (j in 1:ncol(weight.matrix)){
            #For extremely rare variants, do not conduct SKAT
            p.dispersion[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,sum(index.window))),weight=(MAC.window>=MAC.threshold)*weight.matrix[,j]) 
         }  
      }else{
         #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling-based moment matching
         weight.matrix0=(MAC.window>=MAC.threshold)*weight.matrix
         for (j in 1:ncol(weight.matrix)){
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
   # print("1D scan done")

   ###promoter
   if(is.null(G.promoter)){
      warning('no promoter')
      GeneScan3D.Cauchy.promoter=c()
   }else{
      
      #match phenotype id and genotype id
      G.promoter=Matrix(G.promoter[match.index,])
      
      #missing genotype imputation
      G.promoter[G.promoter==-9 | G.promoter==9]=NA
      N_MISS.promoter=sum(is.na(G.promoter))
      MISS.freq.promoter=apply(is.na(G.promoter),2,mean)
      if(N_MISS.promoter>0){
         msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS.promoter/nrow(G.promoter)/ncol(G.promoter))
         warning(msg,call.=F)
         G.promoter=Impute(G.promoter,impute.method)
      }
      
      #MAF filtering
      MAF.promoter<-apply(G.promoter,2,mean)/2 #MAF of nonfiltered variants
      G.promoter[,MAF.promoter>0.5 & !is.na(MAF.promoter)]<-2-G.promoter[,MAF.promoter>0.5 & !is.na(MAF.promoter)]
      MAF.promoter<-apply(G.promoter,2,mean)/2
      s.promoter<-apply(G.promoter,2,sd)
      SNP.index.promoter<-which(MAF.promoter>0 & s.promoter!=0 & !is.na(MAF.promoter)) 
      
      G.promoter<-Matrix(G.promoter[,SNP.index.promoter])
      if(!is.null(Z.promoter)){Z.promoter<-Matrix(Z.promoter[SNP.index.promoter,])}
      
      p_promoter=dim(G.promoter)[2] 
      
      #Obtain p-value for promoter
      if (p_promoter==0){
         warning('0 variant in promoter')
         GeneScan3D.Cauchy.promoter=c()
      }else{
         
         G.window.promoter=Matrix(G.promoter)
         if(!is.null(Z.promoter)){
            Z.window.promoter=Matrix(Z.promoter)
         }else{
            Z.window.promoter=NULL
         }
         
         MAF.window.promoter<-apply(G.window.promoter,2,mean)/2
         MAC.window.promoter<-apply(G.window.promoter,2,sum)
         
         weight.beta_125<-dbeta(MAF.window.promoter,1,25)
         weight.beta_1<-dbeta(MAF.window.promoter,1,1)
         weight.matrix<-cbind(MAC.window.promoter<MAC.threshold,(MAF.window.promoter<MAF.threshold&MAC.window.promoter>=MAC.threshold)*weight.beta_125,(MAF.window.promoter>=MAF.threshold)*weight.beta_1)
         colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
         
         ##adding additional functional scores
         if (!is.null(Z.window.promoter)){
            colnames(Z.window.promoter)<-paste0('MAF<MAF.threshold&MAC>=MAC.threshold&FS',1:ncol(Z.window.promoter))
            weight.matrix<-cbind(weight.matrix,(MAF.window.promoter<MAF.threshold&MAC.window.promoter>=MAC.threshold)*Z.window.promoter)
         }
         weight.matrix<-Matrix(weight.matrix)
         
         if(outcome=='D'){v<-result.null.model$v}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window.promoter))}
         A<-t(G.window.promoter)%*%(v*G.window.promoter)
         B<-t(G.window.promoter)%*%(v*X0)
         C<-solve(t(X0)%*%(v*X0))
         K<-A-B%*%C%*%t(B) #covariance matrix
         
         #SPA gene-based tests
         if(outcome=='D'){ 
            V=diag(K)
            #adjusted variance
            v_tilde=as.vector(S^2)[index.window]/qchisq(p.single.window,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
            #adjusted covariance matrix
            K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
         }
         
         #Burden test
         p.burden.promoter<-matrix(NA,1,ncol(weight.matrix))
         if(resampling==TRUE){
            for (j in 1:ncol(weight.matrix)){
               temp.window.matrix<-weight.matrix[,j]
               X<-as.matrix(G.window.promoter%*%temp.window.matrix)
               p.burden.promoter[,j]<-Get.p.base(X,result.null.model)
            }
         }else{ #do not conduct resampling-based moment matching for large sample size
            for (j in 1:ncol(weight.matrix)){
               if (sum(weight.matrix[,j]!=0)>1){ 
                  #only conduct Burden test for at least 1 variants
                  temp.window.matrix<-weight.matrix[,j]
                  X<-as.matrix(G.window.promoter%*%temp.window.matrix)
                  weights=as.vector(weight.matrix[,j])
                  if(outcome=='D'){ #SPA-adjusted
                     p.burden.promoter[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) 
                  }else{ 
                     #continuous
                     p.burden.promoter[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) 
                  }
               }
            } 
         }
         
         #SKAT test
         score<-as.vector(t(G.window.promoter)%*%Y.res)
         p.dispersion.promoter<-matrix(NA,1,ncol(weight.matrix))
         if(resampling==TRUE){
            re.score<-t(t(G.window.promoter)%*%re.Y.res) #resampling for 1000 times
            for (j in 1:ncol(weight.matrix)){
               #For extremely rare variants, do not conduct SKAT
               p.dispersion.promoter[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,dim(G.window.promoter)[2])),weight=(MAC.window.promoter>=MAC.threshold)*weight.matrix[,j]) 
            }  
         }else{
            #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling-based moment matching
            weight.matrix0=(MAC.window.promoter>=MAC.threshold)*weight.matrix
            for (j in 1:ncol(weight.matrix)){
               if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
                  if(outcome=='D'){ 
                     #binary
                     p.dispersion.promoter[,j]<-Get.p.SKAT_noMA(score,K=K_tilde,window.matrix=as.matrix(rep(1,dim(G.window.promoter)[2])),weight=(MAC.window.promoter>=MAC.threshold)*weight.matrix[,j])
                  }else{ 
                     #continuous
                     p.dispersion.promoter[,j]<-Get.p.SKAT_noMA(score,K=K,window.matrix=as.matrix(rep(1,dim(G.window.promoter)[2])),weight=(MAC.window.promoter>=MAC.threshold)*weight.matrix[,j]) 
                  }
               }
            }
         }
         
         #Single variant score test for all variants in the window
         p.single.promoter<-Get.p(G.window.promoter,result.null.model)
         p.individual1.promoter<-Get.cauchy.scan(p.single.promoter,as.matrix((MAC.window.promoter>=MAC.threshold & MAF.window.promoter<MAF.threshold))) #rare variants
         p.individual2.promoter<-Get.cauchy.scan(p.single.promoter,as.matrix((MAF.window.promoter>=MAF.threshold))) #common and low frequency variants
         p.individual.promoter<-cbind(p.burden.promoter ,p.dispersion.promoter,p.individual1.promoter,p.individual2.promoter);
         colnames(p.individual.promoter)<-c(paste0('burden_',colnames(weight.matrix)),paste0('dispersion_',colnames(weight.matrix)),
                                            'singleCauchy_MAF<MAF.threshold&MAC>=MAC.threshold','singleCauchy_MAF>=MAF.threshold')
         
         #aggregated Cauchy association test
         p.Cauchy.promoter<-as.matrix(apply(p.individual.promoter,1,Get.cauchy))
         test.common<-grep('MAF>=MAF.threshold',colnames(p.individual.promoter))
         p.Cauchy.common.promoter<-as.matrix(apply(p.individual.promoter[,test.common,drop=FALSE],1,Get.cauchy))
         p.Cauchy.rare.promoter<-as.matrix(apply(p.individual.promoter[,-test.common,drop=FALSE],1,Get.cauchy))
         
         GeneScan3D.Cauchy.promoter=c(p.Cauchy.promoter,p.Cauchy.common.promoter,p.Cauchy.rare.promoter)
      }
   }
   # print("promoter scan done")

   ###Obtain p-values for R enhancers
   GeneScan3D.Cauchy.EnhancerAll=c()
   Enhancer_ind=0
   if(R!=0){
      Enhancer_ind=rep(TRUE,R)
      for (r in 1:R){ #Loop for each enhancer
         # print(paste0('Enhancer',r))
         if (r==1){
            G.Enhancer=G.EnhancerAll[,1:cumsum(p_Enhancer)[r],drop=FALSE]
         }else{
            G.Enhancer=G.EnhancerAll[,(cumsum(p_Enhancer)[r-1]+1):cumsum(p_Enhancer)[r],drop=FALSE]
         }
         
         if(!is.null(Z.EnhancerAll)){
            if (r==1){
               Z.Enhancer=Z.EnhancerAll[1:cumsum(p_Enhancer)[r],]
            }else{
               Z.Enhancer=Z.EnhancerAll[(cumsum(p_Enhancer)[r-1]+1):cumsum(p_Enhancer)[r],]
            }
         }else{
            Z.Enhancer=NULL
         }
         
         #individuals ids are matched with genotype
         G.Enhancer=Matrix(G.Enhancer[match.index,])
         if(!is.null(Z.Enhancer)){Z.Enhancer=Matrix(Z.Enhancer)}
         
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
         if(!is.null(Z.Enhancer)){Z.Enhancer<-Matrix(Z.Enhancer[SNP.index.Enhancer,])}
         if(dim(G.Enhancer)[2]<1){
            Enhancer_ind[r]=FALSE
            next
         }else{
            
            G.window.Enhancer=Matrix(G.Enhancer)
            if(!is.null(Z.Enhancer)){Z.window.Enhancer=Matrix(Z.Enhancer)}else{Z.window.Enhancer=NULL}
            
            MAF.window.Enhancer<-apply(G.window.Enhancer,2,mean)/2
            MAC.window.Enhancer<-apply(G.window.Enhancer,2,sum)
            
            weight.beta_125<-dbeta(MAF.window.Enhancer,1,25)
            weight.beta_1<-dbeta(MAF.window.Enhancer,1,1)
            weight.matrix<-cbind(MAC.window.Enhancer<MAC.threshold,(MAF.window.Enhancer<MAF.threshold&MAC.window.Enhancer>=MAC.threshold)*weight.beta_125,(MAF.window.Enhancer>=MAF.threshold)*weight.beta_1)
            colnames(weight.matrix)<-c('MAC<MAC.threshold','MAF<MAF.threshold&MAC>=MAC.threshold&Beta','MAF>=MAF.thresholdBeta')
            
            ##adding additional functional scores
            if (!is.null(Z.window.Enhancer)){
               colnames(Z.window.Enhancer)<-paste0('MAF<MAF.threshold&MAC>=MAC.threshold&FS',1:ncol(Z.window.Enhancer))
               weight.matrix<-cbind(weight.matrix,(MAF.window.Enhancer<MAF.threshold&MAC.window.Enhancer>=MAC.threshold)*Z.window.Enhancer)
            }
            weight.matrix<-Matrix(weight.matrix)
            
            if(outcome=='D'){v<-result.null.model$v}else{v<-rep(as.numeric(var(Y.res)),nrow(G.window.Enhancer))}
            A<-t(G.window.Enhancer)%*%(v*G.window.Enhancer)
            B<-t(G.window.Enhancer)%*%(v*X0)
            C<-solve(t(X0)%*%(v*X0))
            K<-A-B%*%C%*%t(B) #covariance matrix

            #Single variant score test for all variants in the window
            p.single.Enhancer<-Get.p(G.window.Enhancer,result.null.model)

            S.Enhancer=t(G.window.Enhancer)%*%Y.res
            #SPA gene-based tests
            if(outcome=='D'){ 
               V=diag(K)
               #adjusted variance
               # print(str(V))
               v_tilde=as.vector(S.Enhancer^2)/qchisq(p.single.Enhancer,df = 1, ncp = 0, lower.tail = FALSE,log.p = FALSE)
               v_tilde=as.vector(v_tilde)
               #adjusted covariance matrix
               # print(str(v_tilde))
               K_tilde=diag(sqrt(v_tilde/V))%*%K%*%diag(sqrt(v_tilde/V))
            }
            # print("spa")
            ##Burden test
            p.burden.Enhancer<-matrix(NA,1,ncol(weight.matrix))
            if(resampling==TRUE){
               for (j in 1:ncol(weight.matrix)){
                  temp.window.matrix<-weight.matrix[,j]
                  X<-as.matrix(G.window.Enhancer%*%temp.window.matrix)
                  p.burden.Enhancer[,j]<-Get.p.base(X,result.null.model)
               }
            }else{ #do not conduct resampling-based moment matching for large sample size
               for (j in 1:ncol(weight.matrix)){
                  if (sum(weight.matrix[,j]!=0)>1){ 
                     #only conduct Burden test for at least 1 variants
                     temp.window.matrix<-weight.matrix[,j]
                     X<-as.matrix(G.window.Enhancer%*%temp.window.matrix)
                     weights=as.vector(weight.matrix[,j])
                     if(outcome=='D'){ #SPA-adjusted
                        p.burden.Enhancer[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K_tilde%*%t(t(weights))),df=1,lower.tail=F) 
                     }else{ 
                        #continuous
                        p.burden.Enhancer[,j]<-pchisq(as.numeric((t(X)%*%Y.res)^2/weights%*%K%*%t(t(weights))),df=1,lower.tail=F) 
                     }
                  }
               } 
            }
            #SKAT test
            score<-as.vector(t(G.window.Enhancer)%*%Y.res)
            p.dispersion.Enhancer<-matrix(NA,1,ncol(weight.matrix))
            if(resampling==TRUE){
               re.score<-t(t(G.window.Enhancer)%*%re.Y.res) #resampling for 1000 times
               for (j in 1:ncol(weight.matrix)){
                  #For extremely rare variants, do not conduct SKAT
                  p.dispersion.Enhancer[,j]<-Get.p.SKAT(score,re.score,K,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j]) 
               }  
            }else{
               #For extremely rare variants, do not conduct SKAT, change MAC.threshold to 10, do not apply resampling-based moment matching
               weight.matrix0=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix
               for (j in 1:ncol(weight.matrix)){
                  if (sum(weight.matrix[,j]!=0)>1){ #only conduct SKAT test for at least 1 variants
                     if(outcome=='D'){ 
                        #binary
                        p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(score,K=K_tilde,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j])
                     }else{ 
                        #continuous
                        p.dispersion.Enhancer[,j]<-Get.p.SKAT_noMA(score,K=K,window.matrix=as.matrix(rep(1,dim(G.window.Enhancer)[2])),weight=(MAC.window.Enhancer>=MAC.threshold)*weight.matrix[,j]) 
                     }
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
         }
         GeneScan3D.Cauchy.EnhancerAll=rbind(GeneScan3D.Cauchy.EnhancerAll,GeneScan3D.Cauchy.Enhancer)
      } #end of the loop of R enhancers
   }
   # print("enhancer scan done")
   ##Obtain 3D windows and p-values
   #M 1D windows + promoter

   GeneScan3D.window.promoter=data.frame(apply(cbind(GeneScan1D.Cauchy.window[,1],GeneScan3D.Cauchy.promoter[1]),1,Get.cauchy),
                                         apply(cbind(GeneScan1D.Cauchy.window[,2],GeneScan3D.Cauchy.promoter[2]),1,Get.cauchy),
                                         apply(cbind(GeneScan1D.Cauchy.window[,3],GeneScan3D.Cauchy.promoter[3]),1,Get.cauchy))
   colnames(GeneScan3D.window.promoter)=c('all','common','rare')
   
   #M 1D windows + promoter + Enhancer r, r=1, ..., R
   GeneScan3D.window.EnhancerAll=c()
   if(R!=0){
      for (r in 1:dim(GeneScan3D.Cauchy.EnhancerAll)[1]){
         GeneScan3D.window.enhancer=data.frame(apply(cbind(GeneScan1D.Cauchy.window[,1],GeneScan3D.Cauchy.promoter[1],GeneScan3D.Cauchy.EnhancerAll[r,1]),1,Get.cauchy),
                                               apply(cbind(GeneScan1D.Cauchy.window[,2],GeneScan3D.Cauchy.promoter[2],GeneScan3D.Cauchy.EnhancerAll[r,2]),1,Get.cauchy),
                                               apply(cbind(GeneScan1D.Cauchy.window[,3],GeneScan3D.Cauchy.promoter[3],GeneScan3D.Cauchy.EnhancerAll[r,3]),1,Get.cauchy))
         colnames(GeneScan3D.window.enhancer)=c('all','common','rare')
         GeneScan3D.window.EnhancerAll=rbind(GeneScan3D.window.EnhancerAll,GeneScan3D.window.enhancer)
      }
   }
   
   GeneScan3D.Cauchy.RE=rbind(GeneScan3D.window.promoter,GeneScan3D.window.EnhancerAll)
   
   GeneScan3D.Cauchy=c(Get.cauchy(GeneScan3D.Cauchy.RE[,1]), Get.cauchy(GeneScan3D.Cauchy.RE[,2]), Get.cauchy(GeneScan3D.Cauchy.RE[,3]))
   
   ###min-p and RE with min-p
   RE.indicator=c(rep(0,M_gene_buffer),rep((1:R)[Enhancer_ind],each=M_gene_buffer))
   
   if(is.infinite(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))){
      RE_minp.all=NA
   }else{
      RE_minp.all=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,1]==min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE))])
   }
   
   if(is.infinite(min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))){
      RE_minp.common=NA
   }else{
      RE_minp.common=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,2]==min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE))])
   }
   
   if(is.infinite(min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))){
      RE_minp.rare=NA
   }else{
      RE_minp.rare=unique(RE.indicator[which(GeneScan3D.Cauchy.RE[,3]==min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE))])
   }
   
   return(list(GeneScan3D.Cauchy.pvalue=GeneScan3D.Cauchy,M=M_gene_buffer,
               minp=c(min(GeneScan3D.Cauchy.RE[,1],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,2],na.rm=TRUE),min(GeneScan3D.Cauchy.RE[,3],na.rm=TRUE)),
               RE_minp=c(RE_minp.all,RE_minp.common,RE_minp.rare)))
}

######### Other functions #########
### p-values calculation
Get.p<-function(X,result.null.model){ 
   #single variant score test: for continuous traits, score^2/v follows chi-square 1
   #for binary traits, we use fastSPA in ScoreTest_SPA function
   X<-as.matrix(X)
   mu<-result.null.model$nullglm$fitted.values;Y.res<-result.null.model$Y-mu
   outcome<-result.null.model$out_type
   if(outcome=='D'){
      p<-ScoreTest_SPA(t(X),result.null.model$Y,result.null.model$X,method=c("fastSPA"),minmac=-Inf)$p.value
   }else{
      v<-rep(as.numeric(var(Y.res)),nrow(X))
      p<-pchisq(as.vector((t(X)%*%Y.res)^2)/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.null.model$X0)%*%result.null.model$inv.X0*t(t(result.null.model$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
   }
   return(as.matrix(p))
}
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
Get.p.base<-function(X,result.null.model){
   X<-as.matrix(X)
   mu<-result.null.model$nullglm$fitted.values;Y.res<-result.null.model$Y-mu
   outcome<-result.null.model$out_type
   if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(X))}
   p<-pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.null.model$X0)%*%result.null.model$inv.X0*t(t(result.null.model$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
   p[is.na(p)]<-NA
   return(p)
}
Get.p.moment<-function(Q,re.Q){ #Q a A*q matrix of test statistics, re.Q a B*q matrix of resampled test statistics
   re.mean<-apply(re.Q,2,mean)
   re.variance<-apply(re.Q,2,var)
   re.kurtosis<-apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
   re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
   re.p<-t(pchisq((t(Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df,lower.tail=F))
   return(re.p)
}
Get.p.SKAT<-function(score,re.score,K,window.matrix,weight){
  
   Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2) #SKAT statistics
   K.temp<-weight*t(weight*K)
   
   #fast implementation by resampling based moment matching
   p0<-Get.p.moment(as.vector(t(score^2)%*%(weight*window.matrix)^2),re.score^2%*%(weight*window.matrix)^2)
   p<-p0
   for(i in which(p0<0.01 |p0>=1)){
      
      temp<-K.temp[window.matrix[,i]!=0,window.matrix[,i]!=0]
      if(sum(temp^2)==0){p[i]<-NA;next}
      
      lambda=eigen(temp,symmetric=T,only.values=T)$values
      temp.p<-SKAT_davies(Q[i],lambda,acc=10^(-6))$Qq
      
      if(length(temp.p)==0 || temp.p > 1 || temp.p <= 0){
         temp.p<-Get_Liu_PVal.MOD.Lambda(Q[i],lambda)
      }
      p[i]<-temp.p
   }
   return(as.matrix(p))
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
##Cauchy
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


# =============================================================================
# PATCH for GeneScan3DKnock.R
# Replace the single function GeneScan3D.KnockoffGeneration
#
# Changes vs original:
#  * New params: save_knockoff, load_knockoff, knockoff_file,
#                knockoff_sample_ids, stage1_only.
#  * Gene_buffer knockoff: save/load with row & column validation.
#    Input to create.MK.AL_gene_buffer is G_gene_buffer_surround
#    (post-QC surround matrix — wider than gene buffer alone).
#  * Enhancer knockoffs: always generated fresh, NOT saved.
#  * stage1_only=TRUE: save gene_buffer knockoff then return NULL.
#  * NULL null model supported when stage1_only=TRUE.
# =============================================================================

GeneScan3D.KnockoffGeneration <- function(
  G_gene_buffer_surround,
  variants_gene_buffer_surround,
  gene_buffer.pos,
  promoter.pos                  = NULL,
  R                             = 2,
  G_EnhancerAll_surround        = NULL,
  variants_EnhancerAll_surround = NULL,
  p_EnhancerAll_surround        = NULL,
  Enhancer.pos                  = NULL,
  p.EnhancerAll                 = NULL,
  Z                             = NULL,
  Z.promoter                    = NULL,
  Z.EnhancerAll                 = NULL,
  window.size                   = c(1000, 5000, 10000),
  MAC.threshold                 = 10,
  MAF.threshold                 = 0.01,
  Gsub.id                       = NULL,
  result.null.model             = NULL,
  M                             = 5,
  save_knockoff                 = FALSE,
  load_knockoff                 = FALSE,
  knockoff_file                 = NULL,
  knockoff_sample_ids           = NULL,
  stage1_only                   = FALSE
) {
  impute.method <- "fixed"

  # ---- Sample matching (supports NULL null model for stage1) ---------------
  if (is.null(result.null.model)) {
    if (!isTRUE(stage1_only))
      stop("result.null.model is NULL but stage1_only is not TRUE.")
    n           <- nrow(G_gene_buffer_surround)
    match.index <- seq_len(n)
  } else {
    mu       <- result.null.model$nullglm$fitted.values
    Y.res    <- result.null.model$Y - mu
    re.Y.res <- result.null.model$re.Y.res
    X0       <- result.null.model$X0
    outcome  <- result.null.model$out_type
    n        <- length(mu)
    if (length(Gsub.id) == 0) {
      match.index <- match(result.null.model$id,
                           seq_len(nrow(G_gene_buffer_surround)))
    } else {
      match.index <- match(result.null.model$id, Gsub.id)
    }
    if (mean(is.na(match.index)) > 0)
      warning(sprintf("Some individuals not matched with genotype. Rate = %f",
                      mean(is.na(match.index))), call. = FALSE)
  }

  # IDs in matched row order
  matched_ids <- if (!is.null(Gsub.id)) Gsub.id[match.index] else match.index
  if(mean(is.na(match.index))>0){
      msg<-sprintf("Some individuals are not matched with genotype. The rate is%f", mean(is.na(match.index)))
      warning(msg,call.=F)
   }
  # ---- QC: gene buffer surround -------------------------------------------
  G_gene_buffer_surround <- Matrix::Matrix(G_gene_buffer_surround[match.index, ])
  G_gene_buffer_surround[G_gene_buffer_surround == -9 |
                         G_gene_buffer_surround ==  9] <- NA
  N_MISS    <- sum(is.na(G_gene_buffer_surround))
  MISS.freq <- apply(is.na(G_gene_buffer_surround), 2, mean)
  if (N_MISS > 0) {
    warning(sprintf("Missing genotype rate = %f. Imputation applied.",
                    N_MISS / nrow(G_gene_buffer_surround) / ncol(G_gene_buffer_surround)),
            call. = FALSE)
    G_gene_buffer_surround <- Impute(G_gene_buffer_surround, impute.method)
  }
  MAF       <- apply(G_gene_buffer_surround, 2, mean) / 2
  G_gene_buffer_surround[, MAF > 0.5 & !is.na(MAF)] <-
    2 - G_gene_buffer_surround[, MAF > 0.5 & !is.na(MAF)]
  MAF       <- apply(G_gene_buffer_surround, 2, mean) / 2
  MAC       <- apply(G_gene_buffer_surround, 2, sum)
  s         <- apply(G_gene_buffer_surround, 2, sd)
  SNP.index <- which(MAF > 0 & s != 0 & !is.na(MAF) & MISS.freq < 0.1)
  if (length(SNP.index) <= 1) {
    warning("Number of variants passing QC in gene buffer surround is <=1", call. = FALSE)
    return(NULL)
  }
  G_gene_buffer_surround               <- Matrix::Matrix(G_gene_buffer_surround[, SNP.index])
  variants_gene_buffer_surround_filter <- variants_gene_buffer_surround[SNP.index]
  colnames(G_gene_buffer_surround)     <-
    extract_position_universal(colnames(G_gene_buffer_surround))

  # Positions within the gene buffer (subset of surround)
  positions_gene_buffer <- variants_gene_buffer_surround_filter[
    variants_gene_buffer_surround_filter <= gene_buffer.pos[2] &
    variants_gene_buffer_surround_filter >= gene_buffer.pos[1]
  ]
  if (length(positions_gene_buffer) == 0) return(NULL)

  # ---- Gene buffer knockoff: save / load / generate -----------------------
  # create.MK.AL_gene_buffer takes the POST-QC surround matrix as input.
  # It returns array [M × n_matched × p_in_gene_buffer].
  G_gene_buffer_knockoff <- .gene_ko_load_or_gen(
    load_knockoff  = load_knockoff,
    save_knockoff  = save_knockoff,
    knockoff_file  = knockoff_file,
    matched_ids    = matched_ids,
    p_expected     = length(positions_gene_buffer),
    gen_fun        = function() {
      create.MK.AL_gene_buffer(
        X                 = G_gene_buffer_surround,   # surround matrix
        pos               = variants_gene_buffer_surround_filter,
        gene_buffer_start = gene_buffer.pos[1],
        gene_buffer_end   = gene_buffer.pos[2],
        M                 = M,
        corr_max          = 0.75,
        maxN.neighbor     = Inf,
        maxBP.neighbor    = 10000,
        corr_base         = 0.05,
        n.AL              = floor(10 * n^(1/3) * log(n)),
        thres.ultrarare   = 25,
        R2.thres          = 0.75
      )
    },
    snp_pos        = positions_gene_buffer
  )
  if (is.null(G_gene_buffer_knockoff)) return(NULL)

  # Genotype matrices for test region
  G_gene_buffer <- G_gene_buffer_surround[,
    variants_gene_buffer_surround_filter %in% positions_gene_buffer
  ]
  G_promoter <- NULL
  if (!is.null(promoter.pos)) {
    positions_promoter <- positions_gene_buffer[
      positions_gene_buffer <= promoter.pos[2] &
      positions_gene_buffer >= promoter.pos[1]
    ]
    G_promoter <- G_gene_buffer[, positions_gene_buffer %in% positions_promoter]
  }

  # Functional annotations
  Z_gene_buffer <- NULL
  if (!is.null(Z)) {
    pos_gb_nf     <- variants_gene_buffer_surround[
      variants_gene_buffer_surround <= gene_buffer.pos[2] &
      variants_gene_buffer_surround >= gene_buffer.pos[1]
    ]
    Z_gene_buffer <- as.matrix(Z[pos_gb_nf %in% positions_gene_buffer, ])
  }
  Z_promoter <- NULL
  if (!is.null(Z.promoter) && !is.null(promoter.pos)) {
    pos_pr_nf  <- variants_gene_buffer_surround[
      variants_gene_buffer_surround <= promoter.pos[2] &
      variants_gene_buffer_surround >= promoter.pos[1]
    ]
    Z_promoter <- as.matrix(Z.promoter[pos_pr_nf %in% positions_promoter, ])
  }

  # ---- Stage 1: done after saving knockoff --------------------------------
  if (isTRUE(stage1_only)) return(invisible(NULL))

  # ---- R enhancers (always generated fresh — not saved) -------------------
  G_EnhancerAll          <- c()
  p_EnhancerAll_out      <- c()
  Z_EnhancerAll_out      <- c()
  G_EnhancerAll_knockoff <- c()

  if (R != 0) {
    for (r in seq_len(R)) {
      if (r == 1) {
        G_Enh_surround   <- G_EnhancerAll_surround[,
          seq_len(cumsum(p_EnhancerAll_surround)[r]), drop = FALSE]
        pos_Enh_surround <- variants_EnhancerAll_surround[
          seq_len(cumsum(p_EnhancerAll_surround)[r])]
      } else {
        G_Enh_surround   <- G_EnhancerAll_surround[,
          (cumsum(p_EnhancerAll_surround)[r - 1] + 1):
           cumsum(p_EnhancerAll_surround)[r], drop = FALSE]
        pos_Enh_surround <- variants_EnhancerAll_surround[
          (cumsum(p_EnhancerAll_surround)[r - 1] + 1):
           cumsum(p_EnhancerAll_surround)[r]]
      }

      # QC: enhancer surround
      G_Enh_surround <- Matrix::Matrix(G_Enh_surround[match.index, ])
      G_Enh_surround[G_Enh_surround == -9 | G_Enh_surround == 9] <- NA
      N_MISS    <- sum(is.na(G_Enh_surround))
      MISS.freq <- apply(is.na(G_Enh_surround), 2, mean)
      if (N_MISS > 0) {
        warning(sprintf("Enhancer %d: missing rate = %f. Imputation applied.", r,
                        N_MISS / nrow(G_Enh_surround) / ncol(G_Enh_surround)),
                call. = FALSE)
        G_Enh_surround <- Impute(G_Enh_surround, impute.method)
      }
      MAF       <- apply(G_Enh_surround, 2, mean) / 2
      G_Enh_surround[, MAF > 0.5 & !is.na(MAF)] <-
        2 - G_Enh_surround[, MAF > 0.5 & !is.na(MAF)]
      MAF       <- apply(G_Enh_surround, 2, mean) / 2
      s         <- apply(G_Enh_surround, 2, sd)
      SNP.index <- which(MAF > 0 & s != 0 & !is.na(MAF) & MISS.freq < 0.1)
      if (length(SNP.index) <= 1) {
        warning(sprintf("Enhancer %d: variants passing QC <=1; skipping.", r), call. = FALSE)
        next
      }
      G_Enh_surround <- Matrix::Matrix(G_Enh_surround[, SNP.index])
      pos_Enh_filter <- pos_Enh_surround[SNP.index]
      colnames(G_Enh_surround) <- extract_position_universal(colnames(G_Enh_surround))

      # Generate enhancer knockoff fresh (NOT saved)
      # BUG FIX: was create.MK.AL_Enhancer(..., M=5) — hardcoded
      G_Enh_knockoff <- create.MK.AL_Enhancer(
        X               = G_Enh_surround,          # surround matrix
        pos             = pos_Enh_filter,
        Enhancer_start  = as.numeric(Enhancer.pos[r, 1]),
        Enhancer_end    = as.numeric(Enhancer.pos[r, 2]),
        M               = M,                       
        corr_max        = 0.75,
        maxN.neighbor   = Inf,
        maxBP.neighbor  = 10000,
        corr_base       = 0.05,
        n.AL            = floor(10 * n^(1/3) * log(n)),
        thres.ultrarare = 25,
        R2.thres        = 0.75
      )

      positions_enhancer <- pos_Enh_filter[
        pos_Enh_filter <= Enhancer.pos[r, 2] &
        pos_Enh_filter >= Enhancer.pos[r, 1]
      ]
      G_enhancer             <- Matrix::Matrix(
        G_Enh_surround[, pos_Enh_filter %in% positions_enhancer])
      G_EnhancerAll          <- cbind(G_EnhancerAll, G_enhancer)
      p_EnhancerAll_out      <- c(p_EnhancerAll_out, length(positions_enhancer))
      G_EnhancerAll_knockoff <- abind::abind(G_EnhancerAll_knockoff, G_Enh_knockoff)

      # Functional annotation
      if (!is.null(Z.EnhancerAll)) {
        if (r == 1) {
          Z_Enh <- as.matrix(Z.EnhancerAll[seq_len(cumsum(p.EnhancerAll)[r]), ])
        } else {
          Z_Enh <- as.matrix(Z.EnhancerAll[
            (cumsum(p.EnhancerAll)[r - 1] + 1):cumsum(p.EnhancerAll)[r], ])
        }
        Z_Enh <- as.matrix(Z_Enh[
          pos_Enh_surround[
            pos_Enh_surround <= Enhancer.pos[r, 2] &
            pos_Enh_surround >= Enhancer.pos[r, 1]
          ] %in% positions_enhancer, ])
        Z_EnhancerAll_out <- rbind(Z_EnhancerAll_out, Z_Enh)
      }
    }
  }

  # ---- Association tests ---------------------------------------------------
  GeneScan3D.Cauchy <- GeneScan3D(
    G              = G_gene_buffer,
    Z              = Z_gene_buffer,
    G.promoter     = G_promoter,
    Z.promoter     = Z_promoter,
    G.EnhancerAll  = G_EnhancerAll,
    Z.EnhancerAll  = Z_EnhancerAll_out,
    R              = R,
    p_Enhancer     = p_EnhancerAll_out,
    window.size    = window.size,
    pos            = positions_gene_buffer,
    MAC.threshold  = MAC.threshold,
    MAF.threshold  = MAF.threshold,
    result.null.model = result.null.model,
    Gsub.id        = Gsub.id[match.index]
  )$GeneScan3D.Cauchy.pvalue

  GeneScan3D.Cauchy_knockoff <- matrix(NA, nrow = M, ncol = 3)
  for (k in seq_len(M)) {
    G_gbk              <- G_gene_buffer_knockoff[k, , ]
    G_prom_k           <- NULL
    if (!is.null(promoter.pos))
      G_prom_k <- G_gbk[, positions_gene_buffer %in% positions_promoter]

    GeneScan3D.Cauchy_knockoff[k, ] <- GeneScan3D(
      G              = G_gbk,
      Z              = Z_gene_buffer,
      G.promoter     = G_prom_k,
      Z.promoter     = Z_promoter,
      G.EnhancerAll  = if (R > 0 && length(G_EnhancerAll_knockoff) > 0)
                         G_EnhancerAll_knockoff[k, , ] else NULL,
      Z.EnhancerAll  = Z_EnhancerAll_out,
      R              = R,
      p_Enhancer     = p_EnhancerAll_out,
      window.size    = window.size,
      pos            = positions_gene_buffer,
      MAC.threshold  = MAC.threshold,
      MAF.threshold  = MAF.threshold,
      result.null.model = result.null.model,
      Gsub.id        = Gsub.id[match.index]
    )$GeneScan3D.Cauchy.pvalue
  }

  return(list(
    GeneScan3D.Cauchy          = GeneScan3D.Cauchy,
    GeneScan3D.Cauchy_knockoff = GeneScan3D.Cauchy_knockoff
  ))
}

GeneScan3DKnock<-function(M=5,p0=GeneScan3DKnock.example$GeneScan3D.original,
                          p_ko=cbind(GeneScan3DKnock.example$GeneScan3D.ko1,
                                     GeneScan3DKnock.example$GeneScan3D.ko2,
                                     GeneScan3DKnock.example$GeneScan3D.ko3,
                                     GeneScan3DKnock.example$GeneScan3D.ko4,
                                     GeneScan3DKnock.example$GeneScan3D.ko5),fdr = 0.1,gene_id=GeneScan3DKnock.example$gene.id){
   
   p=cbind(p0,p_ko)
   #calculate knockoff statistics W, kappa, tau for given original p-value and M knockoff p-values
   T=-log10(p)
   

   W=(T[,1]-apply(T[,2:(M+1)],1,median))*(T[,1]>=apply(T[,2:(M+1)],1,max))
   kappa=apply(T,1,which.max)-1 #max T is from original data (0) or knockoff data (1 to 5)
   tau=apply(T,1,max)-apply(T,1,function(x)median(x[-which.max(x)]))
   Rej.Bound=10000 
   b=order(tau,decreasing=T)
   c_0=kappa[b]==0  #only calculate q-value for kappa=0
   # print("success")
   #calculate ratios for top Rej.Bound tau values
   ratio<-c();temp_0<-0
   for(i in 1:length(b)){
      temp_0<-temp_0+c_0[i]
      temp_1<-i-temp_0
      temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
      ratio<-c(ratio,temp_ratio)
      if(i>Rej.Bound){break}
   }
   # print("success")
   #calculate q value for each gene/window
   qvalue=rep(1,length(tau))
   for(i in 1:length(b)){
      qvalue[b[i]]=min(ratio[i:min(length(b),Rej.Bound)])*c_0[i]+1-c_0[i] #only calculate q-value for kappa=0, q-value for kappa!=0 is 1
      if(i>Rej.Bound){break}
   }
   print(table(qvalue))
   #W statistics threshold
   W.threshold=MK.threshold.byStat(kappa,tau,M=M,fdr=fdr,Rej.Bound=Rej.Bound)
   
   #gene is significant if its q value less or equal than the fdr threshold; OR W>=W.threshold
   gene_sign=as.character(gene_id[which(qvalue<=fdr)])
   
   return(list(W=W,W.threshold=W.threshold,Qvalue=qvalue,gene_sign=gene_sign))
}


######### Other functions #########
#Knockoff generation for gene buffer regions
create.MK.AL_gene_buffer <- function(X=G_gene_buffer_surround,pos,gene_buffer_start,gene_buffer_end,M,corr_max=0.75,maxN.neighbor=Inf,
                                     maxBP.neighbor=100000,corr_base=0.05,n.AL=floor(10*n^(1/3)*log(n)),
                                     thres.ultrarare=25,R2.thres=0.75) {
   
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
   
   X.AL<-w*X[index.AL,] #n.AL samples
   
   sparse.fit<-sparse.cor(X.AL)
   cor.X.AL<-sparse.fit$cor;cov.X.AL<-sparse.fit$cov
   skip.index<-colSums(X.AL!=0)<=thres.ultrarare #skip features that are ultra sparse, permutation will be directly applied to generate knockoffs
   
   Sigma.distance = as.dist(1 - abs(cor.X))
   if(ncol(X)>1){
      fit = hclust(Sigma.distance, method="single") #hierarchical clustering
      corr_max = corr_max
      clusters = cutree(fit, h=1-corr_max)  #variants from two different clusters do not have a correlation greater than 0.75.
   }else{clusters<-1}
   
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
   
   G_gene_buffer_knockoff <- array(0, dim = c(M, nrow(X), length(snps_ind)))
   for (j in 1:M) {
      G_gene_buffer_knockoff[j, ,] <-X_k[[j]][,snps_ind]
   }
   return(G_gene_buffer_knockoff)
}

#Knockoff generation for Enhancer
create.MK.AL_Enhancer <- function(X=G_Enhancer_surround,pos,Enhancer_start,Enhancer_end,M,corr_max=0.75,maxN.neighbor=Inf,
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
   ##only focus on snps within enhancer
   for(k in 1:M){
      X_k[[k]]<-matrix(0,nrow=nrow(X),ncol=ncol(X))
      #X_k[[k]]<-big.matrix(nrow=nrow(X),ncol=ncol(X),init=0,shared=FALSE)
   }
   
   snps_ind=which(pos<=Enhancer_end&pos>=Enhancer_start)
   
   index.exist<-c()
   for (k in unique(clusters[snps_ind])){
      #print(paste0('cluster',k))
      cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
      for(i in which(clusters==k)[which(clusters==k)%in%snps_ind]){ 
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
   
   G_Enhancer_knockoff <- array(0, dim = c(M, nrow(X), length(snps_ind)))
   for (j in 1:M) {
      G_Enhancer_knockoff[j, ,] <-X_k[[j]][,snps_ind]
   }
   return(G_Enhancer_knockoff)
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
#knockoff filter
MK.threshold.byStat<-function (kappa,tau,M,fdr = 0.1,Rej.Bound=10000){
   b<-order(tau,decreasing=T)
   c_0<-kappa[b]==0
   ratio<-c();temp_0<-0
   for(i in 1:length(b)){
      #if(i==1){temp_0=c_0[i]}
      temp_0<-temp_0+c_0[i]
      temp_1<-i-temp_0
      temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
      ratio<-c(ratio,temp_ratio)
      if(i>Rej.Bound){break}
   }
   ok<-which(ratio<=fdr)
   if(length(ok)>0){
      #ok<-ok[which(ok-ok[1]:(ok[1]+length(ok)-1)<=0)]
      return(tau[b][ok[length(ok)]])
   }else{return(Inf)}
}

