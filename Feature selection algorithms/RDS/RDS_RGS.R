# based on random sampling, data with sample column
RDS_RGS <- function(trainset,ntopgene) {
  
  # generate random gene sets
  numgene <- ncol(trainset)-2
  samplesize <- floor(nrow(trainset)/3)
  numsampleset <- 10
  numgeneset<-200
  
  datags <- t(sapply(1:numgeneset, function(x) sample(numgene,50)))
  datass <- t(sapply(1:numsampleset, function(x) sample(nrow(trainset),samplesize)))
  datapvalue<-matrix(0,numgeneset,numsampleset) 
  X <- t(trainset[,-c(1,ncol(trainset))])
  datacli <- data.frame(trainset$days,trainset[,1] %in% deadsamples[,1])
  names(datacli) <- c('time','recur')
  
  # ptm <- proc.time()
  for(kk in 1:numsampleset){
    print(kk)
    for(i in 1:numgeneset){
      data<-X[datags[i,],c(datass[kk,])]
      fa<-fanny(t(data),2,memb.exp=1.1)
      clu<-datacli[datass[kk,],]
      clu$group<-fa$clustering
      surv<-survdiff(Surv(time,recur)~group,data=clu,rho=0)
      datapvalue[i,kk]<-round(1-pchisq(surv$chisq,df=1),7)
    }
  }
  # proc.time() - ptm
  # registerDoMC(cores=10)
  # # ptm <- proc.time()
  # results <- foreach(kk = 1:numsampleset) %dopar% {
  #   for(i in 1:numgeneset){
  #     data<-X[datags[i,],c(datass[kk,])]
  #     fa<-fanny(t(data),2,memb.exp=1.1)
  #     clu<-datacli[datass[kk,],]
  #     clu$group<-fa$clustering     
  #     surv<-survdiff(Surv(time,recur)~group,data=clu,rho=0)
  #     datapvalue[i,kk]<-round(1-pchisq(surv$chisq,df=1),7)
  #   }
  #   datapvalue[,kk]
  # }
  # proc.time() - ptm
  # datapvalue <- matrix(unlist(results),ncol=20)
  
  cuttoffpvalue<-0.1 	#set cuttoff p-value of gene sets (from 0 to 1)
  cuttoffsample<-0.25 	#set ratio of passed sample set/total sample sets  (from 0 to 1)
  
  passgeneset<-c()   	
  
  #load entrez gene id in GO term gene pool
  geneid<-1:numgene
  
  #analysis gene set 
  for(i in 1:numgeneset){
    if((sum(datapvalue[i,]<cuttoffpvalue)/numsampleset)>cuttoffsample){
      tmp<-passgeneset
      passgeneset<-c(tmp,datags[i,])           
    }
  }
  
  geneidfrequency<-matrix(0,numgene,2)
  geneidfrequency[,1]<-geneid
  for(j in 1:numgene){
    geneidfrequency[j,2]<-sum(passgeneset==j)
  }
  
  #rank genes according to frequency and choose the top25 genes as signature
  geneidrank<-as.matrix(geneidfrequency[order(geneidfrequency[,2],decreasing=T),])
  topgenelist<-geneidrank[1:ntopgene,1]
  
  #save top30 genes
  # filetop30gene<-paste(path,"top30genes.txt",sep="")
  # write.table(top30geneid,filetop30gene,sep="\t",row.names=F,col.names=F,quote=F)
  topgenelist
  
  
}
