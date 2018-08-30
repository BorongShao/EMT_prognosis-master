cluster_measure <- function(nothresdata,numcluster,directionadd){
  originalclr <- tsne_km_cl(nothresdata,1:(ncol(nothresdata)-2),0,algorithmstrings[i],numcluster)
  i=1
  ttestclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=2
  lassoclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=3
  netlassoclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=4
  addgfs <- addg_sf[1:numtopnets]
  # addgfs <- addg_sf[1:20]
  moldata <- nothresdata[,-c(1,ncol(nothresdata))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[x])
  newX2 <- sapply(addgfs,function(x) rowSums(newX[,x])/length(x))
  nothresdatanet <- data.frame(nothresdata[,1],newX2,nothresdata$days)
  names(nothresdatanet)[ncol(nothresdatanet)] <- "days"
  addgclr <- tsne_km_cl(nothresdatanet,1:ncol(newX2),0,algorithmstrings[i],numcluster)
  i=5
  netrankclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=6
  stsvmclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=7
  coxclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=8
  coxregclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=9
  rdsclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=10
  survnetfs <- survnet_sf[1:numtopnets]
  # survnetfs <- survnet_sf[1:20]
  newX2 <- sapply(survnetfs,function(x) rowSums(moldata[,x])/length(x))
  nothresdatanet <- data.frame(nothresdata[,1],newX2,nothresdata$days)
  names(nothresdatanet)[ncol(nothresdatanet)] <- "days"
  survnetclr <- tsne_km_cl(nothresdatanet,1:ncol(newX2),0,algorithmstrings[i],numcluster)
  i=11
  allrankclr <- tsne_km_cl(nothresdata,alg10sfs,0,"all",numcluster)
  
  numbers <- list(originalclr,ttestclr,lassoclr,netlassoclr,addgclr,netrankclr,stsvmclr,coxclr,coxregclr,rdsclr,survnetclr,allrankclr)
  features <- list(ttestfs,lassofs,netlassofs, addgfs, netrankfs, stsvmfs, coxfs, coxregfs, rdsfs, survnetfs,alg10sfs)
  list(numbers,features)
}

tsne_km_cl <- function(nothresdata, features, usetsne=1, algorithm,numcluster){
  
  days <- nothresdata$days
  isdead <- nothresdata[,1] %in% deadsamples[,1]
  survobj <- Surv(days, isdead)
  data <- nothresdata[,features+1]
  
  # tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
  # if(usetsne==1)
  #   d_tsne_1 = as.data.frame(tsne_model_1$Y) else
  #     d_tsne_1 = data
  
  d_tsne_1_original=NULL
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)
  d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit_cluster_hierarchical=hclust(dist(scale(data)))
  d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=numcluster))
  
  sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_hierarchical)$chisq
  pval_hierarchical <- 1- pchisq(sdiff,numcluster-1) 
 
  #for saving plot with correct names. 
  
  fit <- survfit(survobj ~ cl_kmeans,data = d_tsne_1_original)
  # filestring <- paste('../Thesis_Results/alllevels/clustering/', ncol(nothresdata)-2, algorithm, '_kmeans.pdf', sep='')
  # pdf(filestring, width = 6.5,height = 5)
  
  #the following five lines are only for plot.
  # data <- data.frame(Data1,Data2)
  # fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  # cl_kmeans = factor(fit_cluster_kmeans$cluster)
  # fit<- survfit(survobj ~ cl_kmeans)
  # sdiff <- survdiff(survobj ~ cl_kmeans)
  
  # ggsurvplot(
  #   fit,                     # survfit object with calculated statistics.
  #   # data = d_tsne_1_original,  # data used to fit survival curves.
  #   risk.table = TRUE,       # show risk table.
  #   pval = paste("p=",format(pval_kmeans,digits=3),sep=""),             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for
  #   # point estimaes of survival curves.
  #   xlim = c(0,2000),        # present narrower X axis, but not affect
  #   # survival estimates.
  #   break.time.by = 500,     # break X axis in time intervals by 500.
  #   ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  #   risk.table.y.text.col = T, # colour risk table text annotations.
  #   risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
    # title = "Lasso"
  # )
  # dev.off()
  # 
  # fit <- survfit(survobj ~ cl_hierarchical,data=d_tsne_1_original)
  # filestring <- paste('../Thesis_Results/alllevels/clustering/', ncol(nothresdata)-2, algorithm, '_hierarchical.pdf', sep='')
  # pdf(filestring, width = 6.5,height = 5)
  # ggsurvplot(
  #   fit,                     # survfit object with calculated statistics.
  #   data = d_tsne_1_original,  # data used to fit survival curves. 
  #   risk.table = TRUE,       # show risk table.
  #   pval = TRUE,             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for 
  #   # point estimaes of survival curves.
  #   xlim = c(0,2000),        # present narrower X axis, but not affect
  #   # survival estimates.
  #   break.time.by = 500,     # break X axis in time intervals by 500.
  #   ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  #   risk.table.y.text.col = T, # colour risk table text annotations.
  #   risk.table.y.text = FALSE # show bars instead of names in text annotations
  #   # in legend of risk table
  # )
  # dev.off()
  
  # log-rank test
  # sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_kmeans)$chisq
  # pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_hierarchical)$chisq
  # pval_hierarchical <- 1- pchisq(sdiff,numcluster-1) 
  c(pval_kmeans,as.numeric(table(d_tsne_1_original$cl_kmeans)),pval_hierarchical,as.numeric(table(d_tsne_1_original$cl_hierarchical)))
}

SNF_clustering <- function(datalist, featurelist, numclusterm ,numf2dt,numf3dt,numf2dtnet, numf3dtnet,directionadd){
  K = 10;		# number of neighbors, usually (10~30)
  alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
  T = 10; 	# Number of Iterations, usually (10~20)
  C = numcluster			
  days <- datalist[[1]]$days
  isdead <- datalist[[1]][,1] %in% deadsamples[,1]
  survobj <- Surv(days, isdead)
  pvaluesmatrix <- matrix(0,7,12)
  for(j in c(1:3,5:9,11)){
    for(i in 1:3){
      data <- datalist[[i]][,featurelist[[i]][[j]]+1]
      Dist1 = dist2(as.matrix(data),as.matrix(data));
      W = affinityMatrix(Dist1, K, alpha)
      group = spectralClustering(W, C);
      sdiff <- survdiff(survobj ~ group)$chisq
      pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    }
    i=4
    # combine GE and DM
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:numf2dt]+1]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W = SNF(list(W1,W2), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
    i=5
    # combine GE and CNA
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:numf2dt]+1]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W = SNF(list(W1,W2), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
    i=6
    # combine DM and CNA
    Data1 <- datalist[[2]][,featurelist[[2]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:numf2dt]+1]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W = SNF(list(W1,W2), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
    i=7
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf3dt]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:numf3dt]+1]
    Data3 <- datalist[[3]][,featurelist[[3]][[j]][1:numf3dt]+1]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    Dist3 = dist2(as.matrix(Data3),as.matrix(Data3));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W3 = affinityMatrix(Dist3, K, alpha)
    W = SNF(list(W1,W2,W3), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  }
  j=4
    moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
    newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[1]][x])
    dataadd1 <- sapply(featurelist[[1]][[4]],function(x) rowSums(newX[,x])/length(x))
    
    moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
    newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[2]][x])
    dataadd2 <- sapply(featurelist[[2]][[4]],function(x) rowSums(newX[,x])/length(x))
    
    moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
    newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[3]][x])
    dataadd3 <- sapply(featurelist[[3]][[4]],function(x) rowSums(newX[,x])/length(x))
    
    dataaddlist <- list(dataadd1,dataadd2,dataadd3)
    
    for(i in 1:3){
      data <- dataaddlist[[i]]
      Dist1 = dist2(as.matrix(data),as.matrix(data));
      W = affinityMatrix(Dist1, K, alpha)
      group = spectralClustering(W, C);
      sdiff <- survdiff(survobj ~ group)$chisq
      pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    }
    i=4
    # combine GE and DM
    Data1 <- dataaddlist[[1]][,1:numf2dtnet]
    Data2 <- dataaddlist[[2]][,1:numf2dtnet]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W = SNF(list(W1,W2), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
    i=5
    # combine GE and CNA
    Data1 <- dataaddlist[[1]][,1:numf2dtnet]
    Data2 <- dataaddlist[[3]][,1:numf2dtnet]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W = SNF(list(W1,W2), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
    i=6
    # combine DM and CNA
    Data1 <- dataaddlist[[2]][,1:numf2dtnet]
    Data2 <- dataaddlist[[3]][,1:numf2dtnet]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W = SNF(list(W1,W2), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
    i=7
    Data1 <- dataaddlist[[1]][,1:numf3dtnet]
    Data2 <- dataaddlist[[2]][,1:numf3dtnet]
    Data3 <- dataaddlist[[3]][,1:numf3dtnet]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    Dist3 = dist2(as.matrix(Data3),as.matrix(Data3));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W3 = affinityMatrix(Dist3, K, alpha)
    W = SNF(list(W1,W2,W3), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
  j=10
  
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  dataadd1 <- sapply(featurelist[[1]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  dataadd2 <- sapply(featurelist[[2]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  dataadd3 <- sapply(featurelist[[3]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    Dist1 = dist2(as.matrix(data),as.matrix(data));
    W = affinityMatrix(Dist1, K, alpha)
    group = spectralClustering(W, C);
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W = SNF(list(W1,W2), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W = SNF(list(W1,W2), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W = SNF(list(W1,W2), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  Dist3 = dist2(as.matrix(Data3),as.matrix(Data3));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W3 = affinityMatrix(Dist3, K, alpha)
  W = SNF(list(W1,W2,W3), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  
  j=12
  featurelist <- list(c(1:(ncol(datalist[[1]])-2)),c(1:(ncol(datalist[[2]])-2)),c(1:(ncol(datalist[[3]])-2)))
  for(i in 1:3){
    data <- datalist[[i]][,featurelist[[i]]+1]
    Dist1 = dist2(as.matrix(data),as.matrix(data));
    W = affinityMatrix(Dist1, K, alpha)
    group = spectralClustering(W, C);
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  }
  i=4
  # combine GE and DM
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W = SNF(list(W1,W2), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  i=5
  # combine GE and CNA
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W = SNF(list(W1,W2), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  i=6
  # combine DM and CNA
  Data1 <- datalist[[2]][,featurelist[[2]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W = SNF(list(W1,W2), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  i=7
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  Data3 <- datalist[[3]][,featurelist[[3]]+1]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  Dist3 = dist2(as.matrix(Data3),as.matrix(Data3));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W3 = affinityMatrix(Dist3, K, alpha)
  W = SNF(list(W1,W2,W3), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  pvaluesmatrix
}

SNF_clustering_wholedata <- function(datalist, numcluster){
  K = 10;		# number of neighbors, usually (10~30)
  alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
  T = 10; 	# Number of Iterations, usually (10~20)
  C = numcluster			
  days <- datalist[[1]]$days
  isdead <- datalist[[1]][,1] %in% deadsamples[,1]
  survobj <- Surv(days, isdead)
  pvaluesmatrix <- matrix(0,7,12)
  for(j in 1:11){
    for(i in 1:3){
      data <- datalist[[i]][,featurelist[[i]][[j]]+1]
      Dist1 = dist2(as.matrix(data),as.matrix(data));
      W = affinityMatrix(Dist1, K, alpha)
      group = spectralClustering(W, C);
      sdiff <- survdiff(survobj ~ group)$chisq
      pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    }
    i=4
    # combine GE and DM
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:10]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:10]+1]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W = SNF(list(W1,W2), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
    i=5
    # combine GE and CNA
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:10]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:10]+1]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W = SNF(list(W1,W2), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
    i=6
    # combine DM and CNA
    Data1 <- datalist[[2]][,featurelist[[2]][[j]][1:10]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:10]+1]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W = SNF(list(W1,W2), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
    
    i=7
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:7]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:7]+1]
    Data3 <- datalist[[3]][,featurelist[[3]][[j]][1:7]+1]
    Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
    Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
    Dist3 = dist2(as.matrix(Data3),as.matrix(Data3));
    W1 = affinityMatrix(Dist1, K, alpha)
    W2 = affinityMatrix(Dist2, K, alpha)
    W3 = affinityMatrix(Dist3, K, alpha)
    W = SNF(list(W1,W2,W3), K, T)
    group = spectralClustering(W, C); 	# the final subtypes information
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  }
  j=12
  featurelist <- list(c(1:(ncol(datalist[[1]])-2)),c(1:(ncol(datalist[[2]])-2)),c(1:(ncol(datalist[[3]])-2)))
  for(i in 1:3){
    data <- datalist[[i]][,featurelist[[i]]+1]
    Dist1 = dist2(as.matrix(data),as.matrix(data));
    W = affinityMatrix(Dist1, K, alpha)
    group = spectralClustering(W, C);
    sdiff <- survdiff(survobj ~ group)$chisq
    pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  }
  i=4
  # combine GE and DM
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W = SNF(list(W1,W2), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  i=5
  # combine GE and CNA
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W = SNF(list(W1,W2), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  i=6
  # combine DM and CNA
  Data1 <- datalist[[2]][,featurelist[[2]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W = SNF(list(W1,W2), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  
  i=7
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  Data3 <- datalist[[3]][,featurelist[[3]]+1]
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  Dist3 = dist2(as.matrix(Data3),as.matrix(Data3));
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  W3 = affinityMatrix(Dist3, K, alpha)
  W = SNF(list(W1,W2,W3), K, T)
  group = spectralClustering(W, C); 	# the final subtypes information
  sdiff <- survdiff(survobj ~ group)$chisq
  pvaluesmatrix[i,j] <- 1- pchisq(sdiff,numcluster-1) 
  pvaluesmatrix
}

# multilevel_kmeans_clustering(datatypes455_common3,listfeatures455,3,10,7,7,5,add_dire455)
multilevel_kmeans_clustering <- function(datalist, featurelist, numcluster, numf2dt,numf3dt,numf2dtnet, numf3dtnet,directionadd){
  days <- datalist[[1]]$days
  isdead <- datalist[[1]][,1] %in% deadsamples[,1]
  survobj <- Surv(days, isdead)
  pvaluesmatrix <- matrix(0,7,12)
  for(j in c(1:3,5:9,11)){
    for(i in 1:3){
      data <- datalist[[i]][,featurelist[[i]][[j]]+1]
      fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
      cl_kmeans = factor(fit_cluster_kmeans$cluster)
      fit <- survfit(survobj ~ cl_kmeans)
      sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
      pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
      # print(pval_kmeans)
      pvaluesmatrix[i,j]  <- pval_kmeans
    }
    i=4
    # combine GE and DM
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=5
    # combine GE and CNA
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=6
    # combine DM and CNA
    Data1 <- datalist[[2]][,featurelist[[2]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=7
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf3dt]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:numf3dt]+1]
    Data3 <- datalist[[3]][,featurelist[[3]][[j]][1:numf3dt]+1]
    data <- data.frame(Data1,Data2,Data3)
    fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  
  j=4
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[1]][x])
  dataadd1 <- sapply(featurelist[[1]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[2]][x])
  dataadd2 <- sapply(featurelist[[2]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[3]][x])
  dataadd3 <- sapply(featurelist[[3]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  j=10
  
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  dataadd1 <- sapply(featurelist[[1]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  dataadd2 <- sapply(featurelist[[2]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  dataadd3 <- sapply(featurelist[[3]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  j=12
  featurelist <- list(c(1:(ncol(datalist[[1]])-2)),c(1:(ncol(datalist[[2]])-2)),c(1:(ncol(datalist[[3]])-2)))
  for(i in 1:3){
    data <- datalist[[i]][,featurelist[[i]]+1]
    fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- datalist[[2]][,featurelist[[2]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  Data3 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=kmeans(scale(data), numcluster,iter.max = 100)  
  cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  pvaluesmatrix
}

multilevel_hierarchical_clustering <- function(datalist, featurelist, numcluster, numf2dt,numf3dt,numf2dtnet, numf3dtnet,directionadd){
  days <- datalist[[1]]$days
  isdead <- datalist[[1]][,1] %in% deadsamples[,1]
  survobj <- Surv(days, isdead)
  pvaluesmatrix <- matrix(0,7,12)
  for(j in c(1:3,5:9,11)){
    for(i in 1:3){
      data <- datalist[[i]][,featurelist[[i]][[j]]+1]
      fit_cluster_kmeans=hclust(dist(scale(data)))
      cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
      fit <- survfit(survobj ~ cl_kmeans)
      sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
      pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
      # print(pval_kmeans)
      pvaluesmatrix[i,j]  <- pval_kmeans
    }
    i=4
    # combine GE and DM
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=5
    # combine GE and CNA
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=6
    # combine DM and CNA
    Data1 <- datalist[[2]][,featurelist[[2]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=7
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf3dt]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:numf3dt]+1]
    Data3 <- datalist[[3]][,featurelist[[3]][[j]][1:numf3dt]+1]
    data <- data.frame(Data1,Data2,Data3)
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  
  j=4
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[1]][x])
  dataadd1 <- sapply(featurelist[[1]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[2]][x])
  dataadd2 <- sapply(featurelist[[2]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[3]][x])
  dataadd3 <- sapply(featurelist[[3]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  j=10
  
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  dataadd1 <- sapply(featurelist[[1]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  dataadd2 <- sapply(featurelist[[2]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  dataadd3 <- sapply(featurelist[[3]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  j=12
  featurelist <- list(c(1:(ncol(datalist[[1]])-2)),c(1:(ncol(datalist[[2]])-2)),c(1:(ncol(datalist[[3]])-2)))
  for(i in 1:3){
    data <- datalist[[i]][,featurelist[[i]]+1]
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- datalist[[2]][,featurelist[[2]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  Data3 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=hclust(dist(scale(data)))
  cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  pvaluesmatrix
}


multilevel_kkmeans_clustering <- function(datalist, featurelist, numcluster, numf2dt,numf3dt,numf2dtnet, numf3dtnet,directionadd){
  days <- datalist[[1]]$days
  isdead <- datalist[[1]][,1] %in% deadsamples[,1]
  survobj <- Surv(days, isdead)
  pvaluesmatrix <- matrix(0,7,12)
  
  for(j in c(1:3,5:9,11)){
    for(i in 1:3){
      data <- datalist[[i]][,featurelist[[i]][[j]]+1]
      fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
      cl_kmeans = fit_cluster_kmeans@.Data
      fit <- survfit(survobj ~ cl_kmeans)
      sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
      pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
      # print(pval_kmeans)
      pvaluesmatrix[i,j]  <- pval_kmeans
    }
    i=4
    # combine GE and DM
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=5
    # combine GE and CNA
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=6
    # combine DM and CNA
    Data1 <- datalist[[2]][,featurelist[[2]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=7
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf3dt]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:numf3dt]+1]
    Data3 <- datalist[[3]][,featurelist[[3]][[j]][1:numf3dt]+1]
    data <- data.frame(Data1,Data2,Data3)
    fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  
  j=4
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[1]][x])
  dataadd1 <- sapply(featurelist[[1]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[2]][x])
  dataadd2 <- sapply(featurelist[[2]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[3]][x])
  dataadd3 <- sapply(featurelist[[3]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  j=10
  
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  dataadd1 <- sapply(featurelist[[1]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  dataadd2 <- sapply(featurelist[[2]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  dataadd3 <- sapply(featurelist[[3]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  
  j=12
  featurelist <- list(c(1:(ncol(datalist[[1]])-2)),c(1:(ncol(datalist[[2]])-2)),c(1:(ncol(datalist[[3]])-2)))
  for(i in 1:3){
    data <- datalist[[i]][,featurelist[[i]]+1]
    fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- datalist[[2]][,featurelist[[2]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  Data3 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=kkmeans(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  pvaluesmatrix
}

multilevel_spec_clustering <- function(datalist, featurelist, numcluster, numf2dt,numf3dt,numf2dtnet, numf3dtnet,directionadd){
  days <- datalist[[1]]$days
  isdead <- datalist[[1]][,1] %in% deadsamples[,1]
  survobj <- Surv(days, isdead)
  pvaluesmatrix <- matrix(0,7,12)
  for(j in c(1:3,5:9,11)){
    for(i in 1:3){
      data <- datalist[[i]][,featurelist[[i]][[j]]+1]
      fit_cluster_kmeans=specc(data.matrix(data),centers=3)
      cl_kmeans = fit_cluster_kmeans@.Data
      fit <- survfit(survobj ~ cl_kmeans)
      sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
      pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
      # print(pval_kmeans)
      pvaluesmatrix[i,j]  <- pval_kmeans
    }
    i=4
    # combine GE and DM
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:10]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:10]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=specc(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=5
    # combine GE and CNA
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:10]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:10]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=specc(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=6
    # combine DM and CNA
    Data1 <- datalist[[2]][,featurelist[[2]][[j]][1:10]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:10]+1]
    data <- data.frame(Data1,Data2)
    fit_cluster_kmeans=specc(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=7
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:7]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:7]+1]
    Data3 <- datalist[[3]][,featurelist[[3]][[j]][1:7]+1]
    data <- data.frame(Data1,Data2,Data3)
    fit_cluster_kmeans=specc(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  
  j=4
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[1]][x])
  dataadd1 <- sapply(featurelist[[1]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[2]][x])
  dataadd2 <- sapply(featurelist[[2]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[3]][x])
  dataadd3 <- sapply(featurelist[[3]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=specc(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  j=10
  
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  dataadd1 <- sapply(featurelist[[1]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  dataadd2 <- sapply(featurelist[[2]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  dataadd3 <- sapply(featurelist[[3]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=specc(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  
  
  j=12
  featurelist <- list(c(1:(ncol(datalist[[1]])-2)),c(1:(ncol(datalist[[2]])-2)),c(1:(ncol(datalist[[3]])-2)))
  for(i in 1:3){
    data <- datalist[[i]][,featurelist[[i]]+1]
    fit_cluster_kmeans=specc(data.matrix(data),centers=3)
    cl_kmeans = fit_cluster_kmeans@.Data
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- datalist[[2]][,featurelist[[2]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  Data3 <- datalist[[3]][,featurelist[[3]]+1]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=specc(data.matrix(data),centers=3)
  cl_kmeans = fit_cluster_kmeans@.Data
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  pvaluesmatrix
}

multilevel_icluster_clustering <- function(datalist, featurelist, numcluster, numf2dt,numf3dt,numf2dtnet, numf3dtnet,directionadd){
  days <- datalist[[1]]$days
  isdead <- datalist[[1]][,1] %in% deadsamples[,1]
  survobj <- Surv(days, isdead)
  pvaluesmatrix <- matrix(0,7,12)
  for(j in c(1:3,5:9,11)){
    for(i in 1:3){
      data <- datalist[[i]][,featurelist[[i]][[j]]+1]
      fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
      cl_kmeans = fit_cluster_kmeans$cluster
      fit <- survfit(survobj ~ cl_kmeans)
      sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
      pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
      # print(pval_kmeans)
      pvaluesmatrix[i,j]  <- pval_kmeans
    }
    i=4
    # combine GE and DM
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:10]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:10]+1]
    fit_cluster_kmeans=iCluster2(list(Data1,Data2),k=3)
    cl_kmeans = fit_cluster_kmeans$cluster
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=5
    # combine GE and CNA
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:10]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:10]+1]
    fit_cluster_kmeans=iCluster2(list(Data1,Data2),k=3)
    cl_kmeans = fit_cluster_kmeans$cluster
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=6
    # combine DM and CNA
    Data1 <- datalist[[2]][,featurelist[[2]][[j]][1:10]+1]
    Data2 <- datalist[[3]][,featurelist[[3]][[j]][1:10]+1]
    fit_cluster_kmeans=iCluster2(list(Data1,Data2),k=3)
    cl_kmeans = fit_cluster_kmeans$cluster
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    
    i=7
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:7]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:7]+1]
    Data3 <- datalist[[3]][,featurelist[[3]][[j]][1:7]+1]
    fit_cluster_kmeans=iCluster2(list(Data1,Data2,Data3),k=3)
    cl_kmeans = fit_cluster_kmeans$cluster
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  
  j=4
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[1]][x])
  dataadd1 <- sapply(featurelist[[1]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[2]][x])
  dataadd2 <- sapply(featurelist[[2]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[3]][x])
  dataadd3 <- sapply(featurelist[[3]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
    cl_kmeans = fit_cluster_kmeans$cluster
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  j=10
  
  moldata <- datalist[[1]][,-c(1,ncol(datalist[[1]]))]
  dataadd1 <- sapply(featurelist[[1]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[2]][,-c(1,ncol(datalist[[2]]))]
  dataadd2 <- sapply(featurelist[[2]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[3]][,-c(1,ncol(datalist[[3]]))]
  dataadd3 <- sapply(featurelist[[3]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2,dataadd3)
  
  for(i in 1:3){
    data <- dataaddlist[[i]]
    fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
    cl_kmeans = fit_cluster_kmeans$cluster
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- dataaddlist[[2]][,1:numf2dtnet]
  Data2 <- dataaddlist[[3]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- dataaddlist[[1]][,1:numf3dtnet]
  Data2 <- dataaddlist[[2]][,1:numf3dtnet]
  Data3 <- dataaddlist[[3]][,1:numf3dtnet]
  data <- data.frame(Data1,Data2,Data3)
  fit_cluster_kmeans=iCluster2(list(data.matrix(data)),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  

  j=12
  featurelist <- list(c(1:(ncol(datalist[[1]])-2)),c(1:(ncol(datalist[[2]])-2)),c(1:(ncol(datalist[[3]])-2)))
  for(i in 1:3){
    data <- datalist[[i]][,featurelist[[i]]+1]
    fit_cluster_kmeans=iCluster2(list(data),k=3)
    cl_kmeans = fit_cluster_kmeans$cluster
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
  }
  i=4
  # combine GE and DM
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  fit_cluster_kmeans=iCluster2(list(Data1,Data2),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=5
  # combine GE and CNA
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  fit_cluster_kmeans=iCluster2(list(Data1,Data2),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=6
  # combine DM and CNA
  Data1 <- datalist[[2]][,featurelist[[2]]+1]
  Data2 <- datalist[[3]][,featurelist[[3]]+1]
  fit_cluster_kmeans=iCluster2(list(Data1,Data2),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  
  i=7
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  Data3 <- datalist[[3]][,featurelist[[3]]+1]
  fit_cluster_kmeans=iCluster2(list(Data1,Data2,Data3),k=3)
  cl_kmeans = fit_cluster_kmeans$cluster
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  pvaluesmatrix
}


set.seed(1)
#kmeans
kmeansresults74_3cluster <- multilevel_kmeans_clustering(datatypes74_common3,listfeatures74,3,10,7,7,5,add_dire74)
kmeansresults74_3cluster <- as.data.frame(kmeansresults74_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(kmeansresults74_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

#kernal kmeans
kkmeansresults74_3cluster <- multilevel_kkmeans_clustering(datatypes74_common3,listfeatures74,3,10,7,7,5,add_dire74)
kkmeansresults74_3cluster <- as.data.frame(kkmeansresults74_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(kkmeansresults74_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

#spectral clustering
specresults74_3cluster <- multilevel_spec_clustering(datatypes74_common3,listfeatures74,3,10,7,7,5,add_dire74)
specresults74_3cluster <- as.data.frame(specresults74_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(specresults74_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

#SNF
SNFresults74_3cluster <- SNF_clustering(datatypes74_common3,listfeatures74,3,10,7,7,5,add_dire74)
SNFresults74_3cluster <- as.data.frame(SNFresults74_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(SNFresults74_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

#iCluster
iClusterresults74_3cluster <- multilevel_icluster_clustering(datatypes74_common3,listfeatures74,3,10,7,7,5,add_dire74)
iClusterresults74_3cluster <- as.data.frame(iClusterresults74_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(iClusterresults74_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")


set.seed(1)
kmeansresults123_3cluster <- multilevel_kmeans_clustering(datatypes123_common3,listfeatures123,3,10,7,7,5,add_dire123)
kmeansresults123_3cluster <- as.data.frame(kmeansresults123_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(kmeansresults123_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

kkmeansresults123_3cluster <- multilevel_kkmeans_clustering(datatypes123_common3,listfeatures123,3,10,7,7,5,add_dire123)
kkmeansresults123_3cluster <- as.data.frame(kkmeansresults123_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(kkmeansresults123_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

specresults123_3cluster <- multilevel_spec_clustering(datatypes123_common3,listfeatures123,3,10,7,7,5,add_dire123)
specresults123_3cluster <- as.data.frame(specresults123_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(specresults123_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

SNFresults123_3cluster <- SNF_clustering(datatypes123_common3,listfeatures123,3,10,7,7,5,add_dire123)
SNFresults123_3cluster <- as.data.frame(SNFresults123_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(SNFresults123_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

iClusterresults123_3cluster <- multilevel_icluster_clustering(datatypes123_common3,listfeatures123,3,10,7,7,5,add_dire123)
iClusterresults123_3cluster <- as.data.frame(iClusterresults123_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(iClusterresults123_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")


set.seed(1)
kmeansresults455_3cluster <- multilevel_kmeans_clustering(datatypes455_common3,listfeatures455,3,10,7,7,5,add_dire455)
kmeansresults455_3cluster <- as.data.frame(kmeansresults455_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(kmeansresults455_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

kkmeansresults455_3cluster <- multilevel_kkmeans_clustering(datatypes455_common3,listfeatures455,3,10,7,7,5,add_dire455)
kkmeansresults455_3cluster <- as.data.frame(kkmeansresults455_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(kkmeansresults455_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

specresults455_3cluster <- multilevel_spec_clustering(datatypes455_common3,listfeatures455,3,10,7,7,5,add_dire455)
specresults455_3cluster <- as.data.frame(specresults455_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(specresults455_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

SNFresults455_3cluster <- SNF_clustering(datatypes455_common3,listfeatures455,3,10,7,7,5,add_dire455)
SNFresults455_3cluster <- as.data.frame(SNFresults455_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(SNFresults455_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

iClusterresults455_3cluster <- multilevel_icluster_clustering(datatypes455_common3,listfeatures455,3,10,7,7,5,add_dire455) 
iClusterresults455_3cluster <- as.data.frame(iClusterresults455_3cluster, row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(iClusterresults455_3cluster) <-  c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble","allemt")

save.image('LUAD_with_rules.RData')
# use SNF, iCluster, kmeans, and kkmeans on all features
# list(allrnaseq,allmethy,allcna)



