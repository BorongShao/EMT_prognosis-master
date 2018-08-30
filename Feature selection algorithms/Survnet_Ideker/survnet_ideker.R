scoring <- function(data,survobj) {
  #p value from univariate cox model
  pvalues <- sapply(data, function(x) coef(summary(coxph(survobj ~ x)))[5])
  # z-score CDF
  zscores <- qnorm(1-pvalues)
  zscores
}

survnet <- function(trainset,ntopnet) {
  data <- trainset[,-c(1,ncol(trainset))]
  survobj <- Surv(trainset$days, trainset[,1] %in% deadsamples[,1])
  zscores <- scoring(data,survobj)
  scores <- c()
  subnetworks <- list()
  nbins <- floor(log2(nrow(data))+1)
  #start from each node
  for(i in 1:ncol(data)){
    snode = i
    nb <- as.numeric(neighbors(LEnet,snode,mode ='all'))
    if(length(nb)==0)
      next
    avg <- sapply(nb, function(x) {
      vec <- c(snode,x)
      sum(zscores[vec])/sqrt(length(vec))
    })
    idx <- nb[which.max(avg)]
    vec <- c(snode,idx)
    if((sum(zscores[vec])/sqrt(length(vec))) >  (sum(zscores[snode])/sqrt(length(snode)))*1.05){
      snode <- c(snode,idx)
      temp <- extd_survnet(data,snode,zscores)
      subnetworks[[length(subnetworks)+1]] <- temp[[1]]
      scores <- c(scores,temp[[2]])
    }
  }
  scores1 <-score_calibrate(subnetworks,zscores,scores)
  num <- ifelse(length(subnetworks)>ntopnet, ntopnet, length(subnetworks))
  subnetidx <- order(scores1,decreasing = TRUE)[1:num]
  list(subnetworks[subnetidx])
}

survnet_verticle <- function(trainset,ntopnet,interedges) {
  data <- trainset[,-c(1,ncol(trainset))]
  survobj <- Surv(trainset$days, trainset[,1] %in% deadsamples[,1])
  zscores <- scoring(data,survobj)
  scores <- c()
  subnetworks <- list()
  nbins <- floor(log2(nrow(data))+1)
  #start from each node
  for(i in 1:ncol(data)){
    snode = i
    nb <- as.numeric(neighbors(LEnet,snode,mode ='all'))
    # add the nodes for interedges
    vnode <- intersect(interedges,nb)
    nodestoadd <- c()
    for(i in 1:length(vnode)){
      if(i>455)
        nodestoadd <- c(nodestoadd,i-455) else
          nodestoadd <- c(nodestoadd,i+455)
    }
    nb <- unique(c(nb,nodestoadd))
    if(length(nb)==0)
      next
    avg <- sapply(nb, function(x) {
      vec <- c(snode,x)
      sum(zscores[vec])/sqrt(length(vec))
    })
    idx <- nb[which.max(avg)]
    vec <- c(snode,idx)
    if((sum(zscores[vec])/sqrt(length(vec))) >  (sum(zscores[snode])/sqrt(length(snode)))*1.05){
      snode <- c(snode,idx)
      temp <- extd_survnet_verticle(data,snode,zscores,interedges)
      subnetworks[[length(subnetworks)+1]] <- temp[[1]]
      scores <- c(scores,temp[[2]])
    }
  }
  scores1 <-score_calibrate(subnetworks,zscores,scores)
  num <- ifelse(length(subnetworks)>ntopnet, ntopnet, length(subnetworks))
  subnetidx <- order(scores1,decreasing = TRUE)[1:num]
  list(subnetworks[subnetidx])
}


survnetseeds <- function(trainset,ntopnet,interedges) {
  data <- trainset[,-c(1,ncol(trainset))]
  survobj <- Surv(trainset$days, trainset[,1] %in% deadsamples[,1])
  zscores <- scoring(data,survobj)
  scores <- c()
  subnetworks <- list()
  nbins <- floor(log2(nrow(data))+1)
  #start from each node
  for(i in interedges){
    snode = i
    nb <- as.numeric(neighbors(LEnet,snode,mode ='all'))
    if(length(nb)==0)
      next
    avg <- sapply(nb, function(x) {
      vec <- c(snode,x)
      sum(zscores[vec])/sqrt(length(vec))
    })
    idx <- nb[which.max(avg)]
    vec <- c(snode,idx)
    if((sum(zscores[vec])/sqrt(length(vec))) >  (sum(zscores[snode])/sqrt(length(snode)))*1.05){
      snode <- c(snode,idx)
      temp <- extd_survnet(data,snode,zscores)
      subnetworks[[length(subnetworks)+1]] <- temp[[1]]
      scores <- c(scores,temp[[2]])
    }
  }
  scores1 <-score_calibrate(subnetworks,zscores,scores)
  num <- ifelse(length(subnetworks)>ntopnet, ntopnet, length(subnetworks))
  subnetidx <- order(scores1,decreasing = TRUE)[1:num]
  list(subnetworks[subnetidx])
}

extd_survnet <- function(data,snode,zscores) {
  nbins <- floor(log2(nrow(data))+1)
  temp <- unlist(sapply(snode,function(x) as.numeric(neighbors(LEnet,x,mode ='all'))))
  nb <- unique(unlist(temp[!(temp %in% snode)]))
  avg <- sapply(nb, function(x) {
    vec <- c(snode,x)
    sum(zscores[vec])/sqrt(length(vec))
  })
  idx <- nb[which.max(avg)]
  vec <- c(snode,idx)
  if((sum(zscores[vec])/sqrt(length(vec))) > (sum(zscores[snode])/sqrt(length(snode)))*1.05){
    snode <- c(snode,idx)
    extd_survnet(data,snode,zscores)
  }
  else{
    return(list(snode,sum(zscores[vec])/sqrt(length(vec))))
  }
}

extd_survnet_verticle <- function(data,snode,zscores,interedges) {
  nbins <- floor(log2(nrow(data))+1)
  temp <- unlist(sapply(snode,function(x) as.numeric(neighbors(LEnet,x,mode ='all'))))
  nb <- unique(unlist(temp[!(temp %in% snode)]))
  # add the nodes for interedges
  vnode <- intersect(interedges,nb)
  nodestoadd <- c()
  for(i in 1:length(vnode)){
    if(i>455)
      nodestoadd <- c(nodestoadd,i-455) else
        nodestoadd <- c(nodestoadd,i+455)
  }
  nb <- unique(c(nb,nodestoadd))
  nb <- nb[!(nb %in% snode)]
  
  avg <- sapply(nb, function(x) {
    vec <- c(snode,x)
    sum(zscores[vec])/sqrt(length(vec))
  })
  idx <- nb[which.max(avg)]
  vec <- c(snode,idx)
  if((sum(zscores[vec])/sqrt(length(vec))) > (sum(zscores[snode])/sqrt(length(snode)))*1.05){
    snode <- c(snode,idx)
    extd_survnet_verticle(data,snode,zscores,interedges)
  }
  else{
    return(list(snode,sum(zscores[vec])/sqrt(length(vec))))
  }
}


score_calibrate <- function(subnetworks, zscores, scores){
  nfeatures <- ncol(data)
  lengthvec <- unique(sapply(subnetworks,length))
  meanvec <- rep(0,length(lengthvec))
  stdvec <- rep(0,length(lengthvec))
  for(i in 1:length(lengthvec)){
    l <- lengthvec[i]
    temp <- replicate(10000, {sum(zscores[sample(nfeatures,l)])/sqrt(l)})
    meanvec[i] <- mean(temp)
    stdvec[i] <- sd(temp)
  }
  for(i in 1:length(scores)){
    t <- which(lengthvec == length(subnetworks[[i]]))
    scores[i] <- (scores[i] - meanvec[t])/stdvec[t]
  }
  scores
}

#use data without thresholds
# y <- LEgeneexp$days
# data <- LEgeneexp[,-c(1,ncol(LEgeneexp))]
# survobj <- Surv(y, LEgeneexp$samplecode %in% deadsamples$samplecode)
# randomdata <- data[,sample(ncol(data))]
# 
# subnetworks <- list()
# scores <- c()
# zscores <- scoring(data,survobj)
# survnet(data)
# subnetworks1 <- subnetworks
# scores1 <- scores
# scores1.1 <-score_calibrate(subnetworks1,zscores, scores1)
# 
# subnetworks <- list()
# scores <- c()
# #recalculate zscores
# zscores <- scoring(randomdata,survobj)
# survnet(randomdata)
# subnetworks2 <- subnetworks
# scores2 <-scores
# scores2.1 <-score_calibrate(subnetworks2,zscores,scores2)
# 





