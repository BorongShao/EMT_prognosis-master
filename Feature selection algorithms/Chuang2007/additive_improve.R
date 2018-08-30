
#experiment with search distance = 1 and 2

additive_subnet <- function(trainset,ntopnet) {
  data <- trainset[,-c(1,ncol(trainset))]
  y <- trainset$days
  add_direction <- sapply(data,function(x) sign(cor(x,y)))
  scores <- c()
  subnetworks <- list()
  nbins <- floor(log2(nrow(data))+1)
  #start from each node
  for(i in 1:ncol(data)){
    snode = i
    # print(i)
    nb <- as.numeric(neighbors(LEnet,snode,mode ='all'))
    if(length(nb)==0)
      next
    avg <- sapply(1:length(nb), function(x) {
      vec <- data[,snode]*add_direction[snode] + data[,nb[x]]*add_direction[nb[x]]
      mi(vec,y,nbins)
    })
    idx <- nb[which.max(avg)]
    vec <- data[,snode]*add_direction[snode] + data[,idx]*add_direction[idx]
    if(mi(vec,y,nbins) > mi(data[,snode]*add_direction[snode],y,nbins)*1.05){
      snode <- c(snode,idx)
      temp <- extd_subnet(data,y,snode)
      subnetworks[[length(subnetworks)+1]] <- temp[[1]]
      scores <- c(scores,temp[[2]])
    }
  }
  # subnetidx <- order(scores,decreasing = TRUE)[1:ntopnet]
  subnetidx <- order(scores,decreasing = TRUE)
  # netrank_addv <-list(subnetworks[subnetidx],scores[subnetidx])
  # list(scores,subnetworks)
  # subnetworks[subnetidx]
  uniquenets <- unique(sapply(subnetworks[subnetidx], function(x) sort(x)))
  num <- 50
  if(length(uniquenets) < 50)
    num <- length(uniquenets)
  # traindata <- sapply(uniquenets,function(x) rowSums(trainset[,(x+1)])/length(x))
  # trainy <- trainset$days
  # cvfit <- cv.glmnet(data.matrix(traindata),trainy,family = "binomial",type.measure = "auc")
  # coefficients <- coef(cvfit,s=c(cvfit$lambda.min))
  # metagenes <- which(coefficients!=0)-1
  # uniquenets[metagenes]
  uniquenets[1:num]
}
# select # features by lasso
extd_subnet <- function(data,y,snode) {
  add_direction <- sapply(data,function(x) sign(cor(x,y)))
  nbins <- floor(log2(nrow(data))+1)
  temp <- unlist(sapply(snode,function(x) as.numeric(neighbors(LEnet,x,mode ='all'))))
  # temp1 <- unlist(sapply(snode,function(x) as.numeric(neighbors(LEnet,x,mode ='all'))))
  # temp2 <- unlist(sapply(temp1,function(x) as.numeric(neighbors(LEnet,x,mode ='all'))))
  # temp <- c(temp1,temp2)
  nb <- unique(unlist(temp[!(temp %in% snode)]))
  veccurrent <- rep(0,nrow(data))
  for(i in snode){
    veccurrent  <- veccurrent + data[,i]*add_direction[i]
  }
  avg <- sapply(1:length(nb), function(x) {
    vec <- veccurrent + data[,nb[x]]*add_direction[nb[x]]  
    mi(vec,y,nbins)
  })  
  idx <- nb[which.max(avg)]
  vec <- veccurrent + data[,idx]*add_direction[idx]
  vecmi <- mi(vec,y,nbins)
  if(vecmi > mi(veccurrent,y,nbins)*1.05){
    snode <- c(snode,idx)
    extd_subnet(data,y,snode)
  }
  else{
    return(list(snode,mi(veccurrent,y,nbins)))
  }
}

additive_subnetseeds <- function(trainset,ntopnet,interedges) {
  data <- trainset[,-c(1,ncol(trainset))]
  y <- trainset$days
  add_direction <- sapply(data,function(x) sign(cor(x,y)))
  scores <- c()
  subnetworks <- list()
  nbins <- floor(log2(nrow(data))+1)
  #start from each node
  for(i in interedges){
    snode = i
    # print(i)
    nb <- as.numeric(neighbors(LEnet,snode,mode ='all'))
    if(length(nb)==0)
      next
    avg <- sapply(1:length(nb), function(x) {
      vec <- data[,snode]*add_direction[snode] + data[,nb[x]]*add_direction[nb[x]]
      mi(vec,y,nbins)
    })
    idx <- nb[which.max(avg)]
    vec <- data[,snode]*add_direction[snode] + data[,idx]*add_direction[idx]
    if(mi(vec,y,nbins) > mi(data[,snode]*add_direction[snode],y,nbins)*1.05){
      snode <- c(snode,idx)
      temp <- extd_subnet(data,y,snode)
      subnetworks[[length(subnetworks)+1]] <- temp[[1]]
      scores <- c(scores,temp[[2]])
    }
  }
  # subnetidx <- order(scores,decreasing = TRUE)[1:ntopnet]
  subnetidx <- order(scores,decreasing = TRUE)
  # netrank_addv <-list(subnetworks[subnetidx],scores[subnetidx])
  # list(scores,subnetworks)
  # subnetworks[subnetidx]
  uniquenets <- unique(sapply(subnetworks[subnetidx], function(x) sort(x)))
  num <- 50
  if(length(uniquenets) < 50)
    num <- length(uniquenets)
  # traindata <- sapply(uniquenets,function(x) rowSums(trainset[,(x+1)])/length(x))
  # trainy <- trainset$days
  # cvfit <- cv.glmnet(data.matrix(traindata),trainy,family = "binomial",type.measure = "auc")
  # coefficients <- coef(cvfit,s=c(cvfit$lambda.min))
  # metagenes <- which(coefficients!=0)-1
  # uniquenets[metagenes]
  uniquenets[1:num]
}

additive_subnet_verticle <- function(trainset,ntopnet,interedges) {
  data <- trainset[,-c(1,ncol(trainset))]
  y <- trainset$days
  add_direction <- sapply(data,function(x) sign(cor(x,y)))
  scores <- c()
  subnetworks <- list()
  nbins <- floor(log2(nrow(data))+1)
  #start from each node
  for(i in 1:ncol(data)){
    snode = i
    # print(i)
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
    avg <- sapply(1:length(nb), function(x) {
      vec <- data[,snode]*add_direction[snode] + data[,nb[x]]*add_direction[nb[x]]
      mi(vec,y,nbins)
    })
    idx <- nb[which.max(avg)]
    vec <- data[,snode]*add_direction[snode] + data[,idx]*add_direction[idx]
    if(mi(vec,y,nbins) > mi(data[,snode]*add_direction[snode],y,nbins)*1.05){
      snode <- c(snode,idx)
      temp <- extd_subnet_verticle(data,y,snode,interedges)
      subnetworks[[length(subnetworks)+1]] <- temp[[1]]
      scores <- c(scores,temp[[2]])
    }
  }
  # subnetidx <- order(scores,decreasing = TRUE)[1:ntopnet]
  subnetidx <- order(scores,decreasing = TRUE)
  # netrank_addv <-list(subnetworks[subnetidx],scores[subnetidx])
  # list(scores,subnetworks)
  # subnetworks[subnetidx]
  uniquenets <- unique(sapply(subnetworks[subnetidx], function(x) sort(x)))
  num <- 50
  if(length(uniquenets) < 50)
    num <- length(uniquenets)
  # traindata <- sapply(uniquenets,function(x) rowSums(trainset[,(x+1)])/length(x))
  # trainy <- trainset$days
  # cvfit <- cv.glmnet(data.matrix(traindata),trainy,family = "binomial",type.measure = "auc")
  # coefficients <- coef(cvfit,s=c(cvfit$lambda.min))
  # metagenes <- which(coefficients!=0)-1
  # uniquenets[metagenes]
  uniquenets[1:num]
}

extd_subnet_verticle <- function(data,y,snode,interedges) {
  add_direction <- sapply(data,function(x) sign(cor(x,y)))
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
  
  veccurrent <- rep(0,nrow(data))
  for(i in snode){
    veccurrent  <- veccurrent + data[,i]*add_direction[i]
  }
  avg <- sapply(1:length(nb), function(x) {
    vec <- veccurrent + data[,nb[x]]*add_direction[nb[x]]  
    mi(vec,y,nbins)
  })  
  idx <- nb[which.max(avg)]
  vec <- veccurrent + data[,idx]*add_direction[idx]
  vecmi <- mi(vec,y,nbins)
  if(vecmi > mi(veccurrent,y,nbins)*1.05){
    snode <- c(snode,idx)
    extd_subnet_verticle(data,y,snode,interedges)
  }
  else{
    return(list(snode,mi(veccurrent,y,nbins)))
  }
}


# y <- LEthresgeneexp$days
# data <- LEthresgeneexp[,-c(1,ncol(LEthresgeneexp))]
# randomdata <- data[,sample(ncol(data))]



# subnetworks <- list()
# scores <- c()
# additive_subnet(randomdata,y)
# subnetworks2 <- subnetworks
# scores2 <-scores

