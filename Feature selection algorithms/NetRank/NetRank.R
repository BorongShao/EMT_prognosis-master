#LEnet
#training data without samplecode, class-based

NetRank <- function(trainset,ntopgene) {
  dvec <- seq(0,0.9,0.1)
  resultsd <- rep(0,length(dvec))
  adjmatrix <- as_adjacency_matrix(LEnet,sparse=FALSE)
  degrees = colSums(adjmatrix)
  degrees[which(degrees==0)] <- 1
  D1 <- diag(1./degrees)
  # try different d values
  for(i in 1:length(dvec)){
    print(i)
    d <- dvec[i]
    # MCCV 
   X <- trainset[,-c(1,ncol(trainset))]    
   y <- trainset$days
    aucvals <- rep(0,20)
    auprvals <- rep(0,20)
    #only count the results when converged, until 100 times
    for(j in 1:20){
      trainindex <- sample(nrow(X),floor(nrow(X)*0.9))
      traininner <- data.frame(X[trainindex,])
      testinner <- data.frame(X[-trainindex,])
      ex = abs(cor(traininner,y[trainindex],method = "spearman"))
      norm_ex = ex/max(ex)
      A = diag(length(ex)) - d*(adjmatrix%*%D1)
      b = (1-d)*norm_ex
      r = qr.solve(A,b)
      idx <- c(order(r,decreasing = TRUE)[1:ntopgene])
      #classification accuracy with selected genes
      svmresult <- svm1(traininner[,idx],y[trainindex],testinner[,idx],y[-trainindex])
      aucvals[j] <- svmresult[1]
      auprvals[j] <- svmresult[2]
    }
    resultsd[i] <- mean(unlist(aucvals))
  }
  #return the genelist with the best d value
  bestd <- dvec[which.max(resultsd)]
  ex = abs(cor(X,y,method = "spearman"))
  norm_ex = ex/max(ex)
  degrees = colSums(adjmatrix)
  degrees[which(degrees==0)] <- 1
  D1 <- diag(1./degrees)
  A = diag(length(ex)) - bestd*(adjmatrix%*%D1)
  b = (1-bestd)*norm_ex
  r = qr.solve(A,b)
  order(r,decreasing = TRUE)[1:ntopgene]
}

# 
# X <- trainset[,-ncol(trainset)]
# y <- trainset$days
# adjmatrix <- as_adjacency_matrix(LEnet)
# ex = abs(cor(X,y,method="spearman"))
# norm_ex = ex/max(ex)
# degrees = colSums(adjmatrix)
# D1 <- diag(1./degrees)
# A = diag(length(ex)) - d*(adjmatrix%*%D1)
# b = (1-d)*norm_ex
# r = qr.solve(A,b)
# idx <- c(order(r,decreasing = TRUE)[1:ntopgene],ncol(trainset))
# resultsd[i] <- SVM_MC_CV(trainset)
