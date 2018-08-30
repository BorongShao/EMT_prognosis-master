# binary class
netLasso <- function(trainset){
  adjmatrix <- get.adjacency(LEnet,sparse = FALSE)
  L <- make.L(adjmatrix,normalize.Laplacian=FALSE)
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset[,ncol(trainset)]+0
  fit.grace <- grace(y, data.matrix(X), L, lambda.L = seq(0,1,0.1), lambda.1=1,lambda.2=0.01,normalize.L=FALSE, K=10, verbose=FALSE)
  nonzero <- which(fit.grace[[2]]!=0)
  nonzero[order(abs(fit.grace[[2]][nonzero]),decreasing = TRUE)]
}

# single_level_methods <- function(thresdata, nothresdata, threstrainidx, nothrestrainidx){
#   ntopgene <- 25
#   ntopnet <- 10
#   trainset <- LEthresgeneexp[threstrainidx,]
#   testset <- LEthresgeneexp[-threstrainidx,]
# 
#   #netowrk lasso
#   beta <- netLasso(trainset)
#   result_lasso <- sgclassvec(beta,testset)
# }
# 
# registerDoMC(cores=10)
# temp <- foreach(x = 1:100)  %dopar% {
#   threstrainidx <- sample(nrow(LEthresgeneexp),floor(nrow(LEthresgeneexp)*0.8))
#   nothrestrainidx <- sample(nrow(LEgeneexp),floor(nrow(LEgeneexp)*0.8))
#   single_level_methods(LEthresgeneexp,LEgeneexp,threstrainidx,nothrestrainidx)
# }

