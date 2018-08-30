stSVM <- function(trainset, ntopgene){
  adjmatrix <- get.adjacency(LEnet,sparse = FALSE)
  dk <- calc.diffusionKernelp(L=adjmatrix, is.adjacency=TRUE,p=2,a=1)
  sca=TRUE			
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset$days
  vals <- sapply(X, function(x) {
    x1 <- x[y==0]
    x2 <- x[y==1]
    t.test(x1,x2)$statistic
  })
  vals <- abs(vals)/max(abs(vals))
  ranks = (t(vals) %*% dk)[1,]
  order(ranks,decreasing = TRUE)[1:ntopgene]
}


calc.diffusionKernelp = function(L, is.adjacency=TRUE, p=3,a=2)
{
  if(missing(L)) stop('You have to provide either a laplace or a adjacency matrix to calculate the diffusion kernel.')
  print("thresh")
  
  if(is.adjacency)
  {    
    dnames <- dimnames(L)
    L = graph.adjacency(L, mode='undirected', diag=FALSE)
    L = graph.laplacian(L, normalized=TRUE)
    dimnames(L) = dnames
  }
  
  n=ncol(L)
  I=diag(n)      
  
  if( p==1) R = a*I-L
  
  else
  {
    R=a*I -L
    for(ii in 2:p)
    {
      R = R %*% (a*I -L)
      ii = ii+1
    }
  }
  
  #KernelDiag <- sqrt(diag(R) + 1e-10)
  #R.norm 	<-  R/(KernelDiag %*% t(KernelDiag))
  #R			<-  R.norm  
  
  colnames(R) =colnames(L)
  rownames(R) =rownames(L)
  R    
}
