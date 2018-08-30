# t-test
# return a list of geneindex based on training data
# data does not have sample column, with class label, return top n genes.
ttestgenes <- function(trainset,ntopgene){
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset[,ncol(trainset)]
  vals <- sapply(X, function(x) {
    x1 <- x[y==0]
    x2 <- x[y==1]
    wilcox.test(x1,x2)$p.value
  })
  order(vals)[1:ntopgene]
}

corrankgenes <- function(trainset,ntopgene){
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset[,ncol(trainset)]
  vals <- sapply(X, function(x) abs(cor(x,y,method = "spearman")))
  order(vals)[1:ntopgene]
}


#return lasso coefficients based on training data
# data does not have sample column, with class label
lassocoefs <- function(trainset){
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset[,ncol(trainset)]
  coef <- NULL
  count <- 0
  while(TRUE){
    count <- count + 1
    cvfit <- cv.glmnet(data.matrix(X),y,family = "binomial",type.measure = "auc")
    coefficients <- as.numeric(coef(cvfit,s=c(cvfit$lambda.min)))
    coef <- coefficients
    if(length(which(coef!=0)) > 1)
      break
    if(count > 50){
      prtin('bad data')
      coef <- rep(0,length(coef))
      break
    }
  }
  nonzero <- which(coef[-1]!=0)
  nonzero[order(abs(coef[nonzero+1]), decreasing = TRUE)]
}

#return lasso coeffieicnets with cox variable based on training data
# data have sample column, 
lassocoxcoefs <- function(trainset){
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset$days
  survobj <- Surv(y, trainset[,1] %in% deadsamples[,1])
  coef <- NULL
  count <- 0
  while(TRUE){
    count <- count + 1
    cvfit <- cv.glmnet(data.matrix(X),survobj,family="cox",maxit=100000,alpha=0.95)
    coefficients <- as.numeric(coef(cvfit,s=c(cvfit$lambda.min)))
    coef <- coefficients
    if(length(which(coef!=0)) > 1)
      break
    if(count > 50){
      coef <- rep(0,length(coef))
      break
    }
  }
  coef
}

#ridge regression
logistic_reg <- function(trainset) {
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset[,ncol(trainset)]
  cvfit <- glm.fit(data.matrix(X),y,family = binomial(link = "logit"),control = list(maxit = 100))
  coefficients <- as.numeric(coef(cvfit,s=c(cvfit$lambda.min)))
  coefficients
}



# rank features based on p values, dataset with sample column
Coxrank <- function(trainset,ntopgene){
  isdead <- trainset[,1] %in% deadsamples[,1]
  days <- trainset$days
  X <- trainset[,-c(1,ncol(trainset))]
  survobj <- Surv(days, isdead)
  pvalues <- sapply(X, function(x) coef(summary(coxph(survobj ~ x)))[5])
  order(pvalues)[1:ntopgene]
}

SVM_MC_CV <- function(data){
  y <- data$days > thresupper
  aucvals <- rep(0,100)
  auprvals <- rep(0,100)
  for(i in 1:100){
    trainindex <- sample(nrow(data),floor(nrow(data)*0.9))
    train <- data.frame(data[trainindex,-c(ncol(data))])
    test <- data.frame(data[-trainindex,-c(ncol(data))])
    rbf <- accusvm_RBF_default(train,test,y,trainindex)
    aucvals[i] <- rbf[1]
    auprvals[i] <- rbf[2]
    print(i)
  }
  # list(aucvals,auprvals)
  mean(aucvals)
}

lasso_MC_CV <- function(data){
  y <- data$days > thresupper
  aucvals <- rep(0,100)
  auprvals <- rep(0,100)
  #only count the results when converged, until 100 times
  for(i in 1:100){
    trainindex <- sample(nrow(data),floor(nrow(data)*0.9))
    train <- data.frame(data[trainindex,-c(1,ncol(data))])
    test <- data.frame(data[-trainindex,-c(1,ncol(data))])
    ls <- FSbylasso(train,test,y,trainindex)
    aucvals[i] <- ls[1]
    auprvals[i] <- ls[2]
    print(i)
  }
  list(aucvals,auprvals)
}

mi <- function (a,b,nbins=2) {
  #disb <- discretize(b,nbins=2)
  re <- NULL
  disb=b
  if(is.vector(a)==FALSE){
    col <- ncol(a)
    re <- rep(1,col)
    for(i in 1:col){
      disa <- discretize(a[,i],disc='equalwidth',nbins=nbins)
      re[i] <- mutinformation(disa,disb)
    }
  }else{
    dis <- infotheo::discretize(a,disc='equalwidth',nbins=nbins)
    re <- mutinformation(dis,b)
  }
  re
}
# trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days
svm1 <- function(traindata,trainy,testdata,testy){
  X <- traindata
  y <- trainy
  # svm.model <- svm(X,as.factor(y),probability=TRUE, gamma = 2^-7, cost = 2)
  # svm.model <- svm(X,as.factor(y),probability=TRUE, kernel="linear")
  svm.model <- svm(X,as.factor(y),probability=TRUE)
  svm.pred <- predict(svm.model, testdata,probability=TRUE)
  if(trainy[1]==TRUE){
    predvec <- attr(svm.pred,"prob")[,1] 
  }else{
    predvec <- attr(svm.pred,"prob")[,2] 
  }
  # pred <- prediction(predvec,as.numeric(testy))
  roc_obj <- roc(as.numeric(testy), attr(svm.pred,"prob")[,1])
  coordinatesacc <- as.numeric(coords(roc_obj, seq(0,1,0.01),input="threshold", ret="accuracy"))
  # coordinatessens <- as.numeric(coords(roc_obj, seq(1,0,-0.01),input="specificity", ret="sensitivity"))
  # auc <- performance(pred,"auc")
  # auc <- roc.curve(scores.class0 = svm.pred, weights.class0 = y)$auc
  aupr <- pr.curve(scores.class0 = predvec, weights.class0 = testy)$auc.integral
  # auc@y.values
  
  list(as.numeric(roc_obj$auc),aupr, coordinatesacc)
}
# traindata,trainy,testdata,testy,testset2
svm_cox <- function(traindata,trainy,testdata,testy,testset2){
  X <- traindata
  y <- trainy
  # svm.model <- svm(X,as.factor(y),probability=TRUE, gamma = 2^-7, cost = 2)
  # svm.model <- svm(X,as.factor(y),probability=TRUE, kernel="linear")
  svm.model <- svm(X,as.factor(y),probability=TRUE)
  svm.pred <- predict(svm.model, testdata,probability=TRUE)
  if(trainy[1]==TRUE){
    predvec <- attr(svm.pred,"prob")[,1] 
  }else{
    predvec <- attr(svm.pred,"prob")[,2] 
  }
  # pred <- prediction(predvec,as.numeric(testy))
  roc_obj <- roc(as.numeric(testy), attr(svm.pred,"prob")[,1])
  coordinatesacc <- as.numeric(coords(roc_obj, seq(0,1,0.01),input="threshold", ret="accuracy"))
  # coordinatessens <- as.numeric(coords(roc_obj, seq(1,0,-0.01),input="specificity", ret="sensitivity"))
  # auc <- performance(pred,"auc")
  # auc <- roc.curve(scores.class0 = svm.pred, weights.class0 = y)$auc
  aupr <- pr.curve(scores.class0 = predvec, weights.class0 = testy)$auc.integral
  # auc@y.values
  svm.pred2 <- predict(svm.model, testset2[,-c(1,ncol(testset2))], probability=TRUE)
  if(trainy[1]==TRUE){
    predvec <- attr(svm.pred2,"prob")[,1] 
  }else{
    predvec <- attr(svm.pred2,"prob")[,2] 
  }
  temp <- sgsurvivalvec(testset2,predvec)
  list(as.numeric(roc_obj$auc),aupr, coordinatesacc, temp[[1]],temp[[2]])
}



randomforest1 <- function(traindata,trainy,testdata,testy){
  ranfor <- randomForest(x=traindata,y=as.factor(trainy))
  rf.pred <- predict(ranfor,testdata,type="prob")
  roc_obj <- roc(as.numeric(testy), rf.pred[,2])
  coordinatesacc <- as.numeric(coords(roc_obj, seq(0,1,0.01),input="threshold", ret="accuracy"))
  coordinatessens <- as.numeric(coords(roc_obj, seq(1,0,-0.01),input="specificity", ret="sensitivity"))
  aupr <- pr.curve(scores.class0 = rf.pred[,2], weights.class0 = testy)$auc.integral
  list(as.numeric(roc_obj$auc),aupr, coordinatesacc,coordinatessens)
}


randomforest_cox <- function(traindata,trainy,testdata,testy,testset2){
  ranfor <- randomForest(x=traindata,y=as.factor(trainy))
  rf.pred <- predict(ranfor,testdata,type="prob")
  roc_obj <- roc(as.numeric(testy), rf.pred[,2])
  coordinatesacc <- as.numeric(coords(roc_obj, seq(0,1,0.01),input="threshold", ret="accuracy"))
  # coordinatessens <- as.numeric(coords(roc_obj, seq(1,0,-0.01),input="specificity", ret="sensitivity"))
  aupr <- pr.curve(scores.class0 = rf.pred[,2], weights.class0 = testy)$auc.integral
  testdata2 <- matrix(unlist(data.matrix(testset2[,-c(1,ncol(testset2))])),ncol=ncol(traindata))

  rf.pred2=predict(ranfor,testdata2,type="prob")[,2]
  temp <- sgsurvivalvec(testset2,rf.pred2)
  list(as.numeric(roc_obj$auc),aupr, coordinatesacc, temp[[1]],temp[[2]])
}


logisticreg1 <- function(traindata,trainy,testdata,testy){
  temp <- data.frame(traindata,trainy)
  fit <- glm(trainy ~ ., data = temp, family = binomial(link = "logit"),control = list(maxit = 100))
  testdata <- as.data.frame(testdata)
  names(testdata) <- names(temp)[-ncol(temp)]
  lr.pred=predict(fit,testdata,type="response")
  roc_obj <- roc(as.numeric(testy), lr.pred)
  coordinatesacc <- as.numeric(coords(roc_obj, seq(0,1,0.01),input="threshold", ret="accuracy"))
  coordinatessens <- as.numeric(coords(roc_obj, seq(1,0,-0.01),input="specificity", ret="sensitivity"))
  aupr <- pr.curve(scores.class0 = lr.pred, weights.class0 = testy)$auc.integral
  list(as.numeric(roc_obj$auc),aupr, coordinatesacc,coordinatessens)
}

logisticreg_cox <- function(traindata,trainy,testdata,testy,testset2){
  temp <- data.frame(traindata,trainy)
  fit <- glm(trainy ~ ., data = temp, family = binomial(link = "logit"),control = list(maxit = 100))
  testdata <- as.data.frame(testdata)
  names(testdata) <- names(temp)[-ncol(temp)]
  lr.pred=predict(fit,testdata,type="response")
  roc_obj <- roc(as.numeric(testy), lr.pred)
  coordinatesacc <- as.numeric(coords(roc_obj, seq(0,1,0.01),input="threshold", ret="accuracy"))
  # coordinatessens <- as.numeric(coords(roc_obj, seq(1,0,-0.01),input="specificity", ret="sensitivity"))
  aupr <- pr.curve(scores.class0 = lr.pred, weights.class0 = testy)$auc.integral
  
  lr.pred2=predict(fit,testset2[,-c(1,ncol(testset2))],type="response")
  temp <- sgsurvivalvec(testset2,lr.pred2)
  list(as.numeric(roc_obj$auc),aupr, coordinatesacc, temp[[1]],temp[[2]])
}


deeplearning1 <- function(traindata,trainy,testdata,testy){
  temp1 <- data.frame(traindata,trainy)
  write.csv(temp1, file="training.csv",row.names = FALSE)
  traindata <- h2o.importFile("training.csv")
  temp2 <- data.frame(testdata,testy)
  names(temp2) <- names(temp1)
  write.csv(temp2, file="testing.csv",row.names = FALSE)
  testdata <- h2o.importFile("testing.csv")
  dl <- h2o.deeplearning(x=1:118, y="trainy", training_frame = traindata, quiet_mode = TRUE, hidden=c(60,30))
  dn.pred <- as.vector(h2o.predict(dl, testdata)[[3]])
  roc_obj <- roc(as.vector(testdata$trainy), dn.pred)
  coordinatesacc <- as.numeric(coords(roc_obj, seq(0,1,0.01),input="threshold", ret="accuracy"))
  coordinatessens <- as.numeric(coords(roc_obj, seq(1,0,-0.01),input="specificity", ret="sensitivity"))
  aupr <- pr.curve(scores.class0 = dn.pred, weights.class0 = testy)$auc.integral
  list(as.numeric(roc_obj$auc),aupr, coordinatesacc,coordinatessens)
}

deeplearning2 <- function(succeedfolds){
  data <- mx.symbol.Variable("data")
  fc1 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=60)
  act1 <- mx.symbol.Activation(fc1, name="relu1", act_type="sigmoid")
  fc2 <- mx.symbol.FullyConnected(act1, name="fc2", num_hidden=30)
  act2 <- mx.symbol.Activation(fc2, name="relu2", act_type="sigmoid")
  fc3 <- mx.symbol.FullyConnected(act2, name="fc3", num_hidden=1)
  softmax <- mx.symbol.SoftmaxOutput(fc3, name="sm")
  nnmodel <- mx.model.FeedForward.create(symbol = softmax, X=data.matrix(trainset[,-c(1,ncol(trainset))]), y=trainset$days,learning.rate=0.005,initializer = mx.init.uniform(0.01))
  preds <- predict(nnmodel, data.matrix(testset[,-c(1,ncol(testset))]))
}

accusvm_RBF_default <- function(trainsvm,testsvm,y,trainindex){
  svm.model <- svm(trainsvm,as.factor(y[trainindex]),probability=TRUE)
  svm.pred <- predict(svm.model, testsvm,probability=TRUE)
  pred <- prediction(attr(svm.pred,"prob")[,1],as.numeric(y[-trainindex]))
  # auc <- colAUC(attr(svm.pred,"prob")[,1],y[-trainindex])
  # auc <- as.numeric(performance(pred,"auc")@y.values)
  roc_obj <- roc(as.numeric(y[-trainindex]), attr(svm.pred,"prob")[,1])
  auc_obj <- as.numeric(roc_obj$auc)
  prec <- unlist(performance(pred,"prec")@y.values)
  rec <- unlist(performance(pred,"rec")@y.values)
  if(is.na(prec[1])) 
    prec[1] <- 1
  height = (prec[-1]+prec[-length(prec)])/2
  width = diff(rec)
  aupr <- sum(height*width)
  c(auc_obj,aupr)
}

normalize_data <- function(x,index) 
{
  #column wise normalization
  if(index==0)
    data.norm <- x
  
  #column wise normalization
  if(index==1)
    data.norm <- x/t(matrix(rep(colSums(x)),nrow = dim(x)[2],ncol = dim(x)[1]))
  
  #row wise norm normalization
  if(index==2)
    data.norm <- x/matrix(rep(rowSums(x)),nrow = dim(x)[1],ncol = dim(x)[2])
  
  #logarithm
  if(index==3){
    if(min(x) <0){
      x <- sapply(x,function(x) x-min(x))
      data.norm <- log10(1+x)
    }else{
      data.norm <- log10(1+x)
    }
  }
  
  #Zscore
  if(index==4)
    data.norm <- scale(x)
  
  #Zscore transposed
  if(index==5)
    data.norm <- t(scale(t(x)))
  
  #quantile normalization
  if(index==6)
    data.norm <- quantileNorm(x)
  
  #quantile normalization transposed
  if(index==7)
    data.norm <- t(quantileNorm(t(x)))
  
  #Sum the abs value of the min number
  if(index==8)
    data.norm <- x + abs(min(x))
  
  #Pareto scaling
  if(index==9) {
    x_center <- x-t(matrix(rep(colMeans(x)),nrow = dim(x)[2],ncol = dim(x)[1]))
    data.norm <- x_center/matrix(rep(t(sqrt(apply(x_center,2,sd))),each=dim(x_center)[1]),ncol = dim(x_center)[2], nrow = dim(x_center)[1])
  }
  
  #square root
  if(index==10){
    if(min(x) < 0){
      x <- sapply(x,function(x) x-min(x))
      data.norm <- sqrt(x)
    }else{
      data.norm <- sqrt(x)
    }
  }
  
  if(index==11){
    data.norm <- log(x+1)
  }
  
  if(index==12){
    for(n in 1:ncol(x)){
      x[,n] <- rank(x[,n])
    }
    data.norm <- x
  }
  data.norm
}

quantileNorm <- function(x) {
  x_rank <- apply(x,2,rank,ties.method="min")
  x_sorted <- data.frame(apply(x,2,sort))
  x_mean <- apply(x_sorted, 1, mean)
  
  indexToMean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  x_final <- apply(x_rank, 2, indexToMean, my_mean=x_mean)
  return(x_final)
}

clinical_svm <- function() {
  matchday <- match(clinicalnona$samplecode,samplesthres$samplecode)
  days <- samplesthres$days[matchday]
  data <- cbind(clinicalnona[,-1],days)[complete.cases(days),]
  alivestatus <- clinicalnona$samplecode %in% deadsamples[,1]
  X <- data[,-ncol(data)]
  y <- data$days >= thresupper
  
  nruns <- 300
  flds1 <- sapply(1:(nruns/10), function(x) createFolds(y, k = 10, list = TRUE, returnTrain = TRUE))
  resultlist <- list()
  
  for(i in 1:nruns){
    trainindex <- flds1[[i]]
    train <- data.frame(X[trainindex,])
    test <- data.frame(X[-trainindex,])
    clinical <- list(train,test)
    svmresult <- svm1(train,y[trainindex],test,y[-trainindex])
    rnfresult <- randomforest1(train,y[trainindex],test,y[-trainindex])
    resultlist[[length(resultlist)+1]] <-  list(svmresult,rnfresult)
  }
  resultlist
}

clinical_svm2 <- function() {
  matchday <- match(clinicalnona$samplecode,samplesthres$samplecode)
  days <- samplesthres$days[matchday]
  data <- cbind(clinicalnona[,-1],days)
  alivestatus <- clinicalnona$samplecode %in% deadsamples[,1]
  X <- data[,-ncol(data)]
  y <- data$days >= thresupper
  
  nruns <- 1000
  flds1 <- sapply(1:100, function(x) createFolds(y, k = 10, list = TRUE, returnTrain = TRUE))
  accu_clinical_all_RBF <- rep(0,nruns)
  aupr_clinical_all_RBF <- rep(0,nruns)
  
  for(i in 1:nruns){
    trainindex <- flds1[[i]]
    train <- data.frame(X[trainindex,])
    test <- data.frame(X[-trainindex,])
    clinical <- list(train,test)
    rbf <- svm1(train,y[trainindex],test,y[-trainindex])
    accu_clinical_all_RBF[i] <- rbf[1]
    aupr_clinical_all_RBF[i] <- rbf[2]
  }
  
  list(accu_clinical_all_RBF,aupr_clinical_all_RBF)
}

getclinicaldata <- function(thresupper, threslower, rmoutlier){
  
  setwd('/home/bzfshaob/Desktop/Experiments/LUAD')
  
  alivesamples <- read.csv(file='Data/clinical/alivesamples.txt',header=TRUE, sep=' ')
  deadsamples <- read.csv(file='Data/clinical/deadsamples.txt',header=TRUE, sep=' ')
  allsamples <- rbind(deadsamples,alivesamples)
  alivesensor <- alivesamples[alivesamples$days >= thresupper,]
  #alivesensor <- alivesamples
  samples <- rbind(deadsamples,alivesensor)
  samplesthres <- samples[samples$days >= thresupper | samples$days <threslower, ]
  
  clinical <- read.csv(file='Data/clinical/clinical_variables_522.csv',header=TRUE,sep=',')
  clinical_se <- clinical[clinical$samplecode %in% samplesthres$samplecode, ]
  clivars <- clinical_se[,c(1,3,12,10,17,18,19,20,21)]
  clivars$pathologic_t <- replace(clivars$pathologic_t,clivars$pathologic_t=="tx",NA)
  clivars$pathologic_n <- replace(clivars$pathologic_n,clivars$pathologic_n=="nx",NA)
  clinicalnona <- clivars[complete.cases(clivars),]
  # assign levels to each clinical variable
  #pathologic_t: change levels to t1-t4
  clinicalnona$pathologic_t <- replace(clinicalnona$pathologic_t, clinicalnona$pathologic_t=="t1a", "t1")
  clinicalnona$pathologic_t <-replace(clinicalnona$pathologic_t, clinicalnona$pathologic_t=='t1b', 't1')
  clinicalnona$pathologic_t <- replace(clinicalnona$pathologic_t, clinicalnona$pathologic_t=='t2a', 't2')
  clinicalnona$pathologic_t <-replace(clinicalnona$pathologic_t, clinicalnona$pathologic_t=='t2b', 't2')
  clinicalnona$pathologic_t <- factor(clinicalnona$pathologic_t, levels= c("t1","t2","t3","t4"))
  #pathologic_n: change levels to n0 to n3
  clinicalnona$pathologic_n <- factor(clinicalnona$pathologic_n, levels= c("n0","n1","n2","n3"))
  #gender: split to two variables
  #clinicalnona$female <- clinicalnona$gender=='female'
  clinicalnona$male <- clinicalnona$gender=='male'
  clinicalnona$gender <- NULL
  #pathologic stage: change to stage 1-4
  pastage <- as.character(clinicalnona$pathologic_stage)
  pastage <- replace(pastage,pastage=="stage ia","i")
  pastage <- replace(pastage,pastage=="stage ib","i")
  pastage <- replace(pastage,pastage=="stage iia","ii")
  pastage <- replace(pastage,pastage=="stage iib","ii")
  pastage <- replace(pastage,pastage=="stage iiia","iii")
  pastage <- replace(pastage,pastage=="stage iiib","iii")
  pastage <- replace(pastage,pastage=="stage iv","iv")
  pastage <- replace(pastage,pastage=="stage i","i")
  pastage <- replace(pastage,pastage=="stage ii","ii")
  clinicalnona$pathologic_stage <- factor(pastage, levels = c('i','ii','iii','iv'))
  #tobacco smoking history: integer 1 to 5
  clinicalnona$tobacco_smoking_history <- factor(clinicalnona$tobacco_smoking_history,levels = c('1','2','3','4','5'))
  #clinicalnona$age_at_initial_pathologic_diagnosis: integer
  #radiation therapy: split to two variables
  clinicalnona$radiationthyyes <- clinicalnona$radiation_therapy=='yes'
  #clinicalnona$radiationthyno <- clinicalnona$radiation_therapy=='no'
  clinicalnona$radiation_therapy <- NULL
  #target molecular therapy: split to two variables
  clinicalnona$molthyyes <- clinicalnona$targeted_molecular_therapy=='yes'
  #clinicalnona$molthyno <- clinicalnona$targeted_molecular_therapy=='no'
  clinicalnona$targeted_molecular_therapy <- NULL
  #age
  clinicalnona$age_at_initial_pathologic_diagnosis <- 
    cut(clinicalnona$age_at_initial_pathologic_diagnosis,c(30,50,70,90),labels = c("young","middle","old"))
  #change data to numeric
  clinicalnona$pathologic_t <- as.numeric(clinicalnona$pathologic_t)
  clinicalnona$pathologic_n <- as.numeric(clinicalnona$pathologic_n)
  clinicalnona$pathologic_stage <- as.numeric(clinicalnona$pathologic_stage)
  clinicalnona$tobacco_smoking_history <- as.numeric(clinicalnona$tobacco_smoking_history)
  clinicalnona$age_at_initial_pathologic_diagnosis <- as.numeric(clinicalnona$age_at_initial_pathologic_diagnosis)
  clinicalnona[,-1] <- clinicalnona[,-1]+0
  
  commsample <- intersect(clinicalnona$samplecode,rmoutlier[,1])
  days <- samplesthres[match(commsample,samplesthres$samplecode),2] > thresupper
  clinicaldata <- clinicalnona[match(commsample,clinicalnona[,1]),-1]
  moldata <- rmoutlier[match(commsample,rmoutlier[,1]),-c(1,ncol(rmoutlier))]
  setwd('/home/bzfshaob/Desktop/Experiments_largeEMT_correct/')
  
  list(clinicaldata,moldata,days)
}
