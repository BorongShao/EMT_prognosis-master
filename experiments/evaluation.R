# select features based on top ranking ones
domtsinglefeatures <- function(genesidx, trainset,times1){
  y <- trainset$days
  X <- trainset[,genesidx+1]
  newx <- X
  remove <- rep(FALSE,ncol(newx))
  for(j in 1:(ncol(newx)-1)){
    if(remove[j]==FALSE){
      for(i in (j+1):ncol(newx)){
        if(abs(cor(newx[,i],newx[,j],method="spearman")) >= times1*abs(cor(newx[,i],y,method="spearman")))
          remove[i] <- TRUE
      }
    }
  }
  which(remove==FALSE)
}

domtnetfeatures <- function(netlist, trainset, time1){
  
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset$days
  directionadd <- sign(cor(X,y))
  newX <- sapply(1:ncol(X), function(x) X[,x]*directionadd[x])
  traindata <- sapply(netlist,function(x) rowSums(newX[,x])/length(x))
  newx <- traindata
  remove <- rep(FALSE,ncol(newx))
  for(j in 1:(ncol(newx)-1)){
    if(remove[j]==FALSE){
      for(i in (j+1):ncol(newx)){
        if(abs(cor(newx[,i],newx[,j], method = "spearman")) >= time1*abs(cor(newx[,i],y, method = "spearman")))
          remove[i] <- TRUE
      }
    }
  }
  netlist[which(remove==FALSE)]
}

domtnetfeaturesCV <- function(netlist, trainset){
  X <- sapply(netlist,function(x) rowSums(trainset[,(x+1)])/length(x))
  y <- trainset$days
  num_scores <- matrix(0,10,30)
  for (i in 1:10){
    trainidx <- sample(nrow(X),floor(nrow(X)/5))
    traininner <- X[trainidx,]
    testinner <- X[-trainidx,]
    trainy <- y[trainidx]
    testy <- y[-trainidx]
    num_scores[i,] <- sapply(1:30, function(x) svm1(traininner[,1:x],trainy,testinner[,1:x],testy)[1])
  }
  1:which.max(colSums(num_scores))
}

cal_aupr <- function(xx.df) {
  perf  <- performance(xx.df, "prec", "rec")
  xy    <- data.frame(recall=perf@x.values[[1]], precision=perf@y.values[[1]])
  
  # take out division by 0 for lowest threshold
  xy <- subset(xy, !is.nan(xy$precision))
  
  # Designate recall = 0 as precision = x...arbitrary
  xy <- rbind(c(0, 0), xy)
  #xy <- xy[!(rowSums(xy)==0), ]
  
  res   <- trapz(xy$recall, xy$precision)
  res
}

#-------------------------------------binary class without survival models-------------------------#

# single genes, binary class
# input gene index vector, training and testing data without samplecode column
# genesidx_ttest,trainset,testset,2,testset2
sgclass_old <- function(geneindex,trainset,testset,md){
  geneindex <- geneindex + 1
  traindata <- data.matrix(trainset[,geneindex])
  trainy <- trainset$days
  testdata <- data.matrix(testset[,geneindex])
  testy <- testset$days
  if(md==1)
    logisticreg1(as.data.frame(traindata),trainy,as.data.frame(testdata),testy)
  else if(md==2)
    svm1(traindata,trainy,testdata,testy)
  else if(md==3)
    randomforest1(traindata,trainy,testdata,testy)
  else if(md==4)
    deeplearning1(traindata,trainy,testdata,testy)
}

sgclass <- function(geneindex,trainset,testset,md, testset2){
  geneindex <- geneindex + 1
  traindata <- data.matrix(trainset[,geneindex])
  trainy <- trainset$days
  testdata <- data.matrix(testset[,geneindex])
  testy <- testset$days
  temp <- testset2[,geneindex]
  testset2 <- data.frame(testset2[,1],temp,testset2$days)
  names(testset2)[ncol(testset2)] <- "days"
  if(md==1)
    logisticreg1(as.data.frame(traindata),trainy,as.data.frame(testdata),testy)
  else if(md==2)
    svm_cox(traindata,trainy,testdata,testy,testset2)
  else if(md==3)
    randomforest_cox(traindata,trainy,testdata,testy,testset2)
  else if(md==4)
    deeplearning1(traindata,trainy,testdata,testy)
}
sgclass_clinical <- function(geneindex,trainset,testset,md){
  geneindex <- geneindex + 1
  traindata <- data.matrix(trainset[,geneindex])
  trainy <- trainset$days
  testdata <- data.matrix(testset[,geneindex])
  testy <- testset$days
  if(md==1)
    logisticreg1(as.data.frame(traindata),trainy,as.data.frame(testdata),testy)
  else if(md==2)
    svm1(traindata,trainy,testdata,testy)
  else if(md==3)
    randomforest1(traindata,trainy,testdata,testy)
  else if(md==4)
    deeplearning1(traindata,trainy,testdata,testy)
}

# subnetworks, binary class, aggregate  networks
# input subnetwork list of gene indecies, training and testing data without samplecode column
netclassaggre <- function(netlist,trainset,testset,md,testset2){
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset$days
  directionadd <- sign(cor(X,y))
  newX <- sapply(1:ncol(X), function(x) X[,x]*directionadd[x])
  traindata <- sapply(netlist,function(x) rowSums(newX[,x])/length(x))
  trainy <- trainset$days
  
  X2 <- testset[,-c(1,ncol(testset))]
  newXtest <- sapply(1:ncol(X2), function(x) X2[,x]*directionadd[x])
  testdata <- sapply(netlist,function(x) rowSums(newXtest[,x])/length(x))
  testy <- testset$days
  
  X3 <- testset2[,-c(1,ncol(testset2))]
  newXtest2 <- sapply(1:ncol(X3), function(x) X3[,x]*directionadd[x])
  testdata2 <- sapply(netlist,function(x) rowSums(newXtest2[,x])/length(x))
  testset2 <- data.frame(testset2[,1],testdata2, testset2$days)
  names(testset2)[ncol(testset2)] <- "days"
  if(md==1)
    logisticreg1(as.data.frame(traindata),trainy,as.data.frame(testdata),testy)
  else if(md==2)
    svm_cox(traindata,trainy,testdata,testy,testset2)
  else if(md==3)
    randomforest_cox(traindata,trainy,testdata,testy,testset2)
  else if(md==4)
    deeplearning1(traindata,trainy,testdata,testy)
}

netclassaggre_old <- function(netlist,trainset,testset,md){
  X <- trainset[,-c(1,ncol(trainset))]
  y <- trainset$days
  directionadd <- sign(cor(X,y))
  newX <- sapply(1:ncol(X), function(x) X[,x]*directionadd[x])
  traindata <- sapply(netlist,function(x) rowSums(newX[,x])/length(x))
  trainy <- trainset$days
  
  X2 <- testset[,-c(1,ncol(testset))]
  newXtest <- sapply(1:ncol(X2), function(x) X2[,x]*directionadd[x])
  testdata <- sapply(netlist,function(x) rowSums(newXtest[,x])/length(x))
  testy <- testset$days
  
  if(md==1)
    logisticreg1(as.data.frame(traindata),trainy,as.data.frame(testdata),testy)
  else if(md==2)
    svm1(traindata,trainy,testdata,testy)
  else if(md==3)
    randomforest1(traindata,trainy,testdata,testy)
  else if(md==4)
    deeplearning1(traindata,trainy,testdata,testy)
}


netclassaggre_tune <- function(netlist,trainset,testset){
  traindata <- sapply(netlist,function(x) rowSums(trainset[,(x+1)])/length(x))
  trainy <- trainset$days
  testdata <- sapply(netlist,function(x) rowSums(testset[,(x+1)])/length(x))
  testy <- testset$days
  
  tar <- as.factor(trainy)
  gammavec <- rep(0,10)
  costvec <- rep(0,10)
  accuvec <- rep(0,10)
  count = 0
  flds <- createFolds(trainy, k=10)
  for(j in 1:10){
    count <- count + 1
    tuneResult <- tune(svm, train.x = traindata[-flds[[j]],], train.y=tar[-flds[[j]]],
                       validation.x = traindata[flds[[j]],], validation.y = tar[flds[[j]]],
                       ranges = list(gamma = seq(0,0.2,0.01), cost = 2^(-5:5)))
    gammavec[count] <- as.numeric(tuneResult$best.parameters[1])
    costvec[count] <- as.numeric(tuneResult$best.parameters[2])
    accuvec[count] <- tuneResult$best.performance
  }
  index <- which.max(accuvec)
  c <- costvec[index]
  g <- gammavec[index]
  svm.model <- svm(traindata,tar,gamma = g,cost = c, probability=TRUE)
  svm.pred <- predict(svm.model, testdata,probability=TRUE)
  #pred <- prediction(attr(svm.pred,"prob")[,1],as.numeric(y[-trainindex]))
  auc <- colAUC(attr(svm.pred,"prob")[,1],testy)
  #auc <- as.numeric(performance(pred,"auc")@y.values)
  auc
}


# subnetworks, binary class, pcs of networks
# input subnetwork list of gene indecies, training and testing data without samplecode column
netclasspc <- function(netlist,trainset,testset, md){
  PCtrain <- NULL
  PCtest <- NULL
  for(i in 1:length(netlist)){
    trainnetdata <- trainset[,(netlist[[i]]+1)]
    testnetdata <- testset[,(netlist[[i]]+1)]
    train.pc <- prcomp(trainnetdata,center=TRUE,scale. = TRUE)
    PCtrain = cbind(PCtrain, predict(train.pc,trainnetdata)[,1:3])
    PCtest = cbind(PCtest, predict(train.pc,testnetdata)[,1:3])
  }
  traindata <- PCtrain
  trainy <- trainset$days
  testdata <- PCtest
  testy <- testset$days
  if(md==1)
    logisticreg1(as.data.frame(traindata),trainy,as.data.frame(testdata),testy)
  else if(md==2)
    svm1(traindata,trainy,testdata,testy)
  else if(md==3)
    randomforest1(traindata,trainy,testdata,testy)
  else if(md==4)
    deeplearning1(traindata,trainy,testdata,testy)
}

# single genes, binary class, vector 
sgclassvec <- function(beta,testset) {
  testdata <- testset[,-c(1,ncol(testset))]
  testy <- testset$days
  pred <- data.matrix(testdata) %*% beta[-1] + beta[1]
  pred2 <- prediction(pred,testy)
  roc_obj <- roc(as.numeric(testy),as.numeric(pred))
  auc <- as.numeric(roc_obj$auc)
  prec <- unlist(performance(pred2,"prec")@y.values)
  rec <- unlist(performance(pred2,"rec")@y.values)
  if(is.na(prec[1])) 
    prec[1] <- 1
  #plot(rec,prec)
  height = (prec[-1]+prec[-length(prec)])/2
  width = diff(rec)
  aupr <- sum(height*width)
  c(auc,aupr)
}


#-------------------------------------survival object ------------------------------------------#
#intput gene index vector, training and testing data with samplecode column, remember to add 1 to gene index
sgsurvival <- function(geneindex,trainset,testset){
  traindata <- trainset[,geneindex+1]
  trainy <- trainset$days
  testdata <- testset[,geneindex+1]
  testy <- testset$days
  trainsample <- trainset[,1]
  testsample <- testset[,1]
  survobjtrain <- Surv(trainy, trainsample %in% deadsamples[,1])
  survobjtest <- Surv(testy, testsample %in% deadsamples[,1])
  fit <- coxph(survobjtrain ~ ., traindata,iter.max = 20000)
  # fit <- coxph(survobjtrain ~ traindata[,1],iter.max = 20000)
  # prediction
  predvec <- predict(fit, newdata=testdata, type='risk')
  if(length(unique(predvec))==1)
    preddiscrete <- sample(2,length(predvec),replace=TRUE) 
  else
    preddiscrete <- as.numeric(cut(predvec,quantile(predvec,  probs = c(0,0.5,1))))
  
  # log-rank test
  sdiff <- survdiff(survobjtest ~ preddiscrete)$chisq
  pval <- 1- pchisq(sdiff,1) 
  
  #Kaplan-Meier curve
  isdead <- testsample %in% deadsamples[,1]
  fit1 <- survfit(Surv(testy[preddiscrete==1], isdead[preddiscrete==1]) ~ 1)
  fit2 <- survfit(Surv(testy[preddiscrete==2], isdead[preddiscrete==2]) ~ 1)
  # plot(fit1, col = 'blue', xlab = 'Time (days)', ylab = 'Survival Probability')
  # lines(fit2, col = 'red')
  # legend(21,1,c('Group 1 (treatment)', 'Group 2 (placebo)'), col = c('blue','red'), lty = 1)
  # title(main='KM-Curves for Remission Data')
  
  # time dependent AUC
  # roccurve <- survivalROC(testy, isdead, predvec, predict.time=365, cut.values =
  #                           NULL, method = "KM", lambda = NULL, span = NULL, window =
  #                           "symmetric")
  # plot(roccurve$FP, roccurve$TP, type="l", xlim=c(0,1), ylim=c(0,1),
  #      xlab=paste( "FP", "\n", "AUC = ",round(roccurve$AUC,3)),
  #      ylab="TP",main="Mayoscore 4, Method = KM \n Year = 1")
  # abline(0,1)
  # auc <- roccurve$AUC
  
  #AUC valus for many time poinst
  # ptm <- proc.time()
  auclist <- sapply(1:60, function(x) survivalROC(testy, isdead, predvec, predict.time=x*30, cut.values =
                                                    NULL, method = "KM", lambda = NULL, span = NULL, window =
                                                    "symmetric")$AUC)
  # proc.time() - ptm
  # 
  # registerDoMC(cores=2)
  # results <- foreach(i = 1:60) %dopar% {
  #   survivalROC(testy, isdead, predvec, predict.time=i, cut.values = NULL, method = "KM", 
  #               lambda = NULL, span = NULL, window ="symmetric")$AUC
  # }
  # auclist <- unlist(results)
  list(pval,auclist)
  # c(pval,auc)
}


#intput subnetwork list of gene indecies, training and testing data with samplecode column
netsurvivalaggre <- function(netlist,trainset,testset){
  traindata <- sapply(netlist,function(x) rowSums(trainset[,x+1])/length(x))
  trainy <- trainset$days
  testdata <- sapply(netlist,function(x) rowSums(testset[,x+1])/length(x))
  testy <- testset$days
  trainsample <- trainset[,1]
  testsample <- testset[,1]
  survobjtrain <- Surv(trainy, trainsample %in% deadsamples[,1])
  survobjtest <- Surv(testy, testsample %in% deadsamples[,1])
  fit <- coxph(survobjtrain~., data=as.data.frame(traindata),iter.max = 20000)
  # prediction
  predvec <- predict(fit, as.data.frame(testdata), type='risk')
  preddiscrete <- as.numeric(cut(predvec,quantile(predvec,  probs = c(0,0.5,1))))
  
  # log-rank test
  sdiff <- survdiff(survobjtest ~ preddiscrete)$chisq
  pval <- 1- pchisq(sdiff,1) 
  
  #Kaplan-Meier curve
  isdead <- testsample %in% deadsamples[,1]
  fit1 <- survfit(Surv(testy[preddiscrete==1], isdead[preddiscrete==1]) ~ 1)
  fit2 <- survfit(Surv(testy[preddiscrete==2], isdead[preddiscrete==2]) ~ 1)
  plot(fit1, col = 'blue', xlab = 'Time (days)', ylab = 'Survival Probability')
  lines(fit2, col = 'red')
  legend(21,1,c('Group 1 (treatment)', 'Group 2 (placebo)'), col = c('blue','red'), lty = 1)
  title(main='KM-Curves for Remission Data')
  
  # time dependent AUC
  # roccurve <- survivalROC(testy, isdead, predvec, predict.time=365, cut.values =
  #                           NULL, method = "KM", lambda = NULL, span = NULL, window =
  #                           "symmetric")
  # # plot(roccurve$FP, roccurve$TP, type="l", xlim=c(0,1), ylim=c(0,1),
  # #      xlab=paste( "FP", "\n", "AUC = ",round(roccurve$AUC,3)),
  # #      ylab="TP",main="Mayoscore 4, Method = KM \n Year = 1")
  # # abline(0,1)
  # auc <- roccurve$AUC
  
  #AUC valus for many time poinst
  # ptm <- proc.time()
  auclist <- sapply(1:60, function(x) survivalROC(testy, isdead, predvec, predict.time=x*30, cut.values =
                                                    NULL, method = "KM", lambda = NULL, span = NULL, window =
                                                    "symmetric")$AUC)
  list(pval,auclist)
  # c(pval,auc)
}


# subnetworks, pcs of networks
# input subnetwork list of gene indecies, training and testing data with samplecode column
netsurvivalpc <- function(netlist,trainset,testset){
  PCtrain <- NULL
  PCtest <- NULL
  for(i in 1:length(netlist)){
    trainnetdata <- trainset[,netlist[[i]]+1]
    testnetdata <- testset[,netlist[[i]]+1]
    train.pc <- prcomp(trainnetdata,center=TRUE,scale. = TRUE)
    PCtrain = cbind(PCtrain, predict(train.pc,trainnetdata)[,1:3])
    PCtest = cbind(PCtest, predict(train.pc,testnetdata)[,1:3])
  }
  traindata <- PCtrain
  trainy <- trainset$days
  testdata <- PCtest
  testy <- testset$days
  trainsample <- trainset[,1]
  testsample <- testset[,1]
  survobjtrain <- Surv(trainy, trainsample %in% deadsamples[,1])
  survobjtest <- Surv(testy, testsample %in% deadsamples[,1])
  fit <- coxph(survobjtrain~., data=data.frame(traindata),iter.max = 20000)
  # prediction
  predvec <- predict(fit, data.frame(testdata), type='risk')
  preddiscrete <- as.numeric(cut(predvec,quantile(predvec,  probs = c(0,0.5,1))))
  
  # log-rank test
  sdiff <- survdiff(survobjtest ~ preddiscrete)$chisq
  pval <- 1- pchisq(sdiff,1) 
  
  #Kaplan-Meier curve
  isdead <- testsample %in% deadsamples[,1]
  fit1 <- survfit(Surv(testy[preddiscrete==1], isdead[preddiscrete==1]) ~ 1)
  fit2 <- survfit(Surv(testy[preddiscrete==2], isdead[preddiscrete==2]) ~ 1)
  plot(fit1, col = 'blue', xlab = 'Time (days)', ylab = 'Survival Probability')
  lines(fit2, col = 'red')
  legend(21,1,c('Group 1 (treatment)', 'Group 2 (placebo)'), col = c('blue','red'), lty = 1)
  title(main='KM-Curves for Remission Data')
  
  # time dependent AUC
  # roccurve <- survivalROC(testy, isdead, predvec, predict.time=365, cut.values =
  #                           NULL, method = "KM", lambda = NULL, span = NULL, window =
  #                           "symmetric")
  # # plot(roccurve$FP, roccurve$TP, type="l", xlim=c(0,1), ylim=c(0,1),
  # #      xlab=paste( "FP", "\n", "AUC = ",round(roccurve$AUC,3)),
  # #      ylab="TP",main="Mayoscore 4, Method = KM \n Year = 1")
  # # abline(0,1)
  # auc <- roccurve$AUC
  # 
  #AUC valus for many time poinst
  # ptm <- proc.time()
  auclist <- sapply(1:60, function(x) survivalROC(testy, isdead, predvec, predict.time=x*30, cut.values =
                                                    NULL, method = "KM", lambda = NULL, span = NULL, window =
                                                    "symmetric")$AUC)
  # proc.time() - ptm
  # 
  # registerDoMC(cores=2)
  # results <- foreach(i = 1:60) %dopar% {
  #   survivalROC(testy, isdead, predvec, predict.time=i, cut.values = NULL, method = "KM", 
  #               lambda = NULL, span = NULL, window ="symmetric")$AUC
  # }
  # auclist <- unlist(results)
  list(pval,auclist)
  # c(pval,auc)
}


# subnetworks, vector, withsamplecode column
sgsurvivalvec <- function(testset,predvec) {
  testdata <- testset[,-c(1,ncol(testset))]
  testy <- testset$days
  testsample <- testset[,1]
  survobjtest <- Surv(testy, testsample %in% deadsamples[,1])
  if(length(unique(predvec))==1)
    preddiscrete <- sample(2,length(predvec),replace=TRUE) else
    preddiscrete <- as.numeric(cut(predvec,quantile(predvec,  probs = c(0,0.5,1))))
  
  # log-rank test
  sdiff <- survdiff(survobjtest ~ preddiscrete)$chisq
  pval <- 1- pchisq(sdiff,1) 
  
  #Kaplan-Meier curve
  isdead <- testsample %in% deadsamples[,1]
  # fit1 <- survfit(Surv(testy[preddiscrete==1], isdead[preddiscrete==1]) ~ 1)
  # fit2 <- survfit(Surv(testy[preddiscrete==2], isdead[preddiscrete==2]) ~ 1)
  # plot(fit1, col = 'blue', xlab = 'Time (days)', ylab = 'Survival Probability')
  # lines(fit2, col = 'red')
  # legend(2200,1,c('Good prognosis group', 'Poor prognosis group'), col = c('blue','red'), lty = 1)
  # title(main='KM-Curves for Remission Data')
  
  # time dependent AUC
  # roccurve <- survivalROC(testy, isdead, predvec, predict.time=365, cut.values =
  #                           NULL, method = "KM", lambda = NULL, span = NULL, window =
  #                           "symmetric")
  # # plot(roccurve$FP, roccurve$TP, type="l", xlim=c(0,1), ylim=c(0,1),
  # #      xlab=paste( "FP", "\n", "AUC = ",round(roccurve$AUC,3)),
  # #      ylab="TP",main="Mayoscore 4, Method = KM \n Year = 1")
  # # abline(0,1)
  # auc <- roccurve$AUC
  
  #AUC valus for many time poinst
  # ptm <- proc.time()
  auclist <- sapply(1:60, function(x) survivalROC(testy, isdead, predvec, predict.time=x*30, cut.values =
                                                    NULL, method = "KM", lambda = NULL, span = NULL, window =
                                                    "symmetric")$AUC)
  # proc.time() - ptm
  # 
  # registerDoMC(cores=2)
  # results <- foreach(i = 1:60) %dopar% {
  #   survivalROC(testy, isdead, predvec, predict.time=i, cut.values = NULL, method = "KM", 
  #               lambda = NULL, span = NULL, window ="symmetric")$AUC
  # }
  # auclist <- unlist(results)
  list(pval,auclist)
  # c(pval,auc)
}


#for MCL_superpc where the output is not genelist but pcs

survivalpc <- function(pc,testset) {
  testdata <- testset[,-c(1,ncol(testset))]
  testy <- testset$days
  testsample <- testset[,1]
  survobjtest <- Surv(testy, testsample %in% deadsamples[,1])
  pcpvals <- rep(0,ncol(pc))
  for( i in 1:length(pcpvals)){
    preddiscrete <- as.numeric(cut(pc[,i],quantile(pc[,i],  probs = c(0,0.5,1))))
    sdiff <- survdiff(survobjtest ~ preddiscrete)$chisq
    pcpvals[i] <- 1- pchisq(sdiff,1) 
  }
  idx <- which.min(pcpvals)
  predvec <- pc[,idx]
  preddiscrete <- as.numeric(cut(predvec,quantile(predvec,  probs = c(0,0.5,1))))
  
  # log-rank test
  sdiff <- survdiff(survobjtest ~ preddiscrete)$chisq
  pval <- 1- pchisq(sdiff,1) 
  
  #Kaplan-Meier curve
  isdead <- testsample %in% deadsamples[,1]
  fit1 <- survfit(Surv(testy[preddiscrete==1], isdead[preddiscrete==1]) ~ 1)
  fit2 <- survfit(Surv(testy[preddiscrete==2], isdead[preddiscrete==2]) ~ 1)
  plot(fit1, col = 'blue', xlab = 'Time (days)', ylab = 'Survival Probability')
  lines(fit2, col = 'red')
  legend(21,1,c('Group 1 (treatment)', 'Group 2 (placebo)'), col = c('blue','red'), lty = 1)
  title(main='KM-Curves for Remission Data')
  
  # time dependent AUC
  # roccurve <- survivalROC(testy, isdead, predvec, predict.time=365, cut.values =
  #                           NULL, method = "KM", lambda = NULL, span = NULL, window =
  #                           "symmetric")
  # # plot(roccurve$FP, roccurve$TP, type="l", xlim=c(0,1), ylim=c(0,1),
  # #      xlab=paste( "FP", "\n", "AUC = ",round(roccurve$AUC,3)),
  # #      ylab="TP",main="Mayoscore 4, Method = KM \n Year = 1")
  # # abline(0,1)
  # auc <- roccurve$AUC
  
  #AUC valus for many time poinst
  # ptm <- proc.time()
  auclist <- sapply(1:60, function(x) survivalROC(testy, isdead, predvec, predict.time=x*30, cut.values =
                                                    NULL, method = "KM", lambda = NULL, span = NULL, window =
                                                    "symmetric")$AUC)
  # proc.time() - ptm
  # 
  # registerDoMC(cores=2)
  # results <- foreach(i = 1:60) %dopar% {
  #   survivalROC(testy, isdead, predvec, predict.time=i, cut.values = NULL, method = "KM", 
  #               lambda = NULL, span = NULL, window ="symmetric")$AUC
  # }
  # auclist <- unlist(results)
  list(pval,auclist)
  # c(pval,auc)
}


