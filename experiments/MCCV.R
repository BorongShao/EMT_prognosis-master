# Monto-Carlo cross validation 
#source('data_network.R')
source('implementations/Hubc/hub_network.R')
source('implementations/stSVM/stSVM.R')
source('evaluation.R')
source('implementations/NoNetwork/no_network.R')
source('implementations/Chuang2007/additive_improve.R')
source('implementations/NetRank/NetRank.R')
# source('implementations/MCL_superpc/superpc.R')
source('implementations/Survnet_Ideker/survnet_ideker.R')
source('implementations/RDS/RDS_RGS.R')
# source('implementations/WGCNA/WGCNA.R')
source('implementations/NetLasso/NetworkLasso.R')

repetition <- 100
ntopgene <- 50
ntopnet <- 50

# sampling <- function(data, percent){
#   np <- sum(data$days>0)
#   nn <- sum(data$days==0)
#   posidx <- sample(which(data$days >0),floor(np*percent))
#   negidx <- sample(which(data$days ==0),floor(nn*percent))
#   sort(c(posidx,negidx))
# }
single_level_methods1 <- function(thresdata, threstrainidx){
  
  #use all features
  # trainset <- allrnaseqfeatures[threstrainidx,]
  # testset <- allrnaseqfeatures[-threstrainidx,]
  # rmfeatures <- which(sapply(trainset[,-c(1,ncol(trainset))], var)<0.1)+1
  # if(length(rmfeatures)!=0){
  #   trainset <- allrnaseqfeatures[threstrainidx,-rmfeatures]
  #   testset <- allrnaseqfeatures[-threstrainidx,-rmfeatures]
  # }else{
  #   trainset <- allrnaseqfeatures[threstrainidx,]
  #   testset <- allrnaseqfeatures[-threstrainidx,]
  # }
  # svm_allfeatures <- svm1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  # lgs_allfeatures <- logisticreg1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  # rnf_allfeatures <- randomforest1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  #use random features
  # trainset <- LEthresgeneexp_rand[threstrainidx,]
  # testset <- LEthresgeneexp_rand[-threstrainidx,]
  # svm_randfeatures <- svm1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  # EMT
  trainset <- thresdata[threstrainidx,]
  testset <- thresdata[-threstrainidx,]
  # use all EMT features
  svm_allemtfeatures <- svm1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  rnf_allemtfeatures <- randomforest1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  # dn_allemtfeatures <- deeplearning1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  
  # rank genes by ttest
  genesidx_ttest <- ttestgenes(trainset,ntopgene)
  # genesidx <- domtsinglefeatures(genesidx2,trainset,5)
  svm_ttest <- sgclass_old(genesidx_ttest,trainset,testset,2)
  rnf_ttest <- sgclass_old(genesidx_ttest,trainset,testset,3)
  # dn_ttest <- sgclass(genesidx,trainset,testset,4)
  ##ridge regression
  # # beta <- logistic_reg(trainset)
  # # result_logicreg <- sgclassvec(beta,testset)
  #lasso
  genesidx_lasso <- lassocoefs(trainset)
  # genesidx <- domtsinglefeatures(genesidx2,trainset,5)
  svm_lasso <- sgclass_old(genesidx_lasso,trainset,testset,2)
  rnf_lasso <- sgclass_old(genesidx_lasso,trainset,testset,3)
  # dn_lasso <- sgclass(genesidx,trainset,testset,4)
  # ---------------------------with network
  #netowrk lasso
  genesidx_netlasso <- netLasso(trainset)
  # genesidx <- domtsinglefeatures(genesidx2,trainset,1)
  svm_netlasso <- sgclass_old(genesidx_netlasso,trainset,testset,2)
  rnf_netlasso <- sgclass_old(genesidx_netlasso,trainset,testset,3)
  # dn_netlasso <- sgclass(genesidx,trainset,testset,4)
  #additive subnetworks
  addinets <- additive_subnet(trainset,ntopnet)
  # addinets <- domtnetfeatures(addinets2,trainset,10)
  svm_aggreadd <-netclassaggre_old(addinets,trainset,testset,2)
  rnf_aggreadd <-netclassaggre_old(addinets,trainset,testset,3)
  # dn_aggreadd <-netclassaggre(addinets,trainset,testset,4)
  # svm_pcadd <-netclasspc(addinets,trainset,testset,2)
  # lgs_pcadd <-netclasspc(addinets,trainset,testset,1)
  # rnf_pcadd <-netclasspc(addinets,trainset,testset,3)
  #NetRank
  genesidx_netrank <- NetRank(trainset,ntopgene)
  # genesidx <- domtsinglefeatures(genesidx2,trainset,5)
  svm_netrank <- sgclass_old(genesidx_netrank,trainset,testset,2)
  rnf_netrank <- sgclass_old(genesidx_netrank,trainset,testset,3)
  # dn_netrank <- sgclass(genesidx,trainset,testset,4)
  
  list(svm_allemtfeatures,rnf_allemtfeatures,svm_ttest,rnf_ttest,genesidx_ttest,svm_lasso,rnf_lasso,genesidx_lasso,svm_netlasso,rnf_netlasso,genesidx_netlasso,
       svm_aggreadd,rnf_aggreadd,addinets,svm_netrank,rnf_netrank,genesidx_netrank,threstrainidx)
}

evaluation_algorithms <- function(thresdata, threstrainidx, nothresdata, nothrestrainidx){

  # EMT
  trainset <- thresdata[threstrainidx,]
  testset <- thresdata[-threstrainidx,]
  trainset2 <- nothresdata[nothrestrainidx,]
  testset2 <- nothresdata[-nothrestrainidx,]
  names(trainset2) <- gsub("-",'.', names(trainset2),fixed=TRUE)
  names(testset2) <- gsub("-",'.', names(testset2),fixed=TRUE)
  
  # use all EMT features
  svm_allemtfeatures <- svm_cox(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days,testset2)
  rnf_allemtfeatures <- randomforest_cox(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days,testset2)
  
  #wilcoxon t test
  genesidx_ttest <- ttestgenes(trainset,ntopgene)
  svm_ttest <- sgclass(genesidx_ttest,trainset,testset,2,testset2)
  rnf_ttest <- sgclass(genesidx_ttest,trainset,testset,3,testset2)
  
  #lasso
  genesidx_lasso <- lassocoefs(trainset)
  svm_lasso <- sgclass(genesidx_lasso,trainset,testset,2,testset2)
  rnf_lasso <- sgclass(genesidx_lasso,trainset,testset,3,testset2)
  # ---------------------------with network
  #netowrk lasso
  genesidx_netlasso <- netLasso(trainset)
  svm_netlasso <- sgclass(genesidx_netlasso,trainset,testset,2,testset2)
  rnf_netlasso <- sgclass(genesidx_netlasso,trainset,testset,3,testset2)
  
  #network smoothed t-statistic
  genesidx_stSVM <- stSVM(trainset,ntopgene)
  svm_stSVM <- sgclass(genesidx_stSVM,trainset,testset,2,testset2)
  rnf_stSVM <- sgclass(genesidx_stSVM,trainset,testset,3,testset2)

  #additive subnetworks
  addinets <- additive_subnet(trainset,ntopnet)
  svm_aggreadd <-netclassaggre(addinets,trainset,testset,2,testset2)
  rnf_aggreadd <-netclassaggre(addinets,trainset,testset,3,testset2)

  #NetRank
  genesidx_netrank <- NetRank(trainset,ntopgene)
  svm_netrank <- sgclass(genesidx_netrank,trainset,testset,2,testset2)
  rnf_netrank <- sgclass(genesidx_netrank,trainset,testset,3,testset2)
 
  #cox 
  genesidx_cox <- Coxrank(trainset2,ntopgene)
  svm_coxfil <- sgclass(genesidx_cox,trainset,testset,2,testset2)
  rnf_coxfil <- sgclass(genesidx_cox,trainset,testset,3,testset2)

  #cox_lasso
  genesidx_coxlasso <- which(lassocoxcoefs(trainset2) != 0)
  svm_coxlasso <- sgclass(genesidx_coxlasso,trainset,testset,2,testset2)
  rnf_coxlasso <- sgclass(genesidx_coxlasso,trainset,testset,3,testset2)
  
  #RDS
  genesidx_RDS <- RDS_RGS(trainset2,ntopgene)
  svm_rds <- sgclass(genesidx_RDS,trainset,testset,2,testset2)
  rnf_rds <- sgclass(genesidx_RDS,trainset,testset,3,testset2)
  
  #Survnet
  survnets <-survnet(trainset2,ntopnet)[[1]]
  svm_survnet <-netclassaggre(survnets,trainset,testset,2,testset2)
  rnf_survnet <-netclassaggre(survnets,trainset,testset,3,testset2)
  
  list(svm_allemtfeatures,genesidx_ttest,svm_ttest,genesidx_lasso,svm_lasso,genesidx_netlasso,svm_netlasso,genesidx_stSVM,svm_stSVM,
       addinets,svm_aggreadd,genesidx_netrank,svm_netrank,
       genesidx_cox,svm_coxfil,genesidx_coxlasso,svm_coxlasso,genesidx_RDS,svm_rds,survnets,svm_survnet,
       rnf_allemtfeatures,rnf_ttest,rnf_lasso,rnf_netlasso,rnf_stSVM,rnf_aggreadd,rnf_netrank,rnf_coxfil,rnf_coxlasso,rnf_rds,rnf_survnet,
       threstrainidx,nothrestrainidx)
  
}

#clinicaldata does not have sample column
# thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20
single_level_methods1_cliplus <- function(clinicaldata, moldata, threstrainidx, days,numtopfs,numtopnets){
  
  # use only clinical features
  clitrainset <- clinicaldata[threstrainidx,]
  traindays <- days[threstrainidx]
  clitestset <- clinicaldata[-threstrainidx,]
  testdays <- days[-threstrainidx]
  svm_cli <- svm1(clitrainset,traindays, clitestset, testdays)
  # lgs_cli <- logisticreg1(clitrainset,traindays, clitestset, testdays)
  rnf_cli <- randomforest1(clitrainset,traindays, clitestset, testdays)
  
  trainset <- cbind(clitrainset,moldata[threstrainidx,])
  testset <- cbind(clitestset,moldata[-threstrainidx,])
  svm_cliall <- svm1(trainset,traindays, testset, testdays)
  # lgs_clittest <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_cliall <- randomforest1(trainset,traindays,testset,testdays)
  
  # rank genes by ttest
  ttestfs <- order(rowSums(ttest_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,ttestfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,ttestfs])
  svm_clittest <- svm1(trainset,traindays, testset, testdays)
  # lgs_clittest <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clittest <- randomforest1(trainset,traindays,testset,testdays)

  #lasso
  lassofs <- order(rowSums(lasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,lassofs])
  testset <- cbind(clitestset,moldata[-threstrainidx,lassofs])
  svm_clilasso <- svm1(trainset,traindays, testset, testdays)
  # lgs_clilasso <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clilasso <- randomforest1(trainset,traindays,testset,testdays)

  #netowrk lasso
  netlassofs <- order(rowSums(netlasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,netlassofs])
  testset <- cbind(clitestset,moldata[-threstrainidx,netlassofs])
  svm_clinetlasso <- svm1(trainset,traindays, testset, testdays)
  # lgs_clinetlasso <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clinetlasso <- randomforest1(trainset,traindays,testset,testdays)
  
  #additive subnetworks
  addgfs <- addg_sf[1:numtopnets]
  directionadd <- sign(cor(moldata[threstrainidx,],traindays))
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[x])
  newX2 <- sapply(addgfs,function(x) rowSums(newX[,x])/length(x))
  trainset <- cbind(clitrainset,newX2[threstrainidx,])
  testset <- cbind(clitestset,newX2[-threstrainidx,])
  svm_cliaddg <- svm1(trainset,traindays, testset, testdays)
  # lgs_cliaddg <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_cliaddg <- randomforest1(trainset,traindays,testset,testdays)
  
  #NetRank
  netrankfs <- order(rowSums(netrank_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,netrankfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,netrankfs])
  svm_clinetrank <- svm1(trainset,traindays, testset, testdays)
  # lgs_clinetrank <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clinetrank <- randomforest1(trainset,traindays,testset,testdays)
  
  #stSVM
  stsvmfs <- order(rowSums(stsvm_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,stsvmfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,stsvmfs])
  svm_clistsvmfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clistsvmfs <- randomforest1(trainset,traindays,testset,testdays)
  
  #Cox
  coxfs <- order(rowSums(cox_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,coxfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,coxfs])
  svm_clicoxfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clicoxfs <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  coxregfs <- order(rowSums(coxlasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,coxregfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,coxregfs])
  svm_clicoxregfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clicoxregfs <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  rdsfs <- order(rowSums(rds_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,rdsfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,rdsfs])
  svm_clirdsfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clirdsfs <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  survnetfs <- survnet_sf[1:numtopnets]
  newX2 <- sapply(survnetfs,function(x) rowSums(moldata[,x])/length(x))
  trainset <- cbind(clitrainset,newX2[threstrainidx,])
  testset <- cbind(clitestset,newX2[-threstrainidx,])
  svm_clisurvnetfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clisurvnetfs <- randomforest1(trainset,traindays,testset,testdays)
  
  list(svm_cli,svm_clittest,svm_clilasso,svm_clinetlasso,svm_clistsvmfs,svm_cliaddg,svm_clinetrank,svm_clicoxfs,svm_clicoxregfs,
       svm_clirdsfs,svm_clisurvnetfs,rnf_cli,rnf_clittest,rnf_clilasso,rnf_clinetlasso,rnf_clistsvmfs,rnf_cliaddg,rnf_clinetrank,rnf_clicoxfs,rnf_clicoxregfs,
       rnf_clirdsfs,rnf_clisurvnetfs,threstrainidx,svm_cliall,rnf_cliall)
}

single_level_methods1_cliplus_2sets <- function(clinicaldata, moldata, threstrainidx, days,numtopfs,numtopnets){
  
  # use only clinical features
  clitrainset <- clinicaldata[threstrainidx,]
  traindays <- days[threstrainidx]
  clitestset <- clinicaldata[-threstrainidx,]
  testdays <- days[-threstrainidx]
  svm_cli <- svm1(clitrainset,traindays, clitestset, testdays)
  # lgs_cli <- logisticreg1(clitrainset,traindays, clitestset, testdays)
  rnf_cli <- randomforest1(clitrainset,traindays, clitestset, testdays)
  
  trainset <- cbind(clitrainset,moldata[threstrainidx,])
  testset <- cbind(clitestset,moldata[-threstrainidx,])
  svm_cliall <- svm1(trainset,traindays, testset, testdays)
  # lgs_clittest <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_cliall <- randomforest1(trainset,traindays,testset,testdays)
  
  # rank genes by ttest
  ttestfs <- order(rowSums(ttest_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,ttestfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,ttestfs])
  svm_clittest <- svm1(trainset,traindays, testset, testdays)
  # lgs_clittest <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clittest <- randomforest1(trainset,traindays,testset,testdays)
  
  #lasso
  lassofs <- order(rowSums(lasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,lassofs])
  testset <- cbind(clitestset,moldata[-threstrainidx,lassofs])
  svm_clilasso <- svm1(trainset,traindays, testset, testdays)
  # lgs_clilasso <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clilasso <- randomforest1(trainset,traindays,testset,testdays)
  
  #netowrk lasso
  netlassofs <- order(rowSums(netlasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,netlassofs])
  testset <- cbind(clitestset,moldata[-threstrainidx,netlassofs])
  svm_clinetlasso <- svm1(trainset,traindays, testset, testdays)
  # lgs_clinetlasso <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clinetlasso <- randomforest1(trainset,traindays,testset,testdays)
  
  #additive subnetworks
  addgfs <- addg_sf[1:numtopnets]
  directionadd <- sign(cor(moldata[threstrainidx,],traindays))
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[x])
  newX2 <- sapply(addgfs,function(x) rowSums(newX[,x])/length(x))
  trainset <- cbind(clitrainset,newX2[threstrainidx,])
  testset <- cbind(clitestset,newX2[-threstrainidx,])
  svm_cliaddg <- svm1(trainset,traindays, testset, testdays)
  # lgs_cliaddg <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_cliaddg <- randomforest1(trainset,traindays,testset,testdays)
  
  #NetRank
  netrankfs <- order(rowSums(netrank_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,netrankfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,netrankfs])
  svm_clinetrank <- svm1(trainset,traindays, testset, testdays)
  # lgs_clinetrank <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clinetrank <- randomforest1(trainset,traindays,testset,testdays)
  
  #stSVM
  stsvmfs <- order(rowSums(stsvm_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,stsvmfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,stsvmfs])
  svm_clistsvmfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clistsvmfs <- randomforest1(trainset,traindays,testset,testdays)
  
  #Cox
  coxfs <- order(rowSums(cox_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,coxfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,coxfs])
  svm_clicoxfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clicoxfs <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  coxregfs <- order(rowSums(coxlasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,coxregfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,coxregfs])
  svm_clicoxregfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clicoxregfs <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  rdsfs <- order(rowSums(rds_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- cbind(clitrainset,moldata[threstrainidx,rdsfs])
  testset <- cbind(clitestset,moldata[-threstrainidx,rdsfs])
  svm_clirdsfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clirdsfs <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  survnetfs <- survnet_sf[1:numtopnets]
  newX2 <- sapply(survnetfs,function(x) rowSums(moldata[,x])/length(x))
  trainset <- cbind(clitrainset,newX2[threstrainidx,])
  testset <- cbind(clitestset,newX2[-threstrainidx,])
  svm_clisurvnetfs <- svm1(trainset,traindays, testset, testdays)
  rnf_clisurvnetfs <- randomforest1(trainset,traindays,testset,testdays)
  

  trainset <- moldata[threstrainidx,]
  testset <- moldata[-threstrainidx,]
  svm_cliall2 <- svm1(trainset,traindays, testset, testdays)
  # lgs_clittest <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_cliall2 <- randomforest1(trainset,traindays,testset,testdays)
  
  # rank genes by ttest
  ttestfs <- order(rowSums(ttest_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- moldata[threstrainidx,ttestfs]
  testset <- moldata[-threstrainidx,ttestfs]
  svm_clittest2 <- svm1(trainset,traindays, testset, testdays)
  # lgs_clittest <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clittest2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #lasso
  lassofs <- order(rowSums(lasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- moldata[threstrainidx,lassofs]
  testset <- moldata[-threstrainidx,lassofs]
  svm_clilasso2 <- svm1(trainset,traindays, testset, testdays)
  # lgs_clilasso <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clilasso2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #netowrk lasso
  netlassofs <- order(rowSums(netlasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- moldata[threstrainidx,netlassofs]
  testset <- moldata[-threstrainidx,netlassofs]
  svm_clinetlasso2 <- svm1(trainset,traindays, testset, testdays)
  # lgs_clinetlasso <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clinetlasso2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #additive subnetworks
  addgfs <- addg_sf[1:numtopnets]
  directionadd <- sign(cor(moldata[threstrainidx,],traindays))
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[x])
  newX2 <- sapply(addgfs,function(x) rowSums(newX[,x])/length(x))
  trainset <- newX2[threstrainidx,]
  testset <- newX2[-threstrainidx,]
  svm_cliaddg2 <- svm1(trainset,traindays, testset, testdays)
  # lgs_cliaddg <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_cliaddg2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #NetRank
  netrankfs <- order(rowSums(netrank_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- moldata[threstrainidx,netrankfs]
  testset <- moldata[-threstrainidx,netrankfs]
  svm_clinetrank2 <- svm1(trainset,traindays, testset, testdays)
  # lgs_clinetrank <- logisticreg1(trainset,traindays,testset,testdays)
  rnf_clinetrank2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #stSVM
  stsvmfs <- order(rowSums(stsvm_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- moldata[threstrainidx,stsvmfs]
  testset <- moldata[-threstrainidx,stsvmfs]
  svm_clistsvmfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clistsvmfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #Cox
  coxfs <- order(rowSums(cox_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- moldata[threstrainidx,coxfs]
  testset <- moldata[-threstrainidx,coxfs]
  svm_clicoxfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clicoxfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  coxregfs <- order(rowSums(coxlasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- moldata[threstrainidx,coxregfs]
  testset <- moldata[-threstrainidx,coxregfs]
  svm_clicoxregfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clicoxregfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  rdsfs <- order(rowSums(rds_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- moldata[threstrainidx,rdsfs]
  testset <- moldata[-threstrainidx,rdsfs]
  svm_clirdsfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clirdsfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  survnetfs <- survnet_sf[1:numtopnets]
  newX2 <- sapply(survnetfs,function(x) rowSums(moldata[,x])/length(x))
  trainset <- newX2[threstrainidx,]
  testset <- newX2[-threstrainidx,]
  svm_clisurvnetfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clisurvnetfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  
  
  list(svm_cli,svm_clittest,svm_clilasso,svm_clinetlasso,svm_clistsvmfs,svm_cliaddg,svm_clinetrank,svm_clicoxfs,svm_clicoxregfs,
       svm_clirdsfs,svm_clisurvnetfs,rnf_cli,rnf_clittest,rnf_clilasso,rnf_clinetlasso,rnf_clistsvmfs,rnf_cliaddg,rnf_clinetrank,rnf_clicoxfs,rnf_clicoxregfs,
       rnf_clirdsfs,rnf_clisurvnetfs,threstrainidx,svm_cliall,rnf_cliall,
       svm_clittest2,svm_clilasso2,svm_clinetlasso2,svm_clistsvmfs2,svm_cliaddg2,svm_clinetrank2,svm_clicoxfs2,svm_clicoxregfs2,
       svm_clirdsfs2,svm_clisurvnetfs2,rnf_clittest2,rnf_clilasso2,rnf_clinetlasso2,rnf_clistsvmfs2,rnf_cliaddg2,rnf_clinetrank2,rnf_clicoxfs2,rnf_clicoxregfs2,
       rnf_clirdsfs2,rnf_clisurvnetfs2,svm_cliall2,rnf_cliall2)
}

evaluation_algorithms_FSF <- function(moldata, threstrainidx,numtopfs,numtopnets){
  
  matrixX <- moldata[,-c(1,ncol(threstrainidx))]
  trainset <- matrixX[threstrainidx,]
  testset <- matrixX[-threstrainidx,]
  
  traindays <- moldata[threstrainidx,ncol(moldata)]
  testdays <- moldata[-threstrainidx,ncol(moldata)]
  
  svm_all2 <- svm1(trainset,traindays, testset, testdays)
  rnf_all2 <- randomforest1(trainset,traindays,testset,testdays)
  
  # rank genes by ttest
  ttestfs <- order(rowSums(ttest_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- matrixX[threstrainidx,ttestfs]
  testset <- matrixX[-threstrainidx,ttestfs]
  svm_clittest2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clittest2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #lasso
  lassofs <- order(rowSums(lasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- matrixX[threstrainidx,lassofs]
  testset <- matrixX[-threstrainidx,lassofs]
  svm_clilasso2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clilasso2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #netowrk lasso
  netlassofs <- order(rowSums(netlasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- matrixX[threstrainidx,netlassofs]
  testset <- matrixX[-threstrainidx,netlassofs]
  svm_clinetlasso2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clinetlasso2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #additive subnetworks
  addgfs <- addg_sf[1:numtopnets]
  directionadd <- sign(cor(matrixX[threstrainidx,],traindays))
  newX <- sapply(1:ncol(matrixX), function(x) matrixX[,x]*directionadd[x])
  newX2 <- sapply(addgfs,function(x) rowSums(newX[,x])/length(x))
  trainset <- newX2[threstrainidx,]
  testset <- newX2[-threstrainidx,]
  svm_cliaddg2 <- svm1(trainset,traindays, testset, testdays)
  rnf_cliaddg2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #NetRank
  netrankfs <- order(rowSums(netrank_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- matrixX[threstrainidx,netrankfs]
  testset <- matrixX[-threstrainidx,netrankfs]
  svm_clinetrank2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clinetrank2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #stSVM
  stsvmfs <- order(rowSums(stsvm_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- matrixX[threstrainidx,stsvmfs]
  testset <- matrixX[-threstrainidx,stsvmfs]
  svm_clistsvmfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clistsvmfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #Cox
  coxfs <- order(rowSums(cox_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- matrixX[threstrainidx,coxfs]
  testset <- matrixX[-threstrainidx,coxfs]
  svm_clicoxfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clicoxfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  coxregfs <- order(rowSums(coxlasso_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- matrixX[threstrainidx,coxregfs]
  testset <- matrixX[-threstrainidx,coxregfs]
  svm_clicoxregfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clicoxregfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  rdsfs <- order(rowSums(rds_sf),decreasing = TRUE)[1:numtopfs]
  trainset <- matrixX[threstrainidx,rdsfs]
  testset <- matrixX[-threstrainidx,rdsfs]
  svm_clirdsfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clirdsfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  #Coxlasso
  survnetfs <- survnet_sf[1:numtopnets]
  newX2 <- sapply(survnetfs,function(x) rowSums(matrixX[,x])/length(x))
  trainset <- newX2[threstrainidx,]
  testset <- newX2[-threstrainidx,]
  svm_clisurvnetfs2 <- svm1(trainset,traindays, testset, testdays)
  rnf_clisurvnetfs2 <- randomforest1(trainset,traindays,testset,testdays)
  
  list(threstrainidx,svm_clittest2,svm_clilasso2,svm_clinetlasso2,svm_clistsvmfs2,svm_cliaddg2,svm_clinetrank2,svm_clicoxfs2,svm_clicoxregfs2,
       svm_clirdsfs2,svm_clisurvnetfs2,rnf_clittest2,rnf_clilasso2,rnf_clinetlasso2,rnf_clistsvmfs2,rnf_cliaddg2,rnf_clinetrank2,rnf_clicoxfs2,rnf_clicoxregfs2,
       rnf_clirdsfs2,rnf_clisurvnetfs2)
}

# a function for the statistics of clustering
cluster_measure <- function(nothresdata,numcluster){
  originalclr <- tsne_km_cl(nothresdata,1:(ncol(nothresdata)-2),0,algorithmstrings[i],numcluster)
  i=1
  ttestclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=2
  lassoclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=3
  netlassoclr <- tsne_km_cl(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster)
  i=4
  addgfs <- addg_sf[1:numtopnets]
  moldata <- nothresdata[,-c(1,ncol(nothresdata))]
  directionadd <- sign(cor(moldata,nothresdata$days))
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
  
  tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
  if(usetsne==1)
    d_tsne_1 = as.data.frame(tsne_model_1$Y) else
      d_tsne_1 = data
  
  d_tsne_1_original=d_tsne_1
  fit_cluster_kmeans=kmeans(scale(d_tsne_1), numcluster,iter.max = 100)  
  d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))
  d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=numcluster))  
  
  #for saving plot with correct names. 
  fit <- survfit(survobj ~ cl_kmeans,data=d_tsne_1_original)
  # filestring <- paste('../Thesis_Results/alllevels/clustering/', ncol(nothresdata)-2, algorithm, '_kmeans.pdf', sep='')
  # pdf(filestring, width = 6.5,height = 5)
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = d_tsne_1_original,  # data used to fit survival curves.
    risk.table = TRUE,       # show risk table.
    pval = paste('p =',format(pval_kmeans,digits=3)),             # show p-value of log-rank test.
    # pval=TRUE,
    conf.int = TRUE,         # show confidence intervals for
    # point estimaes of survival curves.
    xlim = c(0,2000),        # present narrower X axis, but not affect
    # survival estimates.
    break.time.by = 500,     # break X axis in time intervals by 500.
    ggtheme = theme_minimal(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  )
  # dev.off()
  # 
  
  fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))
  d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=numcluster))  
  fit <- survfit(survobj ~ cl_hierarchical,data=d_tsne_1_original)
  # filestring <- paste('../Thesis_Results/alllevels/clustering/', ncol(nothresdata)-2, algorithm, '_hierarchical.pdf', sep='')
  # pdf(filestring, width = 6.5,height = 5)

  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = d_tsne_1_original,  # data used to fit survival curves.
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for
    # point estimaes of survival curves.
    xlim = c(0,2000),        # present narrower X axis, but not affect
    # survival estimates.
    break.time.by = 500,     # break X axis in time intervals by 500.
    ggtheme = theme_minimal(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
   )
  # dev.off()
  
  # log-rank test
  sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,1) 
  sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_hierarchical)$chisq
  pval_hierarchical <- 1- pchisq(sdiff,1) 
  c(pval_kmeans,as.numeric(table(d_tsne_1_original$cl_kmeans)),pval_hierarchical,as.numeric(table(d_tsne_1_original$cl_kmeans)))
}
# single_level_methods2 <- function(nothresdata, nothrestrainidx){
#   
#   #---------------------------------- with survival object-----------------------------------------
#   
#   # data: LEgeneexp and LEmethy
#   # ----------------------------no network
#   trainset <- nothresdata[nothrestrainidx,]
#   testset <- nothresdata[-nothrestrainidx,]
# 
#   #dataset has sample column, be careful with index
#   #rank genes by cox p values
#   genesidx <- Coxrank(trainset,ntopgene)
#   # fcfbgenes <- genesidx[domtsinglefeatures(genesidx,trainset,1)]
#   result_coxrank <- sgsurvival(genesidx,trainset,testset)
# 
#   #lasso
#   beta <- lassocoxcoefs(trainset)
#   result_coxlasso <- sgsurvivalvec(beta,testset)
# 
#   #superpc
#   superpcvec <- superpc_vec(trainset,testset)
#   result_superpc_net <- survivalpc(superpcvec,testset)
#   superpcvec_sg <- superpc_vec_sg(trainset,testset)
#   result_superpc_sg <- survivalpc(superpcvec_sg,testset)
# 
#   #RDS
#   genesidx <- RDS_RGS(trainset,ntopgene)
#   result_RDS <- sgsurvival(genesidx,trainset,testset)
# 
#   #Survnet
#   netrank_survnet <-survnet(trainset,ntopnet)
#   resultaggre_survnet <-netsurvivalaggre(netrank_survnet[[1]],trainset,testset)
#   resultpc_surnet <-netsurvivalpc(netrank_survnet[[1]],trainset,testset)
# 
#   #WGCNA
#   subnets <- WGCNAnet(trainset,ntopnet)
#   resultaggre_wgcna <-netsurvivalaggre(subnets,trainset,testset)
#   resultpc_wgcna <-netsurvivalpc(subnets,trainset,testset)
# 
#   # list(threstrainidx, nothrestrainidx,result_allfeatures,result_ttest,result_logicreg,result_lasso,result_netlasso,resultaggre_addv,resultpc_addv,result_netrank,result_coxrank,result_superpc_net,result_superpc_sg,result_RDS,resultaggre_survnet,resultpc_surnet,resultaggre_wgcna,resultpc_wgcna)
#   list(result_coxrank,result_coxlasso,result_superpc_net,result_superpc_sg,result_RDS,resultaggre_survnet,resultpc_surnet,resultaggre_wgcna,resultpc_wgcna)
#   
# }
# 
# 
# # LEthresgeneexp <- LEthresgeneexp[sample(nrow(LEthresgeneexp)),]
# # LEgeneexp <- LEgeneexp[sample(nrow(LEgeneexp)),]
# X <- LEthresgeneexp[,-c(1,ncol(LEthresgeneexp))]
# LEthresgeneexp[,-c(1,ncol(LEthresgeneexp))] <- normalize_data(normalize_data(X,2),12)
# X <- LEgeneexp[,-c(1,ncol(LEgeneexp))]
# LEgeneexp[,-c(1,ncol(LEgeneexp))] <- normalize_data(normalize_data(X,2),12)
# 
# flds1 <- sapply(1:10, function(x) createFolds(LEthresgeneexp$days, k = 10, list = TRUE, returnTrain = TRUE))
# registerDoMC(cores=10)
# MCCV_Geneexp1 <- foreach(x = 1:length(flds1),.errorhandling='remove')  %dopar% {
#                                threstrainidx <- sampling(rmoutlier,0.9)
#                               single_level_methods1(rmoutlier,threstrainidx)
# }
# 
# nothres_flds <- sapply(1:10, function(x) createFolds(LEgeneexp$days, k = 10, list = TRUE, returnTrain = TRUE))
# MCCV_Geneexp2 <- foreach(x = 1:length(nothres_flds),.errorhandling='remove')  %dopar% {
#   nothrestrainidx <-  nothres_flds[[x]]
#   single_level_methods2(LEgeneexp,nothrestrainidx)
# }
# save.image('LUAD_Single_Level_Results100.RData')
# 
# 

# LEthresmethy <- LEthresmethy[sample(nrow(LEthresmethy)),]
# LEmethy <- LEmethy[sample(nrow(LEmethy)),]
# svm1 <- function(traindata,trainy,testdata,testy){
#   X <- traindata
#   y <- trainy
#   svm.model <- svm(X,as.factor(y),probability=TRUE, gamma = 2^-9, cost = 2)
#   svm.pred <- predict(svm.model, testdata,probability=TRUE)
#   pred <- prediction(attr(svm.pred,"prob")[,1],as.numeric(testy))
#   roc_obj <- roc(as.numeric(testy), attr(svm.pred,"prob")[,1])
#   auc_obj <- as.numeric(roc_obj$auc)
#   prec <- unlist(performance(pred,"prec")@y.values)
#   rec <- unlist(performance(pred,"rec")@y.values)
#   if(is.na(prec[1])) 
#     prec[1] <- 1
#   height = (prec[-1]+prec[-length(prec)])/2
#   width = diff(rec)
#   aupr <- sum(height*width)
#   c(auc_obj,aupr)
# }
# flds1 <- sapply(1:10, function(x) createFolds(LEthresmethy$days, k = 10, list = TRUE, returnTrain = TRUE))
# registerDoMC(cores=10)
# MCCV_Methy1 <- foreach(x = 1:length(flds1),.errorhandling='remove')  %dopar% {
#   threstrainidx <- flds1[[x]]
#   single_level_methods1(LEthresmethy,threstrainidx)
# }
# flds2 <- sapply(1:10, function(x) createFolds(LEmethy$days, k = 10, list = TRUE, returnTrain = TRUE))
# MCCV_Methy2 <- foreach(x = 1:10,.errorhandling='remove')  %dopar% {
#   nothrestrainidx <- flds2[[x]]
#   single_level_methods2(LEmethy,nothrestrainidx)
# }
# save.image('LUAD_Single_Level_Results100.RData')
# 
# par(mfrow=c(1,1))
# plot(1:200,tempttest3,type="l",main="AUC values on training set VS. # features")
# points(1:200,tempttest3,pch=as.character(1:200))
# # dev.off()
# 
# sapply(seq(0.1,2,0.1), function(x) netclassaggre(subnetworks[domtnetfeatures(subnetworks,trainset,x)],trainset,testset,2))
# sapply(seq(0.1,2,0.1), function(x) sgclass(genesidx[domtsinglefeatures(genesidx,trainset,x)],trainset,testset,2))
# 
# 
# tempttest3 <- foreach(x = 1:30,.errorhandling='remove') %dopar% {
#   fcfbnets <- subnetworks[domtnetfeatures(subnetworks,trainset,1)]
#   subnetidx <- order(scores,decreasing = TRUE)[1:x]
#   netclassaggre(subnetworks[subnetidx],trainset,testset,1)[1]
# }
# plot(1:30,tempttest3,type="l",main="AUC values on testing set VS. # features")
# points(1:30,tempttest3,pch=as.character(1:30))
# 

# save the trainindex to a file for Matlab
# threshold data, take the first element of the list from each run
# gene expression data
# thres_trainidx_geneexp <- sapply(1:length(MCCVresults_Geneexp), function(x) MCCVresults_Geneexp[[x]][[1]])
# write.table(thres_trainidx_geneexp, file="thres_trainidx_geneexp.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, sep=",")
# 
# nothres_trainidx_geneexp <- sapply(1:length(MCCVresults_Geneexp), function(x) MCCVresults_Geneexp[[x]][[2]])
# write.table(nothres_trainidx_geneexp, file="nothres_trainidx_geneexp.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, sep=",")
# 
# #methylation data 
# thres_trainidx_methy <- sapply(1:length(MCCVresults_Methy), function(x) MCCVresults_Methy[[x]][[1]])
# write.table(thres_trainidx_methy, file="thres_trainidx_methy.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, sep=",")
# 
# nothres_trainidx_methy <- sapply(1:length(MCCVresults_Methy), function(x) MCCVresults_Methy[[x]][[2]])
# write.table(nothres_trainidx_methy, file="nothres_trainidx_methy.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, sep=",")

# write dataa to files according to the matlab code, dataset and network info

# read in matlab output and perform analysis.

# check results
# sapply(1:100, function(x) sapply(3:4, function(y) MCCVresults_Geneexp[[x]][[y]][1]))

#order the output so that the same index gets the same value, e.g., auc value, p value. 

