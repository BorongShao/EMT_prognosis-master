repetition <- 300

# registerDoMC(cores=10)
# LUAD_geneexp_flds_700_1400_thres <- sapply(1:(repetition/10), function(x) createFolds(thresdata$days, k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_geneexp_flds_700_1400_nothres <- sapply(1:(repetition/10), function(x) createFolds(nothresdata$days, k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_methy_flds_700_1400_thres <- sapply(1:(repetition/10), function(x) createFolds(thresdata$days, k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_methy_flds_700_1400_nothres <- sapply(1:(repetition/10), function(x) createFolds(nothresdata$days, k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_cna_flds_700_1400_thres <- sapply(1:(repetition/10), function(x) createFolds(thresdata$days, k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_cna_flds_700_1400_nothres <- sapply(1:(repetition/10), function(x) createFolds(nothresdata$days, k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_geneexp_cliflds700_1400 <- sapply(1:(repetition/10), function(x) createFolds(thresdata2[[3]], k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_methy_cliflds700_1400 <- sapply(1:(repetition/10), function(x) createFolds(thresdata2[[3]], k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_cna_cliflds700_1400 <- sapply(1:(repetition/10), function(x) createFolds(thresdata2[[3]], k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_methy_flds_3y_thres <- sapply(1:(repetition/10), function(x) createFolds(thresdata$days, k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_methy_cliflds3y <- sapply(1:(repetition/10), function(x) createFolds(thresdata2[[3]], k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_methy_flds_5y_thres <- sapply(1:(repetition/10), function(x) createFolds(thresdata$days, k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_methy_cliflds5y <- sapply(list, function), function(x) createFolds(thresdata2[[3]], k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_methy_flds_900_1200_thres <- sapply(1:(repetition/10), function(x) createFolds(thresdata$days, k = 10, list = TRUE, returnTrain = TRUE))
# LUAD_methy_cliflds_900_1200 <- sapply(1:(repetition/10), function(x) createFolds(thresdata2[[3]], k = 10, list = TRUE, returnTrain = TRUE))
# algorithmstrings <- c("ttestfs","lassofs","netlassofs", "addgfs", "netrankfs", "stsvmfs", "coxfs", "coxregfs", "rdsfs", "survnetfs","alg10sfs","nofs")
# adddirection_GE74 <- sign(cor(thresdata[,-c(1,ncol(thresdata))],thresdata$days))
# adddirection_GE123 <- sign(cor(thresdata[,-c(1,ncol(thresdata))],thresdata$days))
# adddirection_GE455 <- sign(cor(thresdata[,-c(1,ncol(thresdata))],thresdata$days))
# adddirection_DM74 <- sign(cor(thresdata[,-c(1,ncol(thresdata))],thresdata$days))
# adddirection_DM123 <- sign(cor(thresdata[,-c(1,ncol(thresdata))],thresdata$days))
# adddirection_DM455 <- sign(cor(thresdata[,-c(1,ncol(thresdata))],thresdata$days))
# adddirection_CA74 <- sign(cor(thresdata[,-c(1,ncol(thresdata))],thresdata$days))
# adddirection_CA455 <- sign(cor(thresdata[,-c(1,ncol(thresdata))],thresdata$days))
# add_dire74 <- list(adddirection_GE74, adddirection_DM74,adddirection_CA74)
# add_dire123 <- list(adddirection_GE123, adddirection_DM123,adddirection_CA123)
# add_dire455 <- list(adddirection_GE455, adddirection_DM455,adddirection_CA455)

LEnet <- LEnet74
numfeatures <- 74
thresdata <- read.csv('../Thesis_Results/74network/emt74_geneexp_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/74network/emt74_geneexp_700_1400_nothres.csv')
# LUAD74_geneexp_700_1400<- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_geneexp_flds_700_1400_thres[[x]]
#   nothrestrainidx <- LUAD_geneexp_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }

thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD74_geneexp_700_1400
source('selectedfeatures_clustering.R')
LUAD74_geneexp_cluster2 <- cluster_measure(nothresdata,2,add_dire74[[1]])
LUAD74_geneexp_cluster3 <- cluster_measure(nothresdata,3,add_dire74[[1]])
LUAD74_cligeneexp_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_geneexp_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}


thresdata <- read.csv('../Thesis_Results/74network/emt74_methy_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/74network/emt74_methy_700_1400_nothres.csv')
# LUAD74_methy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_methy_flds_700_1400_thres[[x]]
#   nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }

thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD74_methy_700_1400
source('selectedfeatures_clustering.R')
LUAD74_methy_cluster2 <- cluster_measure(nothresdata,2,add_dire74[[2]])
LUAD74_methy_cluster3 <- cluster_measure(nothresdata,3,add_dire74[[2]])
LUAD74_climethy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}


LEnet <- LEnet74_cna
numfeatures <- 70
thresdata <- read.csv('../Thesis_Results/74network/emt70_cna_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/74network/emt70_cna_700_1400_nothres.csv')
# LUAD74_cna_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_cna_flds_700_1400_thres[[x]]
#   nothrestrainidx <- LUAD_cna_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }

thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD74_cna_700_1400
source('selectedfeatures_clustering.R')
LUAD74_cna_cluster2 <- cluster_measure(nothresdata,2,add_dire74[[3]])
LUAD74_cna_cluster3 <- cluster_measure(nothresdata,3,add_dire74[[3]])
LUAD74_clicna_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_cna_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}



LEnet <- LEnet123
numfeatures <- 123
thresdata <- read.csv('../Thesis_Results/123network/emt123_geneexp_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/123network/emt123_geneexp_700_1400_nothres.csv')
# LUAD123_geneexp_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_geneexp_flds_700_1400_thres[[x]]
#   nothrestrainidx <- LUAD_geneexp_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }

thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD123_geneexp_700_1400
source('selectedfeatures_clustering.R')
LUAD123_geneexp_cluster2 <- cluster_measure(nothresdata,2,add_dire123[[1]])
LUAD123_geneexp_cluster3 <- cluster_measure(nothresdata,3,add_dire123[[1]])
LUAD123_cligeneexp_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_geneexp_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}

thresdata <- read.csv('../Thesis_Results/123network/emt123_methy_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/123network/emt123_methy_700_1400_nothres.csv')
# LUAD123_methy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_methy_flds_700_1400_thres[[x]]
#   nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }
thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD123_methy_700_1400
source('selectedfeatures_clustering.R')
LUAD123_methy_cluster2 <- cluster_measure(nothresdata,2,add_dire123[[2]])
LUAD123_methy_cluster3 <- cluster_measure(nothresdata,3,add_dire123[[2]])
LUAD123_climethy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}

LEnet <- LEnet123_cna
numfeatures <- 117
thresdata <- read.csv('../Thesis_Results/123network/emt117_cna_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/123network/emt117_cna_700_1400_nothres.csv')
# LUAD123_cna_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_cna_flds_700_1400_thres[[x]]
#   nothrestrainidx <- LUAD_cna_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }

thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD123_cna_700_1400
source('selectedfeatures_clustering.R')
LUAD123_cna_cluster2 <- cluster_measure(nothresdata,2,add_dire123[[3]])
LUAD123_cna_cluster3 <- cluster_measure(nothresdata,3,add_dire123[[3]])
LUAD123_clicna_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_cna_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}


# LUAD_folds_clinical_700_1400 <- clinical_svm()
# LUAD_folds_clinical_900_1200 <- clinical_svm()
# LUAD_folds_clinical_500_1500 <- clinical_svm()
# LUAD_folds_clinical_3y <- clinical_svm()

 
LEnet <- LEnet455
numfeatures <- 455
thresdata <- read.csv('../Thesis_Results/455network/emt455_geneexp_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/455network/emt455_geneexp_700_1400_nothres.csv')
# LUAD455_geneexp_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_geneexp_flds_700_1400_thres[[x]]
#   nothrestrainidx <- LUAD_geneexp_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }
thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD455_geneexp_700_1400
source('selectedfeatures_clustering.R')
LUAD455_FSFgeneexp_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_geneexp_flds_700_1400_thres[[x]]
  nothrestrainidx <- LUAD_geneexp_flds_700_1400_nothres[[x]]
  evaluation_algorithms_FSF(thresdata,threstrainidx)
}
LUAD455_geneexp_cluster2 <- cluster_measure(nothresdata,2,add_dire455[[1]])
LUAD455_geneexp_cluster3 <- cluster_measure(nothresdata,3,add_dire455[[1]])
LUAD455_cligeneexp_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_geneexp_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}

thresdata <- read.csv('../Thesis_Results/455network/emt455_methy_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/455network/emt455_methy_700_1400_nothres.csv')
# LUAD455_methy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_methy_flds_700_1400_thres[[x]]
#   nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }
thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD455_methy_700_1400
source('selectedfeatures_clustering.R')
LUAD455_FSFmethy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_flds_700_1400_thres[[x]]
  nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
  evaluation_algorithms_FSF(thresdata,threstrainidx,20,10)
}
LUAD455_methy_cluster2 <- cluster_measure(nothresdata,2,add_dire455[[2]])
LUAD455_methy_cluster3 <- cluster_measure(nothresdata,3,add_dire455[[2]])
LUAD455_climethy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}



LEnet <- LEnet455_cna
numfeatures <- 445
thresdata <- read.csv('../Thesis_Results/455network/emt445_cna_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/455network/emt445_cna_700_1400_nothres.csv')
# LUAD455_cna_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_cna_flds_700_1400_thres[[x]]
#   nothrestrainidx <- LUAD_cna_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }

thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD455_cna_700_1400
source('selectedfeatures_clustering.R')
LUAD455_cna_cluster2 <- cluster_measure(nothresdata,2,add_dire455[[3]])
LUAD455_cna_cluster3 <- cluster_measure(nothresdata,3,add_dire455[[3]])
LUAD455_clicna_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_cna_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}

# thresdata <- read.csv('../Thesis_Results/74network/emt74_geneexp_700_1400.csv')
# tokeep1 <- match(thresdata[,1],LErnaseq_all[,1])
# tokeep2 <- match(thresdata[,1],LEmirna_all[,1])
# allrnaseq <- data.frame(thresdata[,1],LErnaseq_all[tokeep1,-1],LEmirna_all[tokeep2,-1],thresdata$days)
# rmfeatures <- which(sapply(allrnaseq[,-c(1,ncol(allrnaseq))], var)<0.5)+1
# allrnaseq <- allrnaseq[,-rmfeatures]
# names(allrnaseq)[ncol(allrnaseq)] <- "days"
# LUADlevel1all_700_1400<- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_geneexp_flds_700_1400_thres[[x]]
#   trainset <- allrnaseq[threstrainidx,]
#   testset <- allrnaseq[-threstrainidx,]
#   geneindex <- lassocoefs(trainset) + 1
#   traindata <- data.matrix(trainset[,geneindex])
#   trainy <- trainset$days
#   testdata <- data.matrix(testset[,geneindex])
#   testy <- testset$days
#   s1 <- svm1(traindata,trainy,testdata,testy)
#   s2 <- svm1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
#   r1 <- randomforest1(traindata,trainy,testdata,testy)
#   r2 <- randomforest1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
#   list(s1,s2,r1,r2,geneindex)
# }
# 
# 
# thresdata <- read.csv('../Thesis_Results/74network/emt74_methy_700_1400.csv')
# tokeep <- match(thresdata[,1],LEmethy_all[,1])
# allmethy <- data.frame(thresdata[,1],LEmethy_all[tokeep,-1],thresdata$days)
# names(allmethy)[ncol(allmethy)] <- "days"
# allmethy <- allmethy[,sapply(allmethy,function(x) sum(is.na(x))==0)]
# LUADmethyall_700_1400<- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_methy_flds_700_1400_thres[[x]]
#   trainset <- allmethy[threstrainidx,]
#   testset <- allmethy[-threstrainidx,]
#   geneindex <- lassocoefs(trainset) + 1
#   traindata <- data.matrix(trainset[,geneindex])
#   trainy <- trainset$days
#   testdata <- data.matrix(testset[,geneindex])
#   testy <- testset$days
#   s1 <- svm1(traindata,trainy,testdata,testy)
#   s2 <- svm1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
#   r1 <- randomforest1(traindata,trainy,testdata,testy)
#   r2 <- randomforest1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
#   list(s1,s2,r1,r2,geneindex)
# }
# 
# 
# thresdata <- read.csv('../Thesis_Results/74network/emt70_cna_700_1400.csv')
# tokeep <- match(thresdata[,1],LEcna_all[,1])
# allcna <- data.frame(thresdata[,1],LEcna_all[tokeep,-1],thresdata$days)
# names(allcna)[ncol(allcna)] <- "days"
# allcna <- allcna[,sapply(allcna,function(x) sum(is.na(x))==0)]
# LUADcnaall_700_1400<- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_cna_flds_700_1400_thres[[x]]
#   trainset <- allcna[threstrainidx,]
#   testset <- allcna[-threstrainidx,]
#   geneindex <- lassocoefs(trainset) + 1
#   traindata <- data.matrix(trainset[,geneindex])
#   trainy <- trainset$days
#   testdata <- data.matrix(testset[,geneindex])
#   testy <- testset$days
#   s1 <- svm1(traindata,trainy,testdata,testy)
#   s2 <- svm1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
#   r1 <- randomforest1(traindata,trainy,testdata,testy)
#   r2 <- randomforest1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
#   list(s1,s2,r1,r2,geneindex)
# }

#use random networks
# take 30 sets of random features, test each of them using 10 folds CV, gather the results and plot. including all algorithms.

#methylation, 74nodes network
# LEnet <- LEnet74
# thresdata<- read.csv('../Thesis_Results/74network/emt74_methy_700_1400.csv')
# # nothresdata <- read.csv('../Thesis_Results/74network/emt74_methy_700_1400_nothres.csv')
# resultlist <- list()
# for(j in 1:5){
#   for(i in 1:30){
#     
#     # take the same samples
#     randidx <- sample(ncol(allmethy)-2,74)+1
#     randthresdata <- allmethy[match(thresdata[,1], allmethy[,1]),c(1,randidx)]
#     randthresdata <- data.frame(randthresdata,thresdata$days)
#     names(randthresdata)[ncol(randthresdata)] <- "days"
#     # randnothresdata <- allmethy[match(nothresdata[,1], allmethy[,1]),c(1,randidx)]
#     # randnothresdata <- data.frame(nothresdata[,1],randnothresdata,nothresdata$days)
#     # names(randnothresdata)[ncol(randnothresdata)] <- "days"
#     
#     temp <- foreach(x = ((i-1)*10+1):(i*10),.errorhandling='remove')  %dopar% {
#       threstrainidx <-LUAD_methy_flds_700_1400_thres[[x]]
#       # nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
#       single_level_methods1(randthresdata,threstrainidx)
#     }
#     
#     resultlist <- c(resultlist,temp)
#   }
# }
# randomnetwork_74_methy <- resultlist
# 
# 
# #gene expression, 123 nodes network
# LEnet <- LEnet123
# thresdata <- read.csv('../Thesis_Results/123network/emt123_geneexp_700_1400.csv')
# # nothresdata <- read.csv('../Thesis_Results/123network/emt123_geneexp_700_1400_nothres.csv')
# resultlist <- list()
# for(j in 1:5){
#   for(i in 1:30){
#     
#     # take the same samples
#     randidx <- sample(ncol(allrnaseq)-2,123)+1
#     randthresdata <- allrnaseq[match(thresdata[,1], allrnaseq[,1]),c(1,randidx)]
#     randthresdata <- data.frame(randthresdata,thresdata$days)
#     names(randthresdata)[ncol(randthresdata)] <- "days"
#     # randnothresdata <- allrnaseq[match(nothresdata[,1], allrnaseq[,1]),c(1,randidx)]
#     # randnothresdata <- data.frame(nothresdata[,1],randnothresdata,nothresdata$days)
#     # names(randnothresdata)[ncol(randnothresdata)] <- "days"
#     
#     temp <- foreach(x = ((i-1)*10+1):(i*10),.errorhandling='remove')  %dopar% {
#       threstrainidx <-LUAD_geneexp_flds_700_1400_thres[[x]]
#       # nothrestrainidx <- LUAD_geneexp_flds_700_1400_nothres[[x]]
#       single_level_methods1(randthresdata,threstrainidx)
#     }
#     
#     resultlist <- c(resultlist,temp)
#   }
# }
# randomnetwork_123_geneexp <- resultlist
# 
# save.image('LUADplusRNF.RData')

#kernel dens

#experiment with different threshold for one case for each data type, pick the case after the above experiments.

LEnet <- LEnet74
numfeatures <- 74
thresdata <- read.csv('../Thesis_Results/74network/emt74_methy3y.csv')
nothresdata <- read.csv('../Thesis_Results/74network/emt74_methy_700_1400_nothres.csv')
# LUAD74_methy_3y <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_methy_flds_3y_thres[[x]]
#   nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }

thresdata2 <- getclinicaldata(3*365, 3*365,thresdata)
cpy <- LUAD74_methy_3y
source('selectedfeatures_clustering.R')
LUAD74_climethy_3y <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_cliflds3y[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}


thresdata <- read.csv('../Thesis_Results/74network/emt74_methy5y.csv')
# LUAD74_methy_5y <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_methy_flds_5y_thres[[x]]
#   nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }

thresdata2 <- getclinicaldata(1500, 500,thresdata)
cpy <- LUAD74_methy_5y
source('selectedfeatures_clustering.R')
LUAD74_climethy_5y <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_cliflds5y[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}

thresdata <- read.csv('../Thesis_Results/74network/emt74_methy.csv')
# LUAD74_methy_900_1200 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
#   threstrainidx <-LUAD_methy_flds_900_1200_thres[[x]]
#   nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
#   evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
# }

thresdata2 <- getclinicaldata(1200, 900,thresdata)
cpy <- LUAD74_methy_900_1200
source('selectedfeatures_clustering.R')
LUAD74_climethy_900_1200 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_cliflds_900_1200[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}

save.image('LUADplusRNF.RData')


# save.image('LUAD_with_rules2.RData')
# initialization for plots
# auc_all_geneexp <- sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[2]][[1]])
# aupr_all_geneexp <- sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[2]][[2]])
# accu_all_geneexp <- rowSums(sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[2]][[3]]))/length(LUADlevel1all_700_1400)
# 
# auc_allpluslasso_geneexp <- sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[1]][[1]])
# aupr_allpluslasso_geneexp <- sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[1]][[2]])
# accu_allpluslasso_geneexp <- rowSums(sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[1]][[3]]))/length(LUADlevel1all_700_1400)
# 
# 
# auc_all_methy <- sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[2]][[1]])
# aupr_all_methy <- sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[2]][[2]])
# accu_all_methy <- rowSums(sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[2]][[3]]))/length(LUADmethyall_700_1400)
# 
# auc_allpluslasso_methy <- sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[1]][[1]])
# aupr_allpluslasso_methy <- sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[1]][[2]])
# accu_allpluslasso_methy <- rowSums(sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[1]][[3]]))/length(LUADmethyall_700_1400)
# 
# auc_all_cna <- sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[2]][[1]])
# aupr_all_cna <- sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[2]][[2]])
# accu_all_cna <- rowSums(sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[2]][[3]]))/length(LUADcnaall_700_1400)
# 
# auc_allpluslasso_cna <- sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[1]][[1]])
# aupr_allpluslasso_cna <- sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[1]][[2]])
# accu_allpluslasso_cna <- rowSums(sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[1]][[3]]))/length(LUADcnaall_700_1400)
# 
# auc_clinical <- sapply(1:repetition, function(x) LUAD_folds_clinical_700_1400[[x]][[1]][[1]])
# aupr_clinical <- sapply(1:repetition, function(x) LUAD_folds_clinical_700_1400[[x]][[1]][[2]])
# accu_clinical <- rowSums(sapply(1:repetition, function(x) LUAD_folds_clinical_700_1400[[x]][[1]][[3]]))/length(LUAD_folds_clinical_700_1400)
# 
# 
# 
# auc_all_geneexp_rnf <- sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[4]][[1]])
# aupr_all_geneexp_rnf <- sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[4]][[2]])
# accu_all_geneexp_rnf <- rowSums(sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[4]][[3]]))/length(LUADlevel1all_700_1400)
# 
# auc_allpluslasso_geneexp_rnf <- sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[3]][[1]])
# aupr_allpluslasso_geneexp_rnf <- sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[3]][[2]])
# accu_allpluslasso_geneexp_rnf <- rowSums(sapply(1:repetition, function(x) LUADlevel1all_700_1400[[x]][[3]][[3]]))/length(LUADlevel1all_700_1400)
# 
# 
# auc_all_methy_rnf <- sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[4]][[1]])
# aupr_all_methy_rnf <- sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[4]][[2]])
# accu_all_methy_rnf <- rowSums(sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[4]][[3]]))/length(LUADmethyall_700_1400)
# 
# auc_allpluslasso_methy_rnf <- sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[3]][[1]])
# aupr_allpluslasso_methy_rnf <- sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[3]][[2]])
# accu_allpluslasso_methy_rnf <- rowSums(sapply(1:repetition, function(x) LUADmethyall_700_1400[[x]][[3]][[3]]))/length(LUADmethyall_700_1400)
# 
# auc_all_cna_rnf <- sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[4]][[1]])
# aupr_all_cna_rnf <- sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[4]][[2]])
# accu_all_cna_rnf <- rowSums(sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[4]][[3]]))/length(LUADcnaall_700_1400)
# 
# auc_allpluslasso_cna_rnf <- sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[3]][[1]])
# aupr_allpluslasso_cna_rnf <- sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[3]][[2]])
# accu_allpluslasso_cna_rnf <- rowSums(sapply(1:repetition, function(x) LUADcnaall_700_1400[[x]][[3]][[3]]))/length(LUADcnaall_700_1400)
# 
# auc_clinical_rnf <- sapply(1:repetition, function(x) LUAD_folds_clinical_700_1400[[x]][[2]][[1]])
# aupr_clinical_rnf <- sapply(1:repetition, function(x) LUAD_folds_clinical_700_1400[[x]][[2]][[2]])
# accu_clinical_rnf <- rowSums(sapply(1:repetition, function(x) LUAD_folds_clinical_700_1400[[x]][[2]][[3]]))/length(LUAD_folds_clinical_700_1400)
# 

