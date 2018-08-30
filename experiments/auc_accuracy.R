setwd('/home/borong/Desktop/Experiments_largeEMT/')
load('LUAD_Results.RData')
repetition <- 100
# toremove <- which(sapply(LErnaseq_all[,-1], var) < 0.5)+1
# allrnaseqfeatures <- data.frame(LErnaseq_all[match(rmoutlier[,1], LErnaseq_all[,1]),-toremove], rmoutlier$days)
# names(allrnaseqfeatures)[ncol(allrnaseqfeatures)] <- "days"
# allrnaseqfeatures[,-c(1,ncol(allrnaseqfeatures))] <- log2(allrnaseqfeatures[,-c(1,ncol(allrnaseqfeatures))]+1)

registerDoMC(cores=10)
flds1 <- sapply(1:(repetition/10), function(x) createFolds(rmoutlier$days, k = 10, list = TRUE, returnTrain = TRUE))
commondata1 <- flds1
# flds2 <- sapply(1:(repetition/10), function(x) createFolds(LEthresmethy$days, k = 10, list = TRUE, returnTrain = TRUE))
# commondata2 <- flds2
LUAD118_MCCV_rmoutlier_50_50 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-flds1[[x]]
  single_level_methods1(rmoutlier,threstrainidx)
}
save.image('LUAD_Results118.RData')
single_level_methods1 <- function(thresdata, threstrainidx){
  
  # EMT
  trainset <- thresdata[threstrainidx,]
  testset <- thresdata[-threstrainidx,]
  # use all EMT features
  svm_allemtfeatures <- svm1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  lgs_allemtfeatures <- logisticreg1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  rnf_allemtfeatures <- randomforest1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  
  list(svm_allemtfeatures,lgs_allemtfeatures,rnf_allemtfeatures)
}


tempfolds <- c(32,149,105,197,144)
#32 bad svm, 149 good svm
#197 bad rnf, 144 good rnf.

goodfold <- flds1[[87]]
badfold <- flds1[[14]]
idx <- goodfold
trainset <- rmoutlier[idx,]
testset <- rmoutlier[-idx,]
# use all EMT features
# svm1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
X <- trainset[,-c(1,ncol(trainset))]
y <- trainset$days
svm.model <- svm(X,as.factor(y),probability=TRUE)
svm.pred <- predict(svm.model, trainset[,-c(1,ncol(trainset))],probability=TRUE)
if(trainset$days[1]==TRUE){
  predvec <- attr(svm.pred,"prob")[,1] 
}else{
  predvec <- attr(svm.pred,"prob")[,2] 
}
pred <- prediction(predvec,as.numeric(trainset$days))
# roc_obj <- roc(as.numeric(testset$days), attr(svm.pred,"prob")[,1],smooth = TRUE, smooth.method = 'binormal')
plot(performance(pred, measure = "acc"), main='training set')
roc_obj <- roc(as.numeric(trainset$days), predvec)
plot(roc_obj,print.thres=TRUE,print.auc=TRUE)

svm.pred <- predict(svm.model, testset[,-c(1,ncol(testset))],probability=TRUE)
if(trainset$days[1]==TRUE){
  predvec <- attr(svm.pred,"prob")[,1] 
}else{
  predvec <- attr(svm.pred,"prob")[,2] 
}
pred <- prediction(predvec,as.numeric(testset$days))
# roc_obj <- roc(as.numeric(testset$days), attr(svm.pred,"prob")[,1],smooth = TRUE, smooth.method = 'binormal')
plot(performance(pred, measure = "acc"),main="testing set")
roc_obj <- roc(as.numeric(testset$days), predvec)
auc <- performance(pred,measure="auc")
plot(roc_obj,print.thres=TRUE,print.auc=TRUE)

# aupr <- pr.curve(scores.class0 = svm.pred, weights.class0 = y)$auc.integral
aupr <- pr.curve(scores.class0 = roc_obj$predictor, weights.class0 = y,pr.curve=TRUE)
plot((1-aupr$curve[,1]),aupr$curve[,2],ylim=c(0,1))
grid.arrange(plot2, plot2, nrow=1, ncol=2)

tsne_model_1 = Rtsne(as.matrix(trainset[,-c(1,ncol(trainset))]), check_duplicates=FALSE, pca=TRUE, perplexity=5, theta=0.5, dims=2)
d_tsne_1 = as.data.frame(tsne_model_1$Y) 
# ggplot(d_tsne_1, aes(x=V1, y=V2)) +  
#   geom_point(size=0.25) +
#   guides(colour=guide_legend(override.aes=list(size=6))) +
#   xlab("") + ylab("") +
#   ggtitle("t-SNE") +
#   theme_light(base_size=20) +
#   theme(axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   scale_colour_brewer(palette = "Set2")

d_tsne_1_original=d_tsne_1

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_1), 2)  
d_tsne_1_original$training_labels = factor(fit_cluster_kmeans$cluster)

## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))

## setting 3 clusters as output
d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=2))  
#Plotting the cluster models onto t-SNE output

#Now time to plot the result of each cluster model, based on the t-SNE map.

plot_cluster=function(data, var_cluster, palette)  
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=0.4) +
    guides(colour=guide_legend(override.aes=list(size=2))) +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=10) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) 
}


plot_k1=plot_cluster(d_tsne_1_original, "training_labels", "Set1")  
plot_h1=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")

## and finally: putting the plots side by side with gridExtra lib...
library(gridExtra)  
grid.arrange(plot_k1, plot_h1, plot_k2, plot_h2,nrow=2, ncol=2)  



allvars <- sapply(allrnaseqfeatures[,-c(1,ncol(allrnaseqfeatures))], function(x) var(x))
allidx <- order(allvars,decreasing = TRUE)[1:118]+1
allrnaseqfeatures2 <- allrnaseqfeatures[,-c(1,allidx,ncol(allrnaseqfeatures))]

par(mfrow=c(5,2))
tempfolds <- c(95,57,65,43,22)
#32 bad svm, 149 good svm
#197 bad rnf, 144 good rnf.
for(i in 1:length(tempfolds)){
  threstrainidx <- flds1[[tempfolds[i]]]
  trainset <- thresdata[threstrainidx,]
  testset <- thresdata[-threstrainidx,]
  # genesidx <- NetRank(trainset,ntopgene)
  # plot(sapply(1:50, function(x) sgclass(genesidx[1:x],trainset,testset,2)[1]), ylab="AUC", xlab="#features")
  addinets <- additive_subnet(trainset,ntopnet)
  # fcfbnets <- addinets[domtnetfeatures(addinets,trainset,2)]
  # plot(sapply(1:length(addinets), function(x) netclassaggre(addinets[1:x], trainset,trainset,2)[1]),ylab="AUC", xlab="#features")
  plot(1:length(addinets), unlist(sapply(1:length(addinets), function(x) netclassaggre(addinets[1:x], trainset,testset,2)[1])),ylab="AUC", xlab="#features")
  plot(1:length(addinets), unlist(sapply(1:length(addinets), function(x) netclassaggre(addinets[1:x], trainset,testset,3)[1])),ylab="AUC", xlab="#features")
  # 
  # genesidx <- NetRank(trainset,70)
  # plot(sapply(1:length(genesidx), function(x) sgclass(genesidx[1:x],trainset,trainset,2)[1]),ylab="AUC", xlab="#features")
  # plot(sapply(1:length(genesidx), function(x) sgclass(genesidx[1:x],trainset,trainset,3)[1]),ylab="AUC", xlab="#features")
  # netlist <- addinets[c(1:4,6)]
  # traindata <- sapply(netlist,function(x) rowSums(trainset[,(x+1)])/length(x))
  # trainy <- trainset$days
  # testdata <- sapply(netlist,function(x) rowSums(testset[,(x+1)])/length(x))
  # testy <- testset$days
  # # svm1(traindata,trainy,testdata,testy)
  # X <- traindata
  # y <- trainy
  # # svm.model <- svm(X,as.factor(y),probability=TRUE, gamma = 2^-7, cost = 2)
  # # svm.model <- svm(X,as.factor(y),probability=TRUE, kernel="linear")
  # svm.model <- svm(X,as.factor(y),probability=TRUE)
  # svm.pred <- predict(svm.model, testdata,probability=TRUE)
  # if(trainy[1]==TRUE){
  #   predvec <- attr(svm.pred,"prob")[,1] 
  # }else{
  #   predvec <- attr(svm.pred,"prob")[,2] 
  # }
  # pred <- prediction(predvec,as.numeric(testy))
  # roc_obj <- roc(as.numeric(testy), attr(svm.pred,"prob")[,1])
  # plot(roc_obj, print.auc=TRUE)
  # # auc <- performance(pred,"auc")
  # # auc <- roc.curve(scores.class0 = svm.pred, weights.class0 = y)$auc
  # aupr <- pr.curve(scores.class0 = predvec, weights.class0 = y)$auc.integral
  # # auc@y.values
  # c(roc_obj$auc,aupr)
}

splits <- h2o.splitFrame(train.hex, 0.75, seed=1234)
dl <- h2o.deeplearning(x=1:3, y="petal_len",training_frame=splits[[1]],
                         distribution="quantile", quantile_alpha=0.8)
h2o.predict(dl, splits[[2]])
