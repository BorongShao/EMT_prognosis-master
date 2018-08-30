lengthcpy <- length(cpy)
# t-test
ttest_sf <- sapply(1:lengthcpy, function(x) { 
  sf_74 <- rep(0,numfeatures)
  sf_74[cpy[[x]][[2]]] <- 1
  sf_74
})


lasso_sf <- sapply(1:lengthcpy, function(x) { 
  sf_74 <- rep(0,numfeatures)
  sf_74[cpy[[x]][[4]]] <- 1
  sf_74
})

netlasso_sf <- sapply(1:lengthcpy, function(x) { 
  sf_74 <- rep(0,numfeatures)
  sf_74[cpy[[x]][[6]]] <- 1
  sf_74
})

stsvm_sf <- sapply(1:lengthcpy, function(x) { 
  sf_74 <- rep(0,numfeatures)
  sf_74[cpy[[x]][[8]]] <- 1
  sf_74
})

# addg_sf <- sapply(1:lengthcpy, function(x) {
#   sf_74 <- rep(0,numfeatures)
#   temp <- unlist(cpy[[x]][[10]])
#   for(i in 1:length(temp))
#     sf_74[temp[i]] <- sf_74[temp[i]] + 1
#   sf_74
# })

addg_collec <- list()
for(i in 1:lengthcpy) {
  addg_collec <- c(addg_collec,cpy[[i]][[10]])
}
addg_collec_unique <- unique(addg_collec)
addg_sf1 <- rep(0,length(addg_collec_unique))
for(x in 1:lengthcpy){
  print(x)
  temp <- cpy[[x]][[10]]
  for(i in 1:length(temp)){
    for(j in 1:length(addg_collec_unique)){
      if(all.equal(addg_collec_unique[[j]], temp[[i]]) ==TRUE){
        addg_sf1[j] = addg_sf1[j] + 1
        break
      }
    }
  }
}
addg_sf <- addg_collec_unique[order(addg_sf1,decreasing = TRUE)[1:100]]

netrank_sf <- sapply(1:lengthcpy, function(x) { 
  sf_74 <- rep(0,numfeatures)
  sf_74[cpy[[x]][[12]]] <- 1
  sf_74
})

cox_sf <- sapply(1:lengthcpy, function(x) { 
  sf_74 <- rep(0,numfeatures)
  sf_74[cpy[[x]][[14]]] <- 1
  sf_74
})

coxlasso_sf <- sapply(1:lengthcpy, function(x) { 
  sf_74 <- rep(0,numfeatures)
  sf_74[cpy[[x]][[16]]] <- 1
  sf_74
})

rds_sf <- sapply(1:lengthcpy, function(x) { 
  sf_74 <- rep(0,numfeatures)
  sf_74[cpy[[x]][[18]]] <- 1
  sf_74
})

# survnet_sf <- sapply(1:lengthcpy, function(x) {
#   sf_74 <- rep(0,numfeatures)
#   temp <- unlist(cpy[[x]][[20]])
#   for(i in 1:length(temp))
#     sf_74[temp[i]] <- sf_74[temp[i]] + 1
#   sf_74
# })

survnet_collec <- list()
for(i in 1:lengthcpy) {
  survnet_collec <- c(survnet_collec,cpy[[i]][[20]])
}
survnet_collec_unique <- unique(survnet_collec)
survnet_sf1 <- rep(0,length(survnet_collec_unique))
for(x in 1:lengthcpy){
  print(x)
  temp <- cpy[[x]][[20]]
  for(i in 1:length(temp)){
    for(j in 1:length(survnet_collec_unique)){
      if(all.equal(survnet_collec_unique[[j]], temp[[i]]) ==TRUE){
        survnet_sf1[j] = survnet_sf1[j] + 1
        break
      }
    }
  }
}
survnet_sf <- survnet_collec_unique[order(survnet_sf1,decreasing = TRUE)[1:100]]

# 
# frequency <- data.frame(1:numfeatures, rowSums(ttest_sf), rowSums(lasso_sf), rowSums(netlasso_sf), rowSums(addg_sf), rowSums(netrank_sf),
#                         rowSums(stsvm_sf),rowSums(cox_sf),rowSums(coxlasso_sf),rowSums(rds_sf),rowSums(survnet_sf))
# names(frequency) <- c("index","t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet")
# frequency[,-1] <- sapply(frequency[,-1],function(x) x/max(x))
# 
# temp <- thresdata[,c(lassofs+1,ncol(thresdata))]
ttestfs <- order(rowSums(ttest_sf),decreasing = TRUE)[1:numtopfs]
#lasso
lassofs <- order(rowSums(lasso_sf),decreasing = TRUE)[1:numtopfs]
#netowrk lasso
netlassofs <- order(rowSums(netlasso_sf),decreasing = TRUE)[1:numtopfs]
#additive subnetworks
# addgfs <- order(rowSums(addg_sf),decreasing = TRUE)[1:numtopfs]
#NetRank
netrankfs <- order(rowSums(netrank_sf),decreasing = TRUE)[1:numtopfs]
#stSVM
stsvmfs <- order(rowSums(stsvm_sf),decreasing = TRUE)[1:numtopfs]
#Cox
coxfs <- order(rowSums(cox_sf),decreasing = TRUE)[1:numtopfs]
#Coxlasso
coxregfs <- order(rowSums(coxlasso_sf),decreasing = TRUE)[1:numtopfs]
#Coxlasso
rdsfs <- order(rowSums(rds_sf),decreasing = TRUE)[1:numtopfs]
#Coxlasso
# survnetfs <- order(rowSums(survnet_sf),decreasing = TRUE)[1:numtopfs]

temp <- sort(table(c(ttestfs,lassofs,netlassofs, unique(unlist(addg_sf)), netrankfs, stsvmfs, coxfs, coxregfs, rdsfs, 
                     unique(unlist(survnet_sf)))),decreasing = TRUE)[1:numtopfs]
alg10sfs <- as.numeric(names(temp))


 
# 
# net74_level2_freqindex <- order(temp,decreasing = TRUE)[1:num]
# net74_level1_freqindex <- order(temp,decreasing = TRUE)[1:num]
# net74_level3_freqindex <- order(temp,decreasing = TRUE)[1:num]
# temp74level1 <- thresdata74_1[,net74_level1_freqindex + 1]
# temp74level2 <- thresdata74_2[,net74_level2_freqindex + 1]
# temp74level3 <- thresdata74_3[,net74_level3_freqindex + 1]
# 
# 
# tsne_model_1 = Rtsne(as.matrix(temp455level3), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
# d_tsne_1 = as.data.frame(tsne_model_1$Y)
# ggplot(d_tsne_1, aes(x=V1, y=V2)) +  
#   geom_point(size=0.25) +
#   guides(colour=guide_legend(override.aes=list(size=6))) +
#   xlab("") + ylab("") +
#   ggtitle("t-SNE") +
#   theme_light(base_size=20) +
#   theme(axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   scale_colour_brewer(palette = "Set2")
# ## keeping original data
# d_tsne_1_original=d_tsne_1
# 
# ## Creating k-means clustering model, and assigning the result to the data used to create the tsne
# fit_cluster_kmeans=kmeans(scale(d_tsne_1), 3)  
# d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)
# 
# ## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
# fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))
# 
# ## setting 3 clusters as output
# d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=3))  
# 
# d_tsne_1_original$classlabel = factor(thresdata455_3$days)
# 
# plot_cluster=function(data, var_cluster, palette)  
# {
#   ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
#     geom_point(size=0.25) +
#     guides(colour=guide_legend(override.aes=list(size=6))) +
#     xlab("") + ylab("") +
#     ggtitle("") +
#     theme_light(base_size=10) +
#     theme(axis.text.x=element_blank(),
#           axis.text.y=element_blank(),
#           legend.direction = "horizontal", 
#           legend.position = "bottom",
#           legend.box = "horizontal") + 
#     scale_colour_brewer(palette = palette) 
# }
# 
# 
# plot_k3=plot_cluster(d_tsne_1_original, "cl_kmeans", "Set1")  
# plot_h3=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")
# plot_t3=plot_cluster(d_tsne_1_original, "classlabel", "Set2")
# 
# ## and finally: putting the plots side by side with gridExtra lib...
# 
# grid.arrange(plot_k1, plot_h1, plot_t1,plot_k2, plot_h2, plot_t2,plot_k3, plot_h3, plot_t3, ncol=3)  
# 
# 
# net123_level2_freqindex <- order(temp,decreasing = TRUE)[1:num]
# net123_level1_freqindex <- order(temp,decreasing = TRUE)[1:num]
# net123_level3_freqindex <- order(temp,decreasing = TRUE)[1:num]
# temp123level1 <- thresdata123_1[,net123_level1_freqindex + 1]
# temp123level2 <- thresdata123_2[,net123_level2_freqindex + 1]
# temp123level3 <- thresdata123_3[,net123_level3_freqindex + 1]
# 
# 
# #num=20
# net455_level2_freqindex <- order(temp,decreasing = TRUE)[1:num]
# net455_level1_freqindex <- order(temp,decreasing = TRUE)[1:num]
# net455_level3_freqindex <- order(temp,decreasing = TRUE)[1:num]
# temp455level1 <- thresdata455_1[,net455_level1_freqindex + 1]
# temp455level2 <- thresdata455_2[,net455_level2_freqindex + 1]
# temp455level3 <- thresdata455_3[,net455_level3_freqindex + 1]
# 

# plus clnical features
