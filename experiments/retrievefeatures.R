
topfeatures <- 20

ttestfs <- cpyfeatures[[2]][[1]][1:topfeatures]
#lasso
lassofs <- cpyfeatures[[2]][[2]][1:topfeatures]
#netowrk lasso
netlassofs <- cpyfeatures[[2]][[3]][1:topfeatures]
#additive subnetworks
addg_sf <- cpyfeatures[[2]][[4]][1:topfeatures]
#NetRank
netrankfs <- cpyfeatures[[2]][[5]][1:topfeatures]
#stSVM
stsvmfs <-cpyfeatures[[2]][[6]][1:topfeatures]
#Cox
coxfs <- cpyfeatures[[2]][[7]][1:topfeatures]
#Coxlasso
coxregfs <- cpyfeatures[[2]][[8]][1:topfeatures]
#Coxlasso
rdsfs <-cpyfeatures[[2]][[9]][1:topfeatures]
#Coxlasso
survnet_sf <- cpyfeatures[[2]][[10]][1:topfeatures]

temp <- sort(table(c(ttestfs,lassofs,netlassofs, unique(unlist(addg_sf)), netrankfs, stsvmfs, coxfs, coxregfs, rdsfs, 
                     unique(unlist(survnet_sf)))),decreasing = TRUE)[1:numtopfs]
alg10sfs <- as.numeric(names(temp))
nofs <- 1:(ncol(nothresdata)-2)

# features_matlab <- function(nothresdata, network, datatype, index, clusters) {
#   if(index==4 | index==10){
#     addgfs <- addg_sf[1:numtopnets]
#     moldata <- nothresdata[,-c(1,ncol(nothresdata))]
#     directionadd <- sign(cor(moldata,nothresdata$days))
#     newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[x])
#     newX2 <- sapply(addgfs,function(x) rowSums(newX[,x])/length(x))
#     nothresdatanet <- data.frame(nothresdata[,1],newX2,nothresdata$days)
#     names(nothresdatanet)[ncol(nothresdatanet)] <- "days"
#     data <- nothresdatanet[1:ncol(newX2)+1]
#   } else{
#     data <- nothresdata[,get(algorithmstrings[i])+1]
#   }
#   stringfile <- paste('../matlabfiles/', network, '_', datatype,'_',index,".csv")
#   write.table(data.matrix(data), file=stringfile, row.names = FALSE,col.names = FALSE)
# }
# 
# for(i in 1:12){
#   features_matlab(nothresdata, 455, 'level3', i,clusters)
# }

