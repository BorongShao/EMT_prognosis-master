library('Grace')
library('igraph')
library('cluster')
library('caret')
library('survivalROC')
# library('network')
library('infotheo')
library('glmnet')
# library('caTools')
library('e1071')
library('randomForest')
library('survival')
library('parallel')
library('ROCR')
library('PRROC')
library('foreach')
library('doMC')
# library('doSNOW')
library('ggplot2')
library('RColorBrewer')
# library('extrafont')
# library('survcomp')
library('pROC')
# library('superpc')
# library('MCL')
# library('splines')
# library('WGCNA')
# library('IsolationForest')
# library('h2o')
library('mice')
library('gridExtra') 
library('tidyr')
library('survminer')  
library('Rtsne')
library('SNFtool')
library('iCluster')
library('plyr')
library('reshape')

options(stringsAsFactors = FALSE)

# find common molecules 
setwd('/home/borong/Desktop/Experiments/LUAD')
SEnodes <- read.csv(file='LUAD_EMT2/molecules_all_name',header=FALSE)
SEnodes$id <- 1:75
SEnodes <- data.frame(SEnodes$id, SEnodes$V1)
names(SEnodes) <- c('id','name')
SEnodes$type <- c(rep(1,62),rep(2,13))
SEnodes$name <- as.character(SEnodes$name)

#pick out samples according to threshold
thresupper = 1400
threslower = 700
alivesamples <- read.csv(file='Data/clinical/alivesamples.txt',header=TRUE, sep=' ')
deadsamples <- read.csv(file='Data/clinical/deadsamples.txt',header=TRUE, sep=' ')
allsamples <- rbind(deadsamples,alivesamples)
alivesensor <- alivesamples[alivesamples$days >= thresupper,]
#alivesensor <- alivesamples
samples <- rbind(deadsamples,alivesensor)
samplesthres <- samples[samples$days >= thresupper | samples$days <threslower, ]



RNA_Seq_case_all <- read.csv(file='Data/RNA_Seq/all_case_rnaseq_515.txt',header=TRUE, sep='\t')
# RNA_Seq_case_all <- RNA_Seq_case_all[ , colSums(is.na(RNA_Seq_case_all)) == 0]
# RNA_Seq_case_all <- cbind(RNA_Seq_case_all[1],RNA_Seq_case_all[-1][ ,colSums(RNA_Seq_case_all[-1]) != 0])

methy_case_all <- read.csv(file='Data/methylation/all_case_mean_458_methylation.txt',header=TRUE, sep='\t')
# methy_case_all <- methy_case_all[ , colSums(is.na(methy_case_all)) == 0]

# cna_case_all <- read.csv(file='Data/cna/all_case_516_cna.txt',header=TRUE, sep='\t')
# cna_case_all <- cna_case_all[ , colSums(is.na(cna_case_all)) == 0]

luad_names <- intersect(names(RNA_Seq_case_all),names(methy_case_all))

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


setwd('../BRCA/')

RNA_Seq_case_all <- read.csv(file='Data/RNA_Seq/all_case_1093_rnaseq.txt',header=TRUE, sep='\t')
# RNA_Seq_case_all <- RNA_Seq_case_all[ , colSums(is.na(RNA_Seq_case_all)) == 0]
# RNA_Seq_case_all <- cbind(RNA_Seq_case_all[1],RNA_Seq_case_all[-1][ ,colSums(RNA_Seq_case_all[-1]) != 0])
names(RNA_Seq_case_all) <- sapply(names(RNA_Seq_case_all), function(x) unlist(strsplit(x,"[.]"))[1])

methy_case_all <- read.csv(file='Data/methylation/all_case_783_methylation.txt',header=TRUE, sep='\t')
# methy_case_all <- methy_case_all[ , colSums(is.na(methy_case_all)) == 0]

# cna_case_all <- read.csv(file='Data/cna/all_case_1080_cna.txt',header=TRUE, sep='\t')
# cna_case_all <- cna_case_all[ , colSums(is.na(cna_case_all)) == 0]

brca_names <- intersect(names(RNA_Seq_case_all),names(methy_case_all))

# setwd('../LUSC')
# 
# RNA_Seq_case_all <- read.csv(file='Data/RNA_Seq/all_case_501_rnaseq.txt',header=TRUE, sep='\t')
# RNA_Seq_case_all <- RNA_Seq_case_all[ , colSums(is.na(RNA_Seq_case_all)) == 0]
# RNA_Seq_case_all <- cbind(RNA_Seq_case_all[1],RNA_Seq_case_all[-1][ ,colSums(RNA_Seq_case_all[-1]) != 0])
# names(RNA_Seq_case_all) <- sapply(names(RNA_Seq_case_all), function(x) unlist(strsplit(x,"[.]"))[1])
# 
# methy_case_all <- read.csv(file='Data/methylation/all_case_370_methylation.txt',header=TRUE, sep='\t')
# methy_case_all <- methy_case_all[ , colSums(is.na(methy_case_all)) == 0]
# 
# cna_case_all <- read.csv(file='Data/cna/all_case_501_cna.txt',header=TRUE, sep='\t')
# cna_case_all <- cna_case_all[ , colSums(is.na(cna_case_all)) == 0]
# 
# lusc_names <- intersect(intersect(names(RNA_Seq_case_all),names(methy_case_all)), names(cna_case_all))
# 
# setwd('../LGG')
# RNA_Seq_case_all <- read.csv(file='Data/RNA_Seq/all_case_501_rnaseq.txt',header=TRUE, sep='\t')
# RNA_Seq_case_all <- RNA_Seq_case_all[ , colSums(is.na(RNA_Seq_case_all)) == 0]
# RNA_Seq_case_all <- cbind(RNA_Seq_case_all[1],RNA_Seq_case_all[-1][ ,colSums(RNA_Seq_case_all[-1]) != 0])
# names(RNA_Seq_case_all) <- sapply(names(RNA_Seq_case_all), function(x) unlist(strsplit(x,"[.]"))[1])
# 
# methy_case_all <- read.csv(file='Data/methylation/all_case_370_methylation.txt',header=TRUE, sep='\t')
# methy_case_all <- methy_case_all[ , colSums(is.na(methy_case_all)) == 0]
# 
# cna_case_all <- read.csv(file='Data/cna/all_case_501_cna.txt',header=TRUE, sep='\t')
# cna_case_all <- cna_case_all[ , colSums(is.na(cna_case_all)) == 0]
# 
# lgg_names <- intersect(intersect(names(RNA_Seq_case_all),names(methy_case_all)), names(cna_case_all))

#common names among three data types
# nameranges <- intersect(intersect(intersect(brca_names,lgg_names),luad_names),lusc_names)
nameranges <- intersect(brca_names,luad_names)
rm(RNA_Seq_case_all,methy_case_all)

setwd('../../Experiments_largeEMT_correct/')
temp1 <- read.csv('cancer genes/node_table_mirgenes.csv',stringsAsFactors = FALSE,header = TRUE)[,2]
temp2 <- read.csv('cancer genes/node_table_ppi.csv',stringsAsFactors = FALSE,header = TRUE)[,2]
temp3 <- read.csv('cancer genes/node_table_tf.csv',stringsAsFactors = FALSE,header = TRUE)[,2]
temp4 <- read.csv('emt.txt',stringsAsFactors = FALSE,header=FALSE)

addednodes <- unique(c(temp1,temp2,temp3,SEnodes$name))
allnodes <- sort(intersect(addednodes,nameranges))
# write.csv(file='allnodes.csv', allnodes, quote=FALSE, row.names=FALSE, col.names = FALSE)

#find the edges whose nodes are among 'allnodes'
ppiextd <- read.table('interactions/ppi.sif',sep='\t', header = FALSE,stringsAsFactors = FALSE)[,c(1,3)]
tfextd <- read.table('interactions/tf.sif',sep='\t', header = FALSE,stringsAsFactors = FALSE)[,c(1,3)]
mirextd <- read.table('interactions/mir.sif',sep='\t', header = FALSE,stringsAsFactors = FALSE)
ppiextd2 <- ppiextd[apply(ppiextd,1,function(x) (x[1] %in% allnodes & x[2] %in% allnodes)),]
tfextd2 <- tfextd[apply(tfextd,1,function(x) (x[1] %in% allnodes & x[2] %in% allnodes)),]
mirextd2 <- mirextd[apply(mirextd,1,function(x) (x[2] %in% allnodes)),]
emtedges <- read.csv(file='../Experiments/LUAD/LUAD_EMT2/interaction_revised.txt',sep='\t',header=FALSE,stringsAsFactors = FALSE)[,c(1,3)]
emtedges <-  emtedges[apply(emtedges,1,function(x) (x[1] %in% allnodes & x[2] %in% allnodes)),]
names(mirextd2) <- c("V1","V3")

LEedges <- unique(rbind(ppiextd2,tfextd2,mirextd2,emtedges))
# LEedges <- unique(emtedges)
names(LEedges) <- c('from','to')
largeemtnodes <- sort(unique(c(LEedges[,1],LEedges[,2])))
# largeemtnodes <- gsub("\\.","-",largeemtnodes)
largeemtid <- 1:length(largeemtnodes)
LEnodes <- data.frame(largeemtid,largeemtnodes)
names(LEnodes) <- c('id','name')
LEnodes$name <- as.character(LEnodes$name)
LEedges <- data.frame(match(LEedges[,1],LEnodes$name),match(LEedges[,2],LEnodes$name))
LEnet <- graph_from_data_frame(d=LEedges, vertices=LEnodes, directed=F)
LEnet <- simplify(LEnet, remove.multiple = TRUE)

#prediction without network
# rnaseq <- read.csv(file='../Experiments/LUAD/Data/RNA_Seq/case_rnaseq_515.txt',header=TRUE, sep='\t')
# rnaseq_case <- rnaseq[rnaseq[,1] %in% samplesthres[,1], ]
# mirna <- read.csv(file='../Experiments/LUAD/Data/miRNA/case_513_miRNA.txt',header=TRUE, sep='\t')
# names(mirna) <-  gsub(".",'-', names(mirna),fixed=TRUE)
# mirna_case <- mirna[mirna[,1] %in% samplesthres[,1], ]
# mrnaexp <- writedatasets(rnaseq_case,2:69)
toremove <- which(sapply(2:(ncol(LEthresgeneexp)-2), function(x) var(LEthresgeneexp[x])) < 1) +1
rmoutlier <- LEthresgeneexp[,-toremove]
allnodes <- names(rmoutlier)[-c(1,ncol(rmoutlier))]



# cnaindex <- match(nodes$name[1:68],names(cna_case_se))
# proteinindex <- match(nodes$name[1:68],names(protein_case_se))

#-----------------------------LUAD Data--------------------------------------------------------
LErnaseq_all <- read.csv(file='../Experiments/LUAD/Data/RNA_Seq/all_case_rnaseq_515.txt',header=TRUE, sep='\t')
# LErnaseqindex <- match(LEnodes$name[-c(216:231)],names(LErnaseq_all))
LErnaseqindex <- match(V(LEnet)$name[-c(33:45)],names(LErnaseq_all))
LEmirna_all <- read.csv(file='../Experiments/LUAD/Data/miRNA/all_case_513_miRNA.txt',header=TRUE, sep='\t')
names(LEmirna_all) <-  gsub(".",'-', names(LEmirna_all),fixed=TRUE)
names(LEmirna_all) <-  gsub("mir",'miR', names(LEmirna_all),fixed=TRUE)
LEmirnaindex <- match(V(LEnet)$name[33:45],names(LEmirna_all))
# LEmethy_all <- read.csv(file='../Experiments/LUAD/Data/methylation/all_case_mean_458_methylation.txt',header=TRUE, sep='\t')
# names(LEmethy_all) <-  gsub(".",'-', names(LEmethy_all),fixed=TRUE)
methyindex <- match(V(LEnet)$name,names(LEmethy_all))
LEcna_all <- read.csv(file='../Experiments/LUAD/Data/cna/all_case_516_cna.txt',header=TRUE, sep='\t')
#remove the features with more than 50 missing values
tokeep <- sapply(LEcne_all, function(x) sum(is.na(x))<100)
LEcna_all <- LEcna_all[,tokeep]
# rename mir features so that it can be matched.
tomodify <- names(LEcna_all)[grep('MIR',names(LEcna_all))]
tomodify <- sapply(tomodify,tolower)
tomodify <- gsub("mir",'miR-', tomodify,fixed=TRUE)
names(LEcna_all)[grep('MIR',names(LEcna_all))] <- tomodify
# LEcna_all$TRIB3 <- NULL
LEcnaindex <- match(V(LEnet455)$name,names(LEcna_all))
# LEnet74_cna <- subgraph(LEnet74,which(!is.na(LEcnaindex)))
LEnet455_cna <- subgraph(LEnet455,which(!is.na(LEcnaindex)))
LEcnaindex2 <- match(V(LEnet455_cna)$name,names(LEcna_all))
# LEprotein_all <- read.csv(file='../Experiments/LUAD/Data/protein/all_case_protein_365.txt',header=TRUE, sep='\t')
#-------------------------LUAD data with threshold -----------------------------------------

comms <- intersect(LErnaseq_all[,1],LEmirna_all[,1])
comms2 <- comms[comms %in% samplesthres[,1]]
LErnaseq_case <- LErnaseq_all[match(comms2, LErnaseq_all[,1]),c(1,LErnaseqindex)]
LEmirna_case <- LEmirna_all[match(comms2, LEmirna_all[,1]),c(1,LEmirnaindex)]
LEgeneexp <- data.frame(LErnaseq_case,LEmirna_case[,-1])
names(LEgeneexp) <-  gsub(".",'-', names(LEgeneexp),fixed=TRUE)
#gene expression data
X <- log2(LEgeneexp[,c(match(V(LEnet)$name, names(LEgeneexp)))]+1)
days <- samplesthres[match(LEgeneexp[,1], samplesthres[,1]),2] >= thresupper
LEthresgeneexp <- data.frame(LEgeneexp[,1],X,days)
names(LEthresgeneexp) <-  gsub(".",'-', names(LEthresgeneexp),fixed=TRUE)

X <- LEthresgeneexp[,-c(1,ncol(LEthresgeneexp))]
tr<-IsolationTrees(X,ntree = 100, rFactor = 0, rowSamp = T, nRowSamp = 10)
#evaluate anomaly score
as<-AnomalyScore(X,tr)
# show anomaly score
hist(as$outF,20)
# # 8 samples removed
rmoutlier <- LEthresgeneexp[-(which(as$outF > 0.6)),]

# no_outlier <- read.csv('samples.csv')
#remove outlier
# X <- t(LEthresgeneexp[-c(77,136),-c(1,ncol(LEthresgeneexp))])
X <- LEthresgeneexp[,-c(1,ncol(LEthresgeneexp))]
# X <- rmoutlier[,-c(1,ncol(rmoutlier))]
IAC = dist(X)
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
cluster1=hclust(IAC,method="average")
# cutree(cluster1,h=0.2)
plot(cluster1,cex=0.7)
# rmoutlier <- LEthresgeneexp[-which(cutree(cluster1,h=36)!=1),]
# LEthresgeneexp <- LEthresgeneexp[match(unlist(no_outlier), LEthresgeneexp[,1]),]


comms <- LEmethy_all[,1][LEmethy_all[,1] %in% samplesthres[,1]]
LEmethy_case <- LEmethy_all[match(comms, LEmethy_all[,1]),c(1,methyindex)]
names(LEmethy_case) <-  gsub(".",'-', names(LEmethy_case),fixed=TRUE)
# methylation data
X <- scale(LEmethy_case[,c(match(V(LEnet)$name, names(LEmethy_case)))])
days <- samplesthres[match(LEmethy_case[,1], samplesthres[,1]),2] >= thresupper
LEthresmethy <- data.frame(LEmethy_case[,1],X,days)
names(LEthresmethy) <-  gsub(".",'-', names(LEthresmethy),fixed=TRUE)

#remove outlier
X <- t(LEthresmethy[,-c(1,ncol(LEthresmethy))])
IAC = cor(X, use="p")
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
cluster1=hclust(as.dist(1-IAC),method="average")
plot(cluster1,cex=0.7)
    

comms <- LEcna_all[,1][LEcna_all[,1] %in% samplesthres[,1]]
LEcna_case <- LEcna_all[match(comms, LEcna_all[,1]),c(1,LEcnaindex2)]
# methylation data
X <- scale(LEcna_case[,c(match(V(LEnet455_cna)$name, names(LEcna_case)))])
days <- samplesthres[match(LEcna_case[,1], samplesthres[,1]),2] >= thresupper
LEthrescna <- data.frame(LEcna_case[,1],X,days)
names(LEthrescna) <-  gsub(".",'-', names(LEthrescna),fixed=TRUE)
# tempData <- mice(LEthrescna,m=1)
# LEthrescna <- complete(tempData,1)
#impute missing value with mean
toimpute <- which(sapply(LEthrescna, function(x) sum(is.na(x))) > 0)
for(i in 1:length(toimpute)){
  temp <- LEthrescna[,toimpute[i]]
  LEthrescna[,toimpute[i]] <- Hmisc::impute(temp,mean)
}

write.csv(LEthrescna, file='../Thesis_Results/455network/emt445_cna_700_1400.csv', row.names = FALSE)


#--------------------LUAD data without threshold-----------------------------------
comms <- intersect(LErnaseq_all[,1],LEmirna_all[,1])
comms2 <- comms[comms %in% allsamples[,1]]
LErnaseq_case <- LErnaseq_all[match(comms2, LErnaseq_all[,1]),c(1,LErnaseqindex)]
LEmirna_case <- LEmirna_all[match(comms2, LEmirna_all[,1]),c(1,LEmirnaindex)]
LEgeneexp <- data.frame(LErnaseq_case,LEmirna_case[,-1])
names(LEgeneexp) <-  gsub(".",'-', names(LEgeneexp),fixed=TRUE)
#gene expression data
X <- log(LEgeneexp[,c(match(V(LEnet)$name, names(LEgeneexp)))]+1)
days <- allsamples[match(LEgeneexp[,1], allsamples[,1]),2] #>= thresupper
LEgeneexp <- data.frame(LEgeneexp[,1],X,days)
names(LEgeneexp) <-  gsub(".",'-', names(LEgeneexp),fixed=TRUE)

#remove outliers
IAC = cor(X, use="p")
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
cluster1=hclust(as.dist(1-IAC),method="average")
plot(cluster1,cex=0.7)

comms <- LEmethy_all[,1][LEmethy_all[,1] %in% allsamples[,1]]
LEmethy_case <- LEmethy_all[match(comms, LEmethy_all[,1]),c(1,methyindex)]
names(LEmethy_case) <-  gsub(".",'-', names(LEmethy_case),fixed=TRUE)
# methylation data
X <- scale(LEmethy_case[,c(match(V(LEnet)$name, names(LEmethy_case)))])
days <- allsamples[match(LEmethy_case[,1], allsamples[,1]),2] #>= thresupper
LEmethy <- data.frame(LEmethy_case[,1],X,days)
names(LEmethy) <-  gsub(".",'-', names(LEmethy),fixed=TRUE)

comms <- LEcna_all[,1][LEcna_all[,1] %in% allsamples[,1]]
LEcna_case <- LEcna_all[match(comms, LEcna_all[,1]),c(1,LEcnaindex2)]
X <- scale(LEcna_case[,c(match(V(LEnet455_cna)$name, names(LEcna_case)))])
days <- allsamples[match(LEcna_case[,1], allsamples[,1]),2] 
LEcna <- data.frame(LEcna_case[,1],X,days)
names(LEcna) <-  gsub(".",'-', names(LEcna),fixed=TRUE)
#impute missing value with mean
toimpute <- which(sapply(LEcna, function(x) sum(is.na(x))) > 0)
for(i in 1:length(toimpute)){
  temp <- LEcna[,toimpute[i]]
  LEcna[,toimpute[i]] <- Hmisc::impute(temp,mean)
}

write.csv(LEcna, file='../Thesis_Results/455network/emt445_cna_700_1400_nothres.csv', row.names = FALSE)


#remove outlier
X <- LEmethy[,-c(1,ncol(LEmethy))]
IAC = cor(X, use="p")
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
cluster1=hclust(as.dist(1-IAC),method="average")
plot(cluster1,cex=0.7)

# rm(LErnaseq_all,LEmirna_all,LEmethy_all)
temp1 <- c(names(LErnaseq_all)[-1])
temp2 <- names(LEmethy_all[,-1])
temp <- intersect(temp1,temp2)
LEfinodes <- NULL
LEfiedges <- read.csv('../PPI_network/FIsInGene_Pathway_041709.txt', sep='\t')
LEfiedges <- LEfiedges[(LEfiedges[,1] %in% temp & LEfiedges[,2] %in% temp),]
LEfinodes$id <- 1:length(unique(c(LEfiedges[,1],LEfiedges[,2])))
LEfinodes$name <- unique(c(LEfiedges[,1],LEfiedges[,2]))
LEfiedges <- data.frame(match(LEfiedges[,1],LEfinodes$name),match(LEfiedges[,2],LEfinodes$name))
LEfinet <- graph_from_data_frame(d=LEfiedges, vertices=LEfinodes, directed=F)
LEfinet <- simplify(LEfinet, remove.multiple = TRUE)

LEfirnaseqindex <- match(V(LEfinet)$name,names(LErnaseq_all))
LEfimethyindex <- match(V(LEfinet)$name,names(LEmethy_all))

#-------------------------LUAD data with threshold -----------------------------------------

comms <- LErnaseq_all[,1]
comms2 <- comms[comms %in% samplesthres[,1]]
LEfigeneexp <-LErnaseq_all[match(comms2, LErnaseq_all[,1]),c(1,LEfirnaseqindex)]
X <- log2(LEfigeneexp[,c(match(V(LEfinet)$name, names(LEfigeneexp)))]+1)
days <- samplesthres[match(LEfigeneexp[,1], samplesthres[,1]),2] >= thresupper
LEfithresgeneexp <- data.frame(LEfigeneexp[,1],X,days)
names(LEfithresgeneexp) <-  gsub(".",'-', names(LEfithresgeneexp),fixed=TRUE)





#---------------------------------#normalize data with normal samples---------------------------------------------------------------#

LErnaseq_allnormal <- read.csv(file='../Experiments/LUAD/Data/RNA_Seq/all_normal_rnaseq_59.txt',header=TRUE, sep='\t')
# LErnaseqindex <- match(LEnodes$name[-c(216:231)],names(LErnaseq_all))
LErnaseqindex <- match(V(LEnet)$name[-c(33:45)],names(LErnaseq_allnormal))
LEmirna_allnormal <- read.csv(file='../Experiments/LUAD/Data/miRNA/all_normal_46_miRNA.txt',header=TRUE, sep='\t')
names(LEmirna_allnormal) <-  gsub(".",'-', names(LEmirna_allnormal),fixed=TRUE)
names(LEmirna_allnormal) <-  gsub("mir",'miR', names(LEmirna_allnormal),fixed=TRUE)
LEmirnaindex <- match(V(LEnet)$name[33:45],names(LEmirna_allnormal))
LEmethy_allnormal <- read.csv(file='../Experiments/LUAD/Data/methylation/all_normal_32_methylation.txt',header=TRUE, sep='\t')
# names(LEmethy_allnormal) <-  gsub(".",'-', names(LEmethy_allnormal),fixed=TRUE)
methyindex <- match(V(LEnet)$name,names(LEmethy_allnormal))
LEcna_allnormal <- read.csv(file='../Experiments/LUAD/Data/cna/all_normal_577_cna.txt',header=TRUE, sep='\t')
LEcnaindex <- match(V(LEnet)$name[-c(33:45)],names(LEcna_allnormal))
# LEprotein_all <- read.csv(file='../Experiments/LUAD/Data/protein/all_n',header=TRUE, sep='\t')
#-------------------------LUAD data with threshold -----------------------------------------

rnaseq_normal_mean <- sapply(LErnaseq_allnormal[,LErnaseqindex], mean)
rnaseq_normal_sd <- sapply(LErnaseq_allnormal[,LErnaseqindex], sd)
mirna_normal_mean <- sapply(LEmirna_allnormal[,LEmirnaindex],mean)
mirna_normal_sd <- sapply(LEmirna_allnormal[,LEmirnaindex],sd)
level1_mean_vec <- c(rnaseq_normal_mean,mirna_normal_mean)
level1_sd_vec <- c(rnaseq_normal_sd,mirna_normal_sd)
level1_mean <- level1_mean_vec[c(match(V(LEnet)$name, names(level1_mean_vec)))]
level1_sd <- level1_sd_vec[c(match(V(LEnet)$name, names(level1_sd_vec)))]
level1_sd_zeros <- which(level1_sd==0)
level1_sd[level1_sd_zeros] <- 1

X <- LEgeneexp[,c(match(V(LEnet)$name, names(LEgeneexp)))]
X_normal <- sapply(1:ncol(X), function(x) (X[,x]-level1_mean[x])/level1_sd[x])
X_normal <- as.data.frame(X_normal)
names(X_normal) <- names(X)
days <- samplesthres[match(LEgeneexp[,1], samplesthres[,1]),2] >= thresupper
LEthresgeneexp <- data.frame(LEgeneexp[,1],X_normal,days)
names(LEthresgeneexp) <-  gsub(".",'-', names(LEthresgeneexp),fixed=TRUE)

rmoutlier <- LEthresgeneexp


# genes targeted by miRNA
temp <- which(grepl("miR",get.edgelist(LEnet455)[,1]) == TRUE)
mirtargetgenes <- get.edgelist(LEnet455)[temp,2]
#nothreslevel1, 2, 3
nothreslevel1surv <- Surv(nothreslevel1$days, nothreslevel1[,1] %in% deadsamples[,1]) 
nothreslevel2surv <- Surv(nothreslevel2$days, nothreslevel2[,1] %in% deadsamples[,1]) 
nothreslevel3surv <- Surv(nothreslevel3$days, nothreslevel3[,1] %in% deadsamples[,1]) 

level1pvalue455 <-  sapply(nothreslevel1[,-1], function(x) coef(summary(coxph(nothreslevel1surv ~ as.numeric(x))))[5])   #104,54
level2pvalue455 <-  sapply(nothreslevel2[,-1], function(x) coef(summary(coxph(nothreslevel2surv ~ as.numeric(x))))[5])   # 53,13
level3pvalue455 <-  sapply(nothreslevel3[,-1], function(x) coef(summary(coxph(nothreslevel3surv ~ as.numeric(x))))[5])   #36,9
level1pvalue123 <-  sapply(nothreslevel1[,-1], function(x) coef(summary(coxph(nothreslevel1surv ~ as.numeric(x))))[5])   #37,21
level2pvalue123 <-  sapply(nothreslevel2[,-1], function(x) coef(summary(coxph(nothreslevel2surv ~ as.numeric(x))))[5])   # 15,2
level3pvalue123 <-  sapply(nothreslevel3[,-1], function(x) coef(summary(coxph(nothreslevel3surv ~ as.numeric(x))))[5])   #7,3
level1pvalue74 <-  sapply(nothreslevel1[,-1], function(x) coef(summary(coxph(nothreslevel1surv ~ as.numeric(x))))[5])   #15,9
level2pvalue74 <-  sapply(nothreslevel2[,-1], function(x) coef(summary(coxph(nothreslevel2surv ~ as.numeric(x))))[5])   # 13,3
level3pvalue74 <-  sapply(nothreslevel3[,-1], function(x) coef(summary(coxph(nothreslevel3surv ~ as.numeric(x))))[5])   #9,3

mirtargetlevel1pvalue <-  sapply(nothreslevel1[,match(mirtargetgenes,names(nothreslevel1))], function(x) coef(summary(coxph(nothreslevel1surv ~ as.numeric(x))))[5])   #40,18
mirtargetlevel2pvalue <-  sapply(nothreslevel2[,match(mirtargetgenes,names(nothreslevel2))], function(x) coef(summary(coxph(nothreslevel2surv ~ as.numeric(x))))[5])   # 29,5
mirtargetgenes <- get.edgelist(LEnet455_cna)[temp,2]
mirtargetlevel3pvalue <-  sapply(nothreslevel3[,match(mirtargetgenes,names(nothreslevel3))], function(x) coef(summary(coxph(nothreslevel3surv ~ as.numeric(x))))[5])   #21,6

# chi-square test

