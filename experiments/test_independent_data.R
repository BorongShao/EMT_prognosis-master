
#network genes
molecules74 <- read.csv('net74.csv',header=FALSE,stringsAsFactors = FALSE)$V1
molecules123 <- read.csv('net123.csv',header=FALSE,stringsAsFactors = FALSE)$V1
molecules455 <- read.csv('net455.csv',header=FALSE,stringsAsFactors = FALSE)$V1
#clinical variables
norclinical <- read.csv('test/clinical.csv',sep=';',header = TRUE)
norclinical$egfr <- factor(norclinical$egfr)
norclinical$egfr_type <- factor(norclinical$egfr_type)
norclinical$kras <- factor(norclinical$kras)
norclinical$braf <- factor(norclinical$braf)
norclinical$p53 <- factor(norclinical$p53)

# methylation data
normethy <- read.table('methylation/methylation.txt',header = TRUE,stringsAsFactors = FALSE)
methy74 <- read.table('methylation/probes74',header=TRUE,quote = FALSE,stringsAsFactors = FALSE)
methy123 <- read.table('methylation/probes123',header=TRUE, quote = "",stringsAsFactors = FALSE)
methy455 <- read.table('methylation/probes455',header=TRUE,quote = "",stringsAsFactors = FALSE)
names(methy74) = names(methy123)  =names(methy455) = c("probe","genes","regions","chromosome","cpglocate")
methy_sampletitle <- read.table('methylation/sampletitle_methy',header=TRUE,stringsAsFactors = FALSE)
methyprobetable <- getprobetablemethy(methy455)

methyprobetable <- methyprobetable[order(methyprobetable$gene),]   # sort by gene
methyprobetable <- methyprobetable[which(!is.na(match(methyprobetable$gene,molecules455))),]     # remove genes not in the network

tumoridxmethy <- sapply(names(methy_sampletitle), function(x) unlist(strsplit(x,"_"))[2])
tumoridxmethy <- which(tumoridxmethy=="Tumor")
tumorsamplesmethy <- sapply(names(tumoridxmethy), function(x) unlist(strsplit(x,"_"))[1])
normethytumor <- normethy[,c(1,tumoridxmethy)]
normethynormal <- normethy[,-tumoridxmethy]
names(normethytumor) <- c("ID_REF",tumorsamplesmethy)


# univariate cox model, methylation probe level
methyclinical <- norclinical[match(tumorsamplesmethy,norclinical[,1]),]
methysurv <- Surv(methyclinical$time.pfs.months, methyclinical$event.pfs) 
methyprobepvalues <-  apply(normethytumor[,-1], 1, function(x) coef(summary(coxph(methysurv ~ as.numeric(x))))[5])
probegenes05 <- unique(methyprobetable[match(normethytumor[which(methyprobepvalues < 0.05),1],methyprobetable$probe),2]) 
probegenes05 <- probegenes05[-which(is.na(probegenes05))]                  # 342 genes have at least one probe with p < 0.05
probegenes01 <- unique(methyprobetable[match(normethytumor[which(methyprobepvalues < 0.01),1],methyprobetable$probe),2])  
probegenes01 <- probegenes01[-which(is.na(probegenes01))]                  # 171 genes have probes with p < 0.01
sum(!is.na(match(molecules455,probegenes05)))


#methylation gene level, average probes per gene
temp <- as.numeric(as.factor(methyprobetable$gene))
probecolumn <- normethytumor[,1]
methy_avg <- matrix(0,nrow=455,ncol=164)
idx <- unique(temp)
for(i in 1:length(idx)){
  probes <- methyprobetable$probe[temp==idx[i]]
  methy_avg[i,] <- colMeans(normethytumor[match(probes,normethytumor$ID_REF),-1])
}
methy_avg_pvalues <-  apply(methy_avg, 1, function(x) coef(summary(coxph(methysurv ~ as.numeric(x))))[5])  # 89 < 0.05, 33 < 0.01
methydata455 <- as.data.frame(cbind(names(normethytumor)[-1],t(methy_avg)))
data455names <- unique(as.character(as.factor(methyprobetable$gene)))
names(methydata455) <- c("sampleID",data455names)
names(methydata455)[233:245] <- V(LEnet455)$name[232:244]
pfsmonth <- norclinical$time.pfs.months[match(methydata455$sampleID,norclinical$GEO_ID)]
pfsevent <- norclinical$event.pfs[match(methydata455$sampleID,norclinical$GEO_ID)]
methydata455 <- methydata455[,c(1,match(V(LEnet455)$name,names(methydata455)))]
methydata123 <- methydata455[,c(1,match(V(LEnet123)$name,names(methydata455)))]
methydata74 <- methydata455[,c(1,match(V(LEnet74)$name,names(methydata455)))]
methydata74[,-1] <- sapply(methydata74[,-1], function(x) as.numeric(as.character(x)))
methydata123[,-1] <- sapply(methydata123[,-1], function(x) as.numeric(as.character(x)))
methydata455[,-1] <- sapply(methydata455[,-1], function(x) as.numeric(as.character(x)))



#methylation gene level, promoter region per gene
promoterprobe <- which(methyprobetable$region=="TSS1500" | methyprobetable$region=="TSS200")
methyprobepromoter <- methyprobetable[promoterprobe,]
temp <- as.numeric(as.factor(methyprobepromoter$gene))
idx <- unique(temp)
methy_promoter <- matrix(0,nrow=455,ncol=164)
for(i in 1:length(idx)){
  probes <- methyprobepromoter$probe[temp==idx[i]]
  methy_promoter[i,] <- colMeans(normethytumor[match(probes,normethytumor$ID_REF),-1])
}
# CSN2 does not have probes
methy_promoter_pvalues <-  apply(methy_promoter, 1, function(x) coef(summary(coxph(methysurv ~ as.numeric(x))))[5])[-454] #43 < 0.05, 8 < 0.01


#methylation gene level, least correlated probe per gene, using gene_median
temp <- as.numeric(as.factor(methyprobetable$gene))
genenames <- unique(as.character(as.factor(methyprobetable$gene)))
probecolumn <- normethytumor[,1]
methy_cor <- matrix(0,nrow=455,ncol=121)
matching <- match(names(normethytumor)[-1],names(gene_median)[-1])
matched <- matching[!is.na(matching)]
matchingsamples <- normethytumor[,c(1,matched+1)]
idx <- unique(temp)
for(i in 1:length(idx)){
  probes <- methyprobetable$probe[temp==idx[i]]
  genename <- genenames[i]
  expression <- gene_median[match(genename,gene_median$geneID),-1]
  if(!is.na(expression[1])){
    y <- matchingsamples[match(probes,normethytumor$ID_REF),-1]
    datamatridx <- matrix(0,nrow=nrow(y),ncol=ncol(y))
    for(j in 1:ncol(y)){
      datamatridx[,j] <- y[,j]
    }
    correlations <- apply(datamatridx, 1, function(x) {
      cor(x,as.numeric(as.character(unlist(expression))))
    })
    probeid <- which.min(correlations)
    methy_cor[i,] <- datamatridx[probeid,]
  }else{
    methy_cor[i,] <- colMeans(matchingsamples[match(probes,normethytumor$ID_REF),-1])
  }
}
methyclinical2 <- norclinical[match(names(matchingsamples)[-1],norclinical[,1]),]
methysurv2 <- Surv(methyclinical2$time.pfs.months, methyclinical2$event.pfs) 
methy_cor_pvalues <-  apply(methy_cor, 1, function(x) coef(summary(coxph(methysurv2 ~ as.numeric(x))))[5])  # 34 < 0.05, 6 < 0.01


# find the best for testing biomarkers. 



#microarray data
norgeneexp <- read.table('microarray/microarray.txt',header = TRUE,stringsAsFactors = FALSE)
genes74 <- read.table('microarray/74genes',header=FALSE, sep='\t',quote = "",stringsAsFactors = FALSE)
genes123 <- read.table('microarray/123genes',header=FALSE, sep='\t',quote = "",stringsAsFactors = FALSE)
genes455 <- read.table('microarray/455genes',header=FALSE, sep='\t',quote = "",stringsAsFactors = FALSE)
names(genes74) = names(genes123)  =names(genes455) = c("probe","gene","description")
geneprobetable <- genes455[order(genes455$gene),]
tumorsamplesgene <- match(names(norgeneexp),names(normethytumor))
sampleID <- read.table('microarray/microarray_sampleID',header = FALSE,stringsAsFactors = FALSE)
names(sampleID) <- c("GSM","tumor")
names(norgeneexp) <-  c("ID_REF",sapply(sampleID$tumor, function(x) paste("Sample", unlist(strsplit(x,"_"))[2],sep="")))

geneclinical <- norclinical[match(names(norgeneexp)[-1],norclinical[,1]),]
genesurv <- Surv(geneclinical$time.pfs.months, geneclinical$event.pfs) 

# univariate cox model, gene probe level
geneprobepvalues <-  apply(norgeneexp[,-1], 1, function(x) coef(summary(coxph(genesurv ~ as.numeric(x))))[5])     # 72 < 0.05, 20 < 0.01
variance <- apply(norgeneexp[,-1], 1, function(x) var(x))
avgexpression <- apply(norgeneexp[,-1], 1, function(x) mean(x))
temp <- data.frame(norgeneexp[,1],genes455$gene[match(norgeneexp[,1],genes455$probe)],geneprobepvalues,variance,avgexpression)
names(temp) <- c("probe","gene","pvalue","variance","average")
geneselection <- temp[order(temp$gene),]
#keep the genes that is in network
geneselection <- geneselection[which(!is.na(match(geneselection$gene,molecules455))),]


# take the mean for duplicated probes
temp <- as.numeric(as.factor(geneselection$gene))
probecolumn <- norgeneexp[,1]
gene_mean <- matrix(0,nrow=442,ncol=121)
idx <- unique(temp)
for(i in 1:length(idx)){
  probes <- geneselection$probe[temp==idx[i]]
  gene_mean[i,] <- colMeans(norgeneexp[match(probes,norgeneexp$ID_REF),-1])
}
gene_mean_pvalues <-  apply(gene_mean, 1, function(x) coef(summary(coxph(genesurv ~ as.numeric(x))))[5])     # 53 < 0.05, 14 < 0.01


# take the median for duplicated probes
temp <- as.numeric(as.factor(geneselection$gene))
genenames <-  unique(as.character(as.factor(geneselection$gene)))
probecolumn <- norgeneexp[,1]
gene_median <- matrix(0,nrow=442,ncol=121)
idx <- unique(temp)
for(i in 1:length(idx)){
  probes <- geneselection$probe[temp==idx[i]]
  gene_median[i,] <- apply((norgeneexp[match(probes,norgeneexp$ID_REF),-1]),2,function(x) median(x))
}
gene_median <- as.data.frame(cbind(genenames,gene_median))
names(gene_median) <- names(norgeneexp)
names(gene_median)[1] <- "geneID"
gene_median_pvalues <-  apply(gene_median[,-1], 1, function(x) coef(summary(coxph(genesurv ~ as.numeric(x))))[5])     # 51 < 0.05, 16 < 0.01


genedata455 <- as.data.frame(cbind(names(gene_median)[-1],t(gene_median[,-1])))
mrnazeros <- as.data.frame(matrix(0,nrow=nrow(genedata455),ncol=13))
data455names <- unique(as.character(as.factor(geneselection$gene)))
genedata455 <- as.data.frame(cbind(genedata455,mrnazeros))
names(genedata455) <- c("sampleID",data455names,V(LEnet455)$name[232:244])
pfsmonth <- norclinical$time.pfs.months[match(genedata455$sampleID,norclinical$GEO_ID)]
pfsevent <- norclinical$event.pfs[match(genedata455$sampleID,norclinical$GEO_ID)]
genedata455 <- genedata455[,c(1,match(V(LEnet455)$name,names(genedata455)))]
genedata123 <- genedata455[,c(1,match(V(LEnet123)$name,names(genedata455)))]
genedata74 <- genedata455[,c(1,match(V(LEnet74)$name,names(genedata455)))]
genedata74[,-1] <- sapply(genedata74[,-1], function(x) as.numeric(as.character(x)))
genedata123[,-1] <- sapply(genedata123[,-1], function(x) as.numeric(as.character(x)))
genedata455[,-1] <- sapply(genedata455[,-1], function(x) as.numeric(as.character(x)))






getprobetablemethy <- function(ptable) {
  results <-  data.frame(matrix(ncol = 3),stringsAsFactors = FALSE)
  names(results) <- c("probe","gene","region")
  numprobe <- nrow(ptable)
  count <- 0
  for(i in 1:numprobe){
    #print(i)
    probe <- as.character(ptable[i,1])
    genes <-unlist(strsplit(as.character(ptable[i,2]),";"))
    regions <- unlist(strsplit(as.character(ptable[i,3]),";"))
    for(j in 1:length(genes)){
      count <- count+1
      results[count,] <- c(probe,genes[j],regions[j])
    }
  }
  unique(results)
}

getprobetablegene <- function(ptable) {
  results <-  data.frame(matrix(ncol = 3),stringsAsFactors = FALSE)
  names(results) <- c("probe","gene","description")
  numprobe <- nrow(ptable)
  count <- 0
  for(i in 1:numprobe){
    #print(i)
    probe <- as.character(ptable[i,1])
    gene <-as.character(ptable[i,2])
    description <- as.character(ptable[i,3])
    count <- count+1
    results[count,] <- c(probe,gene,description)
    
  }
  unique(results)
}

# list objects by size
# sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = TRUE) [1:30]

# test FSF of 10 algoithms on 3 networks
# single data levels
# gene expression data
cpyfeatures <- LUAD74_geneexp_cluster3
mirrange <- 33:45
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures2.R')
LUAD74_testgene_cluster2 <- cluster_measure_test(genedata74,2,add_dire74[[1]],1,genesurv)
LUAD74_testgene_cluster3 <- cluster_measure_test(genedata74,3,add_dire74[[1]],1,genesurv)

cpyfeatures <- LUAD123_geneexp_cluster3
mirrange <- 70:82
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures2.R')
LUAD123_testgene_cluster2 <- cluster_measure_test(genedata123,2,add_dire123[[1]],1,genesurv)
LUAD123_testgene_cluster3 <- cluster_measure_test(genedata123,3,add_dire123[[1]],1,genesurv)

cpyfeatures <- LUAD455_geneexp_cluster3
mirrange <- 232:244
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures2.R')
LUAD455_testgene_cluster2 <- cluster_measure_test(genedata455,2,add_dire455[[1]],1,genesurv)
LUAD455_testgene_cluster3 <- cluster_measure_test(genedata455,3,add_dire455[[1]],1,genesurv)

LUAD_testgene_cluster2 <- rbind(sapply(LUAD74_testgene_cluster2[[1]], function(x) x[[4]]),
sapply(LUAD123_testgene_cluster2[[1]], function(x) x[[4]]),
sapply(LUAD455_testgene_cluster2[[1]], function(x) x[[4]]))

LUAD_testgene_cluster3 <- rbind(sapply(LUAD74_testgene_cluster3[[1]], function(x) x[[5]]),
                                sapply(LUAD123_testgene_cluster3[[1]], function(x) x[[5]]),
                                sapply(LUAD455_testgene_cluster3[[1]], function(x) x[[5]]))

write.csv(formatC(LUAD_testgene_cluster3, format = "e", digits = 2),file="../../Desktop/ge_c3.csv")
write.csv(formatC(LUAD_testgene_cluster2, format = "e", digits = 2),file="../../Desktop/ge_c2.csv")

#methylation data
cpyfeatures <- LUAD74_methy_cluster3
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
LUAD74_testmethy_cluster2 <- cluster_measure_test(methydata74,2,add_dire74[[2]],0,methysurv)
LUAD74_testmethy_cluster3 <- cluster_measure_test(methydata74,3,add_dire74[[2]],0,methysurv)

cpyfeatures <- LUAD123_methy_cluster3
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
LUAD123_testmethy_cluster2 <- cluster_measure_test(methydata123,2,add_dire123[[2]],0,methysurv)
LUAD123_testmethy_cluster3 <- cluster_measure_test(methydata123,3,add_dire123[[2]],0,methysurv)

cpyfeatures <- LUAD455_methy_cluster3
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
LUAD455_testmethy_cluster2 <- cluster_measure_test(methydata455,2,add_dire455[[2]],0,methysurv)
LUAD455_testmethy_cluster3 <- cluster_measure_test(methydata455,3,add_dire455[[2]],0,methysurv)

LUAD_testmethy_cluster2 <- rbind(sapply(LUAD74_testmethy_cluster2[[1]], function(x) x[[4]]),
                                sapply(LUAD123_testmethy_cluster2[[1]], function(x) x[[4]]),
                                sapply(LUAD455_testmethy_cluster2[[1]], function(x) x[[4]]))

LUAD_testmethy_cluster3 <- rbind(sapply(LUAD74_testmethy_cluster3[[1]], function(x) x[[5]]),
                                sapply(LUAD123_testmethy_cluster3[[1]], function(x) x[[5]]),
                                sapply(LUAD455_testmethy_cluster3[[1]], function(x) x[[5]]))
write.csv(formatC(LUAD_testmethy_cluster3, format = "e", digits = 2),file="../../Desktop/dm_c3.csv")
write.csv(formatC(LUAD_testmethy_cluster2, format = "e", digits = 2),file="../../Desktop/dm_c2.csv")

# combined FSF

datalist121nor74 <- list(genedata74,methydata74[match(genedata74$sampleID,methydata74$sampleID),])
datalist121nor123 <- list(genedata123,methydata123[match(genedata123$sampleID,methydata123$sampleID),])
datalist121nor455 <- list(genedata455,methydata455[match(genedata455$sampleID,methydata455$sampleID),])
level1_74_features <- getselectedfeatures(LUAD74_geneexp_cluster3)
level2_74_features <- getselectedfeatures(LUAD74_methy_cluster3)
featurelistnor74 <- list(level1_74_features,level2_74_features)

level1_123_features <- getselectedfeatures(LUAD123_geneexp_cluster3)
level2_123_features <- getselectedfeatures(LUAD123_methy_cluster3)
featurelistnor123 <- list(level1_123_features,level2_123_features)

level1_455_features <- getselectedfeatures(LUAD455_geneexp_cluster3)
level2_455_features <- getselectedfeatures(LUAD455_methy_cluster3)
featurelistnor455 <- list(level1_455_features,level2_455_features)

dmmatchdata74 <- methydata74[match(genedata74$sampleID,methydata74$sampleID),]
dmmatchdata123 <- methydata123[match(genedata123$sampleID,methydata123$sampleID),]
dmmatchdata455 <- methydata455[match(genedata455$sampleID,methydata455$sampleID),]

#hierarchical clustering
combined74nor_3cluster <- clustering_test(datalist121nor74,featurelistnor74,3,10,7,7,5,add_dire74,0)
combined123nor_3cluster <- clustering_test(datalist121nor123,featurelistnor123,3,10,7,7,5,add_dire123,0)
combined455nor_3cluster <- clustering_test(datalist121nor455,featurelistnor455,3,10,7,7,5,add_dire455,0)

combined_3cluster <- rbind(combined74nor_3cluster[[1]],combined123nor_3cluster[[1]],combined455nor_3cluster[[1]])

combined74nor_2cluster <- clustering_test(datalist121nor74,featurelistnor74,2,10,7,7,5,add_dire74,0)
combined123nor_2cluster <- clustering_test(datalist121nor123,featurelistnor123,2,10,7,7,5,add_dire123,0)
combined455nor_2cluster <- clustering_test(datalist121nor455,featurelistnor455,2,10,7,7,5,add_dire455,0)
combined_2cluster <- rbind(combined74nor_2cluster[[1]],combined123nor_2cluster[[1]],combined455nor_2cluster[[1]])
write.csv(formatC(combined_3cluster, format = "e", digits = 2),file="../../Desktop/combine_c3.csv")
write.csv(formatC(combined_2cluster, format = "e", digits = 2),file="../../Desktop/combine_c2.csv")

#kmeans
# clustering_test(datalist121nor74,featurelistnor74,3,10,7,7,5,add_dire74,1)
# clustering_test(datalist121nor123,featurelistnor123,3,10,7,7,5,add_dire123,1)
# clustering_test(datalist121nor455,featurelistnor455,3,10,7,7,5,add_dire455,1)

# FSF from multiplex
nothresdata455 <- cbind(genedata455,methydata455[match(genedata455$sampleID,methydata455$sampleID),-1])
nothresdata123 <- cbind(genedata123,methydata123[match(genedata123$sampleID,methydata123$sampleID),-1])
nothresdata74 <- cbind(genedata74,methydata74[match(genedata74$sampleID,methydata74$sampleID),-1])

  
cpyfeatures <- LUAD455_gedm_cluster2_1
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor455cluster2_1 <- sapply(cluster_measure_test(nothresdata455,2,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor455cluster3_1 <- sapply(cluster_measure_test(nothresdata455,3,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD455_gedm_cluster2_2
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
# this one is used for plotting.
nor455cluster2_2 <- sapply(cluster_measure_test(nothresdata455,2,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor455cluster3_2 <- sapply(cluster_measure_test(nothresdata455,3,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD455_gedm_cluster2_all
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor455cluster2_all <- sapply(cluster_measure_test(nothresdata455,2,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor455cluster3_all <- sapply(cluster_measure_test(nothresdata455,3,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD455ge_cluster2
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor455cluster2_ge <- sapply(cluster_measure_test(genedata455,2,add_dire455[[1]],geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor455cluster3_ge <- sapply(cluster_measure_test(genedata455,3,add_dire455[[1]],geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD455dm_cluster2 
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor455cluster2_dm <- sapply(cluster_measure_test(dmmatchdata455,2,add_dire455[[2]],geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor455cluster3_dm <- sapply(cluster_measure_test(dmmatchdata455,3,add_dire455[[2]],geneindicator=0,genesurv)[[1]],function(x) x[5:8])


cpyfeatures <- LUAD455_gedm_cluster2_scoring1 
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor455cluster2_scoring1 <- sapply(cluster_measure_test(nothresdata455,2,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor455cluster3_scoring1 <- sapply(cluster_measure_test(nothresdata455,3,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])


cpyfeatures <- LUAD455_gedm_cluster2_scoring5
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor455cluster2_scoring5 <- sapply(cluster_measure_test(nothresdata455,2,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor455cluster3_scoring5 <- sapply(cluster_measure_test(nothresdata455,3,c(add_dire455[[1]],add_dire455[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])


cpyfeatures <- LUAD74_gedm_cluster2_1
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor74cluster2_1 <- sapply(cluster_measure_test(nothresdata74,2,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor74cluster3_1 <- sapply(cluster_measure_test(nothresdata74,3,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD74_gedm_cluster2_2
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor74cluster2_2 <- sapply(cluster_measure_test(nothresdata74,2,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor74cluster3_2 <- sapply(cluster_measure_test(nothresdata74,3,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD74_gedm_cluster2_all
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor74cluster2_all <- sapply(cluster_measure_test(nothresdata74,2,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor74cluster3_all <- sapply(cluster_measure_test(nothresdata74,3,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD74ge_cluster2
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor74cluster2_ge <- sapply(cluster_measure_test(genedata74,2,add_dire74[[1]],geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor74cluster3_ge <- sapply(cluster_measure_test(genedata74,3,add_dire74[[1]],geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD74dm_cluster2 
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor74cluster2_dm <- sapply(cluster_measure_test(dmmatchdata74,2,add_dire74[[2]],geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor74cluster3_dm <- sapply(cluster_measure_test(dmmatchdata74,3,add_dire74[[2]],geneindicator=0,genesurv)[[1]],function(x) x[5:8])


cpyfeatures <- LUAD74_gedm_cluster2_scoring1 
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor74cluster2_scoring1 <- sapply(cluster_measure_test(nothresdata74,2,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor74cluster3_scoring1 <- sapply(cluster_measure_test(nothresdata74,3,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])


cpyfeatures <- LUAD74_gedm_cluster2_scoring5
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor74cluster2_scoring5 <- sapply(cluster_measure_test(nothresdata74,2,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor74cluster3_scoring5 <- sapply(cluster_measure_test(nothresdata74,3,c(add_dire74[[1]],add_dire74[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])


cpyfeatures <- LUAD123_gedm_cluster2_1
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor123cluster2_1 <- sapply(cluster_measure_test(nothresdata123,2,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor123cluster3_1 <- sapply(cluster_measure_test(nothresdata123,3,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD123_gedm_cluster2_2
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor123cluster2_2 <- sapply(cluster_measure_test(nothresdata123,2,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor123cluster3_2 <- sapply(cluster_measure_test(nothresdata123,3,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD123_gedm_cluster2_all
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor123cluster2_all <- sapply(cluster_measure_test(nothresdata123,2,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor123cluster3_all <- sapply(cluster_measure_test(nothresdata123,3,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD123ge_cluster2
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor123cluster2_ge <- sapply(cluster_measure_test(genedata123,2,add_dire123[[1]],geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor123cluster3_ge <- sapply(cluster_measure_test(genedata123,3,add_dire123[[1]],geneindicator=0,genesurv)[[1]],function(x) x[5:8])

cpyfeatures <- LUAD123dm_cluster2 
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor123cluster2_dm <- sapply(cluster_measure_test(dmmatchdata123,2,add_dire123[[2]],geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor123cluster3_dm <- sapply(cluster_measure_test(dmmatchdata123,3,add_dire123[[2]],geneindicator=0,genesurv)[[1]],function(x) x[5:8])


cpyfeatures <- LUAD123_gedm_cluster2_scoring1 
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor123cluster2_scoring1 <- sapply(cluster_measure_test(nothresdata123,2,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor123cluster3_scoring1 <- sapply(cluster_measure_test(nothresdata123,3,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])


cpyfeatures <- LUAD123_gedm_cluster2_scoring5
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
source('../../Desktop/Experiments_largeEMT_correct/retrievefeatures.R')
nor123cluster2_scoring5 <- sapply(cluster_measure_test(nothresdata123,2,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[4:6])
nor123cluster3_scoring5 <- sapply(cluster_measure_test(nothresdata123,3,c(add_dire123[[1]],add_dire123[[2]]),geneindicator=0,genesurv)[[1]],function(x) x[5:8])

multiplex_nor74_2cluster <- list(nor74cluster2_1,nor74cluster2_2,nor74cluster2_all,nor74cluster2_dm,nor74cluster2_ge,nor74cluster2_scoring1,nor74cluster2_scoring5)
multiplex_nor74_3cluster <- list(nor74cluster3_1,nor74cluster3_2,nor74cluster3_all,nor74cluster3_dm,nor74cluster3_ge,nor74cluster3_scoring1,nor74cluster3_scoring5)
multiplex_nor123_2cluster <- list(nor123cluster2_1,nor123cluster2_2,nor123cluster2_all,nor123cluster2_dm,nor123cluster2_ge,nor123cluster2_scoring1,nor123cluster2_scoring5)
multiplex_nor123_3cluster <- list(nor123cluster3_1,nor123cluster3_2,nor123cluster3_all,nor123cluster3_dm,nor123cluster3_ge,nor123cluster3_scoring1,nor123cluster3_scoring5)
multiplex_nor455_2cluster <- list(nor455cluster2_1,nor455cluster2_2,nor455cluster2_all,nor455cluster2_dm,nor455cluster2_ge,nor455cluster2_scoring1,nor455cluster2_scoring5)
multiplex_nor455_3cluster <- list(nor455cluster3_1,nor455cluster3_2,nor455cluster3_all,nor455cluster3_dm,nor455cluster3_ge,nor455cluster3_scoring1,nor455cluster3_scoring5)

temp <- sapply(multiplex_nor74_3cluster, function(x) sapply(1:12, function(y) x[1,y]))[,c(1,2,3,6,7)]
nor74multiplex3c <- apply(temp,1,min)
temp <- sapply(multiplex_nor123_3cluster, function(x) sapply(1:12, function(y) x[1,y]))[,c(1,2,3,6,7)]
nor123multiplex3c <- apply(temp,1,min)
temp <- sapply(multiplex_nor455_3cluster, function(x) sapply(1:12, function(y) x[1,y]))[,c(1,2,3,6,7)]
nor455multiplex3c <- apply(temp,1,min)

temp <- sapply(multiplex_nor74_2cluster, function(x) sapply(1:12, function(y) x[1,y]))[,c(1,2,3,6,7)]
nor74multiplex2c <- apply(temp,1,min)
temp <- sapply(multiplex_nor123_2cluster, function(x) sapply(1:12, function(y) x[1,y]))[,c(1,2,3,6,7)]
nor123multiplex2c <- apply(temp,1,min)
temp <- sapply(multiplex_nor455_2cluster, function(x) sapply(1:12, function(y) x[1,y]))[,c(1,2,3,6,7)]
nor455multiplex2c <- apply(temp,1,min)

normultiplex2c <- rbind(nor74multiplex2c,nor123multiplex2c,nor455multiplex2c)
row.names(normultiplex2c) <- c("74 nodes","123 nodes", "455 nodes")
colnames(normultiplex2c) <- c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble","all")
write.csv(formatC(normultiplex2c, format = "e", digits = 2),file="../../Desktop/multiplex_c2.csv")

normultiplex3c <- rbind(nor74multiplex3c,nor123multiplex3c,nor455multiplex3c)
row.names(normultiplex3c) <- c("74 nodes","123 nodes", "455 nodes")
colnames(normultiplex3c) <- c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble","all")
write.csv(formatC(normultiplex3c, format = "e", digits = 2),file="../../Desktop/multiplex_c3.csv")

# pvalues for FSFs, for network features, alltogether and indicidual subnets. 
cluster_measure_test <- function(nothresdata,numcluster,directionadd,geneindicator=0,survobj){
  if(geneindicator==0){
    originalclr <- tsne_km_cl_test(nothresdata,1:(ncol(nothresdata)-2),0,algorithmstrings[i],numcluster,survobj)
  }else{
    temp <- 1:(ncol(nothresdata)-2)
    originalclr <- tsne_km_cl_test(nothresdata,temp[-mirrange],0,algorithmstrings[i],numcluster,survobj)
  }
    
  i=1
  ttestclr <- tsne_km_cl_test(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=2
  lassoclr <- tsne_km_cl_test(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=3
  netlassoclr <- tsne_km_cl_test(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=4
  addgfs <- addg_sf[1:numtopnets]
  # addgfs <- addg_sf[1:20]
  moldata <- nothresdata[,-1]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[x])
  newX2 <- sapply(addgfs, function(x) {if(length(x) > 1) {rowSums(newX[,x])/length(x)} else newX[,x]})
  nothresdatanet <- data.frame(nothresdata[,1],newX2)
  addgclr <- tsne_km_cl_test(nothresdatanet,1:ncol(newX2),0,algorithmstrings[i],numcluster,survobj)
  i=5
  netrankclr <- tsne_km_cl_test(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=6
  stsvmclr <- tsne_km_cl_test(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=7
  coxclr <- tsne_km_cl_test(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=8
  coxregclr <- tsne_km_cl_test(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=9
  rdsclr <- tsne_km_cl_test(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=10
  survnetfs <- survnet_sf[1:numtopnets]
  # survnetfs <- survnet_sf[1:20]
  newX2 <- sapply(survnetfs,function(x) {if(length(x) > 1) {rowSums(moldata[,x])/length(x)} else newX[,x]})
  nothresdatanet <- data.frame(nothresdata[,1],newX2)
  survnetclr <- tsne_km_cl_test(nothresdatanet,1:ncol(newX2),0,algorithmstrings[i],numcluster,survobj)
  i=11
  allrankclr <- tsne_km_cl_test(nothresdata,alg10sfs,0,"all",numcluster,survobj)
  
  numbers <- list(ttestclr,lassoclr,netlassoclr,addgclr,netrankclr,stsvmclr,coxclr,coxregclr,rdsclr,survnetclr,allrankclr,originalclr)
  features <- list(ttestfs,lassofs,netlassofs, addgfs, netrankfs, stsvmfs, coxfs, coxregfs, rdsfs, survnetfs,alg10sfs)
  list(numbers,features)
}

tsne_km_cl_test <- function(nothresdata, features, usetsne=1, algorithm,numcluster,survobj){

  data <- nothresdata[,features+1]
  
  # tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
  # if(usetsne==1)
  #   d_tsne_1 = as.data.frame(tsne_model_1$Y) else
  #     d_tsne_1 = data
  
  d_tsne_1_original=NULL
  fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
  d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit_cluster_hierarchical=hclust(dist(scale(data)),method="average")
  d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=numcluster))  
  
  #for saving plot with correct names. 
  
  # fit <- survfit(survobj ~ cl_kmeans,data=d_tsne_1_original)
  # # filestring <- paste('../Thesis_Results/alllevels/clustering/', ncol(nothresdata)-2, algorithm, '_kmeans.pdf', sep='')
  # # pdf(filestring, width = 6.5,height = 5)
  # ggsurvplot(
  #   fit,                     # survfit object with calculated statistics.
  #   #data = d_tsne_1_original,  # data used to fit survival curves.
  #   risk.table = TRUE,       # show risk table.
  #   pval = paste("p=",format(pval_kmeans,digits=3),sep=""),             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for
  #   # point estimaes of survival curves.
  #   xlim = c(0,2000),        # present narrower X axis, but not affect
  #   # survival estimates.
  #   break.time.by = 500,     # break X axis in time intervals by 500.
  #   ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  #   risk.table.y.text.col = T, # colour risk table text annotations.
  #   risk.table.y.text = FALSE # show bars instead of names in text annotations
  #   # in legend of risk table
  #   # title = "Lasso"
  # )
  # dev.off()
  # 
  # fit <- survfit(survobj ~ cl_hierarchical,data=d_tsne_1_original)
  # filestring <- paste('../Thesis_Results/alllevels/clustering/', ncol(nothresdata)-2, algorithm, '_hierarchical.pdf', sep='')
  # pdf(filestring, width = 6.5,height = 5)
  # ggsurvplot(
  #   fit,                     # survfit object with calculated statistics.
  #   data = d_tsne_1_original,  # data used to fit survival curves.
  #   risk.table = TRUE,       # show risk table.
  #   pval = TRUE,             # show p-value of log-rank test.
  #   # pval = TRUE,             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for
  #   # point estimaes of survival curves.
  #   xlim = c(0,90),        # present narrower X axis, but not affect
  #   # survival estimates.
  #   break.time.by = 20,     # break X axis in time intervals by 500.
  #   ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  #   risk.table.y.text.col = T, # colour risk table text annotations.
  #   risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  # )
  # dev.off()
  
  # log-rank test
  sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_hierarchical)$chisq
  pval_hierarchical <- 1- pchisq(sdiff,numcluster-1) 
  c(pval_kmeans,as.numeric(table(d_tsne_1_original$cl_kmeans)),pval_hierarchical,as.numeric(table(d_tsne_1_original$cl_hierarchical)))
}

cluster_measure_test_withcli <- function(nothresdata,numcluster,directionadd,geneindicator=0,survobj){
  if(geneindicator==0){
    originalclr <- tsne_km_cl_test_withcli(nothresdata,1:(ncol(nothresdata)-2),0,algorithmstrings[i],numcluster,survobj)
  }else{
    temp <- 1:(ncol(nothresdata)-2)
    originalclr <- tsne_km_cl_test_withcli(nothresdata,temp[-mirrange],0,algorithmstrings[i],numcluster,survobj)
  }
  
  i=1
  ttestclr <- tsne_km_cl_test_withcli(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=2
  lassoclr <- tsne_km_cl_test_withcli(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=3
  netlassoclr <- tsne_km_cl_test_withcli(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=4
  addgfs <- addg_sf[1:numtopnets]
  # addgfs <- addg_sf[1:20]
  moldata <- nothresdata[,-1]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[x])
  newX2 <- sapply(addgfs, function(x) {if(length(x) > 1) {rowSums(newX[,x])/length(x)} else newX[,x]})
  nothresdatanet <- data.frame(nothresdata[,1],newX2)
  addgclr <- tsne_km_cl_test_withcli(nothresdatanet,1:ncol(newX2),0,algorithmstrings[i],numcluster,survobj)
  i=5
  netrankclr <- tsne_km_cl_test_withcli(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=6
  stsvmclr <- tsne_km_cl_test_withcli(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=7
  coxclr <- tsne_km_cl_test_withcli(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=8
  coxregclr <- tsne_km_cl_test_withcli(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=9
  rdsclr <- tsne_km_cl_test_withcli(nothresdata,get(algorithmstrings[i]),0,algorithmstrings[i],numcluster,survobj)
  i=10
  survnetfs <- survnet_sf[1:numtopnets]
  # survnetfs <- survnet_sf[1:20]
  newX2 <- sapply(survnetfs,function(x) {if(length(x) > 1) {rowSums(moldata[,x])/length(x)} else newX[,x]})
  nothresdatanet <- data.frame(nothresdata[,1],newX2)
  survnetclr <- tsne_km_cl_test_withcli(nothresdatanet,1:ncol(newX2),0,algorithmstrings[i],numcluster,survobj)
  i=11
  allrankclr <- tsne_km_cl_test_withcli(nothresdata,alg10sfs,0,"all",numcluster,survobj)
  
  numbers <- list(ttestclr,lassoclr,netlassoclr,addgclr,netrankclr,stsvmclr,coxclr,coxregclr,rdsclr,survnetclr,allrankclr,originalclr)
  features <- list(ttestfs,lassofs,netlassofs, addgfs, netrankfs, stsvmfs, coxfs, coxregfs, rdsfs, survnetfs,alg10sfs)
  list(numbers,features)
}

tsne_km_cl_test_withcli <- function(nothresdata, features, usetsne=1, algorithm,numcluster,survobj){
  
  data1 <- nothresdata[,features+1]
  data2 <- norclinical[match(nothresdata[,1],norclinical$GEO_ID),-c(1,7,13,14)]
  data <- data.frame(data1,data2)
  # tsne_model_1 = Rtsne(as.matrix(data), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
  # if(usetsne==1)
  #   d_tsne_1 = as.data.frame(tsne_model_1$Y) else
  #     d_tsne_1 = data
  
  # d_tsne_1_original=NULL
  # fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
  # d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)
  fit_cluster_hierarchical=hclust(daisy((data)))
  d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=numcluster))  
  
  #for saving plot with correct names. 
  
  # fit <- survfit(survobj ~ cl_kmeans,data=d_tsne_1_original)
  # # filestring <- paste('../Thesis_Results/alllevels/clustering/', ncol(nothresdata)-2, algorithm, '_kmeans.pdf', sep='')
  # # pdf(filestring, width = 6.5,height = 5)
  # ggsurvplot(
  #   fit,                     # survfit object with calculated statistics.
  #   #data = d_tsne_1_original,  # data used to fit survival curves.
  #   risk.table = TRUE,       # show risk table.
  #   pval = paste("p=",format(pval_kmeans,digits=3),sep=""),             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for
  #   # point estimaes of survival curves.
  #   xlim = c(0,2000),        # present narrower X axis, but not affect
  #   # survival estimates.
  #   break.time.by = 500,     # break X axis in time intervals by 500.
  #   ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  #   risk.table.y.text.col = T, # colour risk table text annotations.
  #   risk.table.y.text = FALSE # show bars instead of names in text annotations
  #   # in legend of risk table
  #   # title = "Lasso"
  # )
  # dev.off()
  # 
  # fit <- survfit(survobj ~ cl_hierarchical,data=d_tsne_1_original)
  # filestring <- paste('../Thesis_Results/alllevels/clustering/', ncol(nothresdata)-2, algorithm, '_hierarchical.pdf', sep='')
  # pdf(filestring, width = 6.5,height = 5)
  # ggsurvplot(
  #   fit,                     # survfit object with calculated statistics.
  #   data = d_tsne_1_original,  # data used to fit survival curves.
  #   risk.table = TRUE,       # show risk table.
  #   pval = TRUE,             # show p-value of log-rank test.
  #   # pval = TRUE,             # show p-value of log-rank test.
  #   conf.int = TRUE,         # show confidence intervals for
  #   # point estimaes of survival curves.
  #   xlim = c(0,90),        # present narrower X axis, but not affect
  #   # survival estimates.
  #   break.time.by = 20,     # break X axis in time intervals by 500.
  #   ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  #   risk.table.y.text.col = T, # colour risk table text annotations.
  #   risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  # )
  # dev.off()
  
  # log-rank test
  # sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_kmeans)$chisq
  # pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  sdiff <- survdiff(survobj ~ d_tsne_1_original$cl_hierarchical)$chisq
  pval_hierarchical <- 1- pchisq(sdiff,numcluster-1) 
  c(pval_hierarchical,as.numeric(table(d_tsne_1_original$cl_hierarchical)))
}

#clustering_test(datalist121nor455,featurelistnor455,2,10,7,7,5,add_dire455,0)
clustering_test <- function(datalist, featurelist, numcluster, numf2dt,numf3dt,numf2dtnet, numf3dtnet,directionadd,kmeansindicator=0){

  survobj <- genesurv
  pvaluesmatrix <- matrix(0,3,12)
  sizematrix1 <- matrix(0,numcluster,12)
  sizematrix2 <- matrix(0,numcluster,12)
  sizematrix3 <- matrix(0,numcluster,12)
  sizelist <- list(sizematrix1,sizematrix2,sizematrix3)
  for(j in c(1:3,5:9,11)){
    for(i in 1:2){
      data <- datalist[[i]][,featurelist[[i]][[j]]+1]
      if(kmeansindicator ==0){
        fit_cluster_kmeans=hclust(dist(scale(data)))
        cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
      }else{
        fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
        cl_kmeans = factor(fit_cluster_kmeans$cluster)
      }
      fit <- survfit(survobj ~ cl_kmeans)
      sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
      pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
      # print(pval_kmeans)
      pvaluesmatrix[i,j]  <- pval_kmeans
      sizelist[[i]][,j] <- as.numeric(table(cl_kmeans))
    }
    i=3
    # combine GE and DM
    Data1 <- datalist[[1]][,featurelist[[1]][[j]][1:numf2dt]+1]
    Data2 <- datalist[[2]][,featurelist[[2]][[j]][1:numf2dt]+1]
    data <- data.frame(Data1,Data2)
    if(kmeansindicator ==0){
      fit_cluster_kmeans=hclust(dist(scale(data)))
      cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    }else{
      fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
      cl_kmeans = factor(fit_cluster_kmeans$cluster)
    }
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    sizelist[[i]][,j] <- as.numeric(table(cl_kmeans))
  }
  
  j=4
  moldata <- datalist[[1]][,-1]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[1]][x])
  dataadd1 <- sapply(featurelist[[1]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  moldata <- datalist[[2]][,-1]
  newX <- sapply(1:ncol(moldata), function(x) moldata[,x]*directionadd[[2]][x])
  dataadd2 <- sapply(featurelist[[2]][[4]],function(x) rowSums(newX[,x])/length(x))
  
  dataaddlist <- list(dataadd1,dataadd2)
  
  for(i in 1:2){
    data <- dataaddlist[[i]]
    if(kmeansindicator ==0){
      fit_cluster_kmeans=hclust(dist(scale(data)))
      cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    }else{
      fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
      cl_kmeans = factor(fit_cluster_kmeans$cluster)
    }
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    sizelist[[i]][,j] <- as.numeric(table(cl_kmeans))
  }
  i=3
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  if(kmeansindicator ==0){
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  }else{
    fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
  }
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  sizelist[[i]][,j] <- as.numeric(table(cl_kmeans))
  
  
  
  j=10
  
  moldata <- datalist[[1]][,-1]
  dataadd1 <- sapply(featurelist[[1]][[10]],function(x) rowSums(moldata[,x])/length(x))
  
  moldata <- datalist[[2]][,-1]
  dataadd2 <- sapply(featurelist[[2]][[10]],function(x) rowSums(moldata[,x])/length(x))

  
  dataaddlist <- list(dataadd1,dataadd2)
  
  for(i in 1:2){
    data <- dataaddlist[[i]]
    if(kmeansindicator ==0){
      fit_cluster_kmeans=hclust(dist(scale(data)))
      cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    }else{
      fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
      cl_kmeans = factor(fit_cluster_kmeans$cluster)
    }
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    sizelist[[i]][,j] <- as.numeric(table(cl_kmeans))
  }
  i=3
  # combine GE and DM
  Data1 <- dataaddlist[[1]][,1:numf2dtnet]
  Data2 <- dataaddlist[[2]][,1:numf2dtnet]
  data <- data.frame(Data1,Data2)
  if(kmeansindicator ==0){
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  }else{
    fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
  }
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  sizelist[[i]][,j] <- as.numeric(table(cl_kmeans))
 
  
  j=12
  featurelist <- list(c(1:(ncol(datalist[[1]])-1)),c(1:(ncol(datalist[[2]])-1)))
  for(i in 1:2){
    data <- datalist[[i]][,featurelist[[i]]+1]
    if(kmeansindicator ==0){
      fit_cluster_kmeans=hclust(dist(scale(data)))
      cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
    }else{
      fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
      cl_kmeans = factor(fit_cluster_kmeans$cluster)
    }
    fit <- survfit(survobj ~ cl_kmeans)
    sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
    pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
    # print(pval_kmeans)
    pvaluesmatrix[i,j]  <- pval_kmeans
    sizelist[[i]][,j] <- as.numeric(table(cl_kmeans))
  }
  i=3
  # combine GE and DM
  Data1 <- datalist[[1]][,featurelist[[1]]+1]
  Data2 <- datalist[[2]][,featurelist[[2]]+1]
  data <- data.frame(Data1,Data2)
  if(kmeansindicator ==0){
    fit_cluster_kmeans=hclust(dist(scale(data)))
    cl_kmeans =  factor(cutree(fit_cluster_kmeans, k=numcluster))
  }else{
    fit_cluster_kmeans=kmeans(data, numcluster,iter.max = 100)  
    cl_kmeans = factor(fit_cluster_kmeans$cluster)
  }
  fit <- survfit(survobj ~ cl_kmeans)
  sdiff <- survdiff(survobj ~ cl_kmeans)$chisq
  pval_kmeans <- 1- pchisq(sdiff,numcluster-1) 
  # print(pval_kmeans)
  pvaluesmatrix[i,j]  <- pval_kmeans
  sizelist[[i]][,j] <- as.numeric(table(cl_kmeans))

  list(pvaluesmatrix,sizelist)
}

save.image('test/testdata.RData')