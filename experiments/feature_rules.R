#level1: gene expression, level2: DNA methylation, level3: CNA
#level4: protein level5: clinical 
# commonsamples <- intersect(intersect(datatypes74[[1]][,1],datatypes74[[2]][,1]),datatypes74[[3]][,1])
# commonsamples_nothres <- intersect(intersect(level1_74_nothres[,1],level2_74_nothres[,1]),level3_74_nothres[,1])
# commonsamples_nothres_early <- intersect(intersect(level1_earlystage[,1],level2_earlystage[,1]), level3_earlystage[,1])
# ----------------------------prepare features---------------------------------------------#
getselectedfeatures <- function(cpyfeatures) {
  topfeatures <- 20
  numtopnets <- 10
  ttestfs <- cpyfeatures[[2]][[1]][1:topfeatures]
  #lasso
  lassofs <- cpyfeatures[[2]][[2]][1:topfeatures]
  #netowrk lasso
  netlassofs <- cpyfeatures[[2]][[3]][1:topfeatures]
  #additive subnetworks
  addgfs <- cpyfeatures[[2]][[4]][1:numtopnets]
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
  # coxlassotemp <- length(unique(unlist(cpyfeatures[[2]][[10]][1:numtopnets])))
  # numtofs <- if(topfeatures > coxlassotemp) coxlassotemp else topfeatures
  survnetfs <-cpyfeatures[[2]][[10]][1:numtopnets]
  
  temp <- sort(table(c(ttestfs,lassofs,netlassofs, unique(unlist(addgfs)), netrankfs, stsvmfs, coxfs, coxregfs, rdsfs, 
                       unique(unlist(survnetfs)))),decreasing = TRUE)[1:topfeatures]
  alg10sfs <- as.numeric(names(temp))
  list(ttestfs,lassofs,netlassofs,addgfs,netrankfs,stsvmfs,coxfs,coxregfs,rdsfs,survnetfs,alg10sfs)
}


#----------------------------prepare data-------------------------------------------------#
level1_74 <- read.csv('../Thesis_Results/74network/emt74_geneexp_700_1400.csv')
level2_74 <-  read.csv('../Thesis_Results/74network/emt74_methy_700_1400.csv')
level3_74 <- read.csv('../Thesis_Results/74network/emt70_cna_700_1400.csv')
level1_123 <- read.csv('../Thesis_Results/123network/emt123_geneexp_700_1400.csv')
level2_123 <- read.csv('../Thesis_Results/123network/emt123_methy_700_1400.csv')
level3_123 <- read.csv('../Thesis_Results/123network/emt117_cna_700_1400.csv')
level1_455 <-  read.csv('../Thesis_Results/455network/emt455_geneexp_700_1400.csv')
level2_455 <- read.csv('../Thesis_Results/455network/emt455_methy_700_1400.csv')
level3_455 <-  read.csv('../Thesis_Results/455network/emt445_cna_700_1400.csv')

level1_74_nothres <- read.csv('../Thesis_Results/74network/emt74_geneexp_700_1400_nothres.csv')
level2_74_nothres <-  read.csv('../Thesis_Results/74network/emt74_methy_700_1400_nothres.csv')
level3_74_nothres <- read.csv('../Thesis_Results/74network/emt70_cna_700_1400_nothres.csv')
level1_123_nothres <- read.csv('../Thesis_Results/123network/emt123_geneexp_700_1400_nothres.csv')
level2_123_nothres <- read.csv('../Thesis_Results/123network/emt123_methy_700_1400_nothres.csv')
level3_123_nothres <- read.csv('../Thesis_Results/123network/emt117_cna_700_1400_nothres.csv')
level1_455_nothres <-  read.csv('../Thesis_Results/455network/emt455_geneexp_700_1400_nothres.csv')
level2_455_nothres <- read.csv('../Thesis_Results/455network/emt455_methy_700_1400_nothres.csv')
level3_455_nothres <-  read.csv('../Thesis_Results/455network/emt445_cna_700_1400_nothres.csv')

dataproteins <- read.csv('../Experiments/LUAD/Data/protein/all_case_protein_365.txt', sep='\t')
rowmatch <- which(!is.na(match(samplesthres[,1],dataproteins[,1])))
colmatch <-match(V(LEnet74)$name,names(dataproteins))
colmatch <- colmatch[!is.na(colmatch)]
level4_74 <- dataproteins[rowmatch,c(1,colmatch)]
colmatch <- match(V(LEnet123)$name,names(dataproteins))
colmatch <- colmatch[!is.na(colmatch)]
level4_123 <- dataproteins[rowmatch,c(1,colmatch)]
colmatch <- match(V(LEnet455)$name,names(dataproteins))
colmatch <- colmatch[!is.na(colmatch)]
level4_455 <- dataproteins[rowmatch,c(1,colmatch)]
level5 <- clinicalnona

datatypes74 <- list(level1_74,level2_74,level3_74,level4_74,level5)
datatypes123 <- list(level1_123,level2_123,level3_123,level4_123,level5)
datatypes455 <- list(level1_455,level2_455,level3_455,level4_455,level5)
# datatypes74_common3 <- list(level1_74[match(commonsamples,level1_74[,1]),],level2_74[match(commonsamples,level2_74[,1]),],
#                             level3_74[match(commonsamples,level3_74[,1]),])
# datatypes123_common3 <- list(level1_123[match(commonsamples,level1_123[,1]),],level2_123[match(commonsamples,level2_123[,1]),],
#                             level3_123[match(commonsamples,level3_123[,1]),])
# datatypes455_common3 <- list(level1_455[match(commonsamples,level1_455[,1]),],level2_455[match(commonsamples,level2_455[,1]),],
#                             level3_455[match(commonsamples,level3_455[,1]),])

datatypes74_common3 <- list(level1_74_nothres[match(commonsamples_nothres,level1_74_nothres[,1]),],level2_74_nothres[match(commonsamples_nothres,level2_74_nothres[,1]),],
                            level3_74_nothres[match(commonsamples_nothres,level3_74_nothres[,1]),])
datatypes123_common3 <- list(level1_123_nothres[match(commonsamples_nothres,level1_123_nothres[,1]),],level2_123_nothres[match(commonsamples_nothres,level2_123_nothres[,1]),],
                             level3_123_nothres[match(commonsamples_nothres,level3_123_nothres[,1]),])
datatypes455_common3 <- list(level1_455_nothres[match(commonsamples_nothres,level1_455_nothres[,1]),],level2_455_nothres[match(commonsamples_nothres,level2_455_nothres[,1]),],
                             level3_455_nothres[match(commonsamples_nothres,level3_455_nothres[,1]),])
earlys_datatypes455_common3 <- list(level1_455_nothres[match(commonsamples_nothres_early,level1_455_nothres[,1]),],level2_455_nothres[match(commonsamples_nothres_early,level2_455_nothres[,1]),],
                             level3_455_nothres[match(commonsamples_nothres_early,level3_455_nothres[,1]),])

level1_74_features <- getselectedfeatures(LUAD74_geneexp_cluster3)
level2_74_features <- getselectedfeatures(LUAD74_methy_cluster3)
level3_74_features <- getselectedfeatures(LUAD74_cna_cluster3)
listfeatures74 <- list(level1_74_features,level2_74_features,level3_74_features)

level1_123_features <- getselectedfeatures(LUAD123_geneexp_cluster3)
level2_123_features <- getselectedfeatures(LUAD123_methy_cluster3)
level3_123_features <- getselectedfeatures(LUAD123_cna_cluster3)
listfeatures123 <- list(level1_123_features,level2_123_features,level3_123_features)

level1_455_features <- getselectedfeatures(LUAD455_geneexp_cluster3)
level2_455_features <- getselectedfeatures(LUAD455_methy_cluster3)
level3_455_features <- getselectedfeatures(LUAD455_cna_cluster3)
listfeatures455 <- list(level1_455_features,level2_455_features,level3_455_features)
# molnames <- c("CDH1","CDH2","SNAI1","SNAI2","TP53","ZEB1","ZEB2")
# mirnanames <- c("miR-141","miR-200a","miR-200b","miR-200c","miR-34a","miR-34b","miR-34c")
# ruleidx1 <- match(molnames,names(datatypes[[1]]))
# ruleidx2 <- match(mirnanames,names(datatypes[[2]]))
# ruleidx3 <- c(match(c(molnames,mirnanames),names(datatypes[[4]])))

#random variables
commmolecules <- intersect(intersect(names(RNA_Seq_case_all),names(cna_case_all)),names(methy_case_all))
set.seed(1)
random84 <- sample(commmolecules, 84)
random68 <- sample(random84, 68)
randmir16 <- sample(names(miRNA_case_all),16)
commproteins <- intersect(commmolecules, names(protein_case_all))
protein_randidx <- match(sample(commproteins,4),names(protein_case_all))

RNA_Seq_case_se_rand <- RNA_Seq_case_all[allrnaseq[,1] %in% samplesthres$samplecode, 
                                         c(1,match(random68,names(RNA_Seq_case_all)))]
methy_case_se_rand <- methy_case_all[methy_case_all$samplecode %in% samplesthres$samplecode, 
                                     c(1,match(random84,names(methy_case_all)))]
# cna_case_se_rand <- cna_case_all[cna_case_all$samplecode %in% samplesthres$samplecode, 
# c(1,match(random68,names(cna_case_all)))]
cna_case_se_rand <- cna_case_all[match(cna_case_se$samplecode,cna_case_all$samplecode),c(1,match(random68,names(cna_case_all)))]
protein_case_se_rand <- protein_case_all[protein_case_all$samplecode %in% samplesthres$samplecode, 
                                         c(1,protein_randidx)]

emtdatalist <- function(datatypes, featureslist) {
  datatypes[[1]][,-1] <- scale(datatypes[[1]][,-1])
  datatypes[[2]][,-1] <- scale(datatypes[[2]][,-1])
  datatypes[[3]][,-1] <- scale(datatypes[[3]][,-1])
  datatypes[[4]][,-1] <- scale(datatypes[[4]][,-1])
  # datatypes[[5]][,-1] <- scale(datatypes[[5]][,-1])
  for(i in 2:ncol(datatypes[[5]]))
    datatypes[[5]][,i] <- as.factor(datatypes[[5]][,i])
  datalist <- list()
  #add feature index below
  datalist[[1]] <- data.frame(datatypes[[1]][,1],(datatypes[[1]][,featureslist[[1]]+1]>0))
  datalist[[2]] <- data.frame(datatypes[[2]][,1],(datatypes[[2]][,featureslist[[2]]+1]>0))
  datalist[[3]] <- data.frame(datatypes[[3]][,1],(datatypes[[3]][,featureslist[[3]]+1]>0))
  molnames <- unique(c(names(datatypes[[1]][,-1])[featureslist[[1]]],
                names(datatypes[[2]][,-1])[featureslist[[2]]],
                names(datatypes[[3]][,-1])[featureslist[[3]]]))
  ptidx <- match(molnames,names(datatypes[[4]]))
  ptidx <- ptidx[!is.na(ptidx)]
  datalist[[4]] <- data.frame(datatypes[[4]][,1],(datatypes[[4]][,ptidx]>0))
  datalist[[5]] <- datatypes[[5]]
  
  datalist[[1]] <- datalist[[1]][datalist[[1]][,1] %in% samplesthres$samplecode,]
  datalist[[2]] <- datalist[[2]][datalist[[2]][,1] %in% samplesthres$samplecode,]
  datalist[[3]] <- datalist[[3]][datalist[[3]][,1] %in% samplesthres$samplecode,]
  datalist[[4]] <- datalist[[4]][datalist[[4]][,1] %in% samplesthres$samplecode,]
  datalist[[5]] <- datalist[[5]][datalist[[5]][,1] %in% samplesthres$samplecode,]
  datalist
}

randomdatalist <- function() {
  datatypes <- list(RNA_Seq_case_se_rand,miRNA_case_se_rand,cna_case_se_rand,methy_case_se_rand,protein_case_se_rand,clinicalnona)
  datatypes[[1]][,-1] <- scale(datatypes[[1]][,-1])
  datatypes[[2]][,-1] <- scale(datatypes[[2]][,-1])
  datatypes[[3]][,-1] <- scale(datatypes[[3]][,-1])
  datatypes[[4]][,-1] <- scale(datatypes[[4]][,-1])
  datatypes[[5]][,-1] <- scale(datatypes[[5]][,-1])
  for(i in 2:ncol(datatypes[[6]]))
    datatypes[[6]][,i] <- as.factor(datatypes[[6]][,i])
  datalist <- list()
  datalist[[1]] <- data.frame(datatypes[[1]][,1],(datatypes[[1]][,ruleidx1]>0))
  datalist[[2]] <- data.frame(datatypes[[2]][,1],(datatypes[[2]][,ruleidx2]>0))
  datalist[[3]] <- data.frame(datatypes[[3]][,1],(datatypes[[3]][,ruleidx1]>0))
  datalist[[4]] <- data.frame(datatypes[[4]][,1],(datatypes[[4]][,ruleidx3]>0))
  datalist[[5]] <- data.frame(datatypes[[5]][,1],(datatypes[[5]][,-1]>0))
  datalist[[6]] <- datatypes[[6]]
  
  datalist[[1]] <- datalist[[1]][datalist[[1]][,1] %in% samplesthres$samplecode,]
  datalist[[2]] <- datalist[[2]][datalist[[2]][,1] %in% samplesthres$samplecode,]
  datalist[[3]] <- datalist[[3]][datalist[[3]][,1] %in% samplesthres$samplecode,]
  datalist[[4]] <- datalist[[4]][datalist[[4]][,1] %in% samplesthres$samplecode,]
  datalist[[5]] <- datalist[[5]][datalist[[5]][,1] %in% samplesthres$samplecode,]
  datalist[[6]] <- datalist[[6]][datalist[[6]][,1] %in% samplesthres$samplecode,]
  datalist
  datalist
}

simumatrix <- function(data) {
  row <- nrow(data)
  col <- ncol(data)
  matrix <- matrix(0,row,col)
  for(i in 1:col){
    idx <- sample(row,floor(row/2))
    matrix[,i][idx] <- 1 
  }
  matrix>0
}

simudatalist <- function() {
  datatypes <- list(RNA_Seq_case_se[,c(1,rnaseqindex)],miRNA_case_se[,c(1,mirnaindex)],cna_case_se[,c(1,cnaindex)],methy_case_se[,c(1,methyindex)],protein_case_se[,c(1,4,5,19,21)],clinicalnona)
  datalist <- list()
  datalist[[1]] <- data.frame(datatypes[[1]][,1],simumatrix(datatypes[[1]][,ruleidx1]))
  datalist[[2]] <- data.frame(datatypes[[2]][,1],simumatrix(datatypes[[2]][,ruleidx2]))
  datalist[[3]] <- data.frame(datatypes[[3]][,1],simumatrix(datatypes[[3]][,ruleidx1]))
  datalist[[4]] <- data.frame(datatypes[[4]][,1],simumatrix(datatypes[[4]][,ruleidx3]))
  datalist[[5]] <- data.frame(datatypes[[5]][,1],simumatrix(datatypes[[5]][,-1]))
  datalist[[6]] <- datatypes[[6]]
  datalist
}

#type 1 for emt features, type 2 for random features, type 3 for simulated features
datacombiner <- function(vec,type,datatypes,featureslist) {
  if(type == 1) 
    datalist <- emtdatalist(datatypes,featureslist)
  if(type == 2)
    datalist <- randomdatalist()
  if(type == 3)
    datalist <- simudatalist()
  data <- NULL
  days <- NULL
  #return data and days
  # replace samplesthres$samplecode, with common samples code
  if(length(vec)==1){
    samples <- datalist[[vec]][,1]
    matchday <- match(samples,samplesthres$samplecode)
    days <- samplesthres$days[matchday]
    data <- datalist[[vec]]
  }else if (length(vec)==2) {
    samples <- intersect(datalist[[vec[1]]][,1],datalist[[vec[2]]][,1])
    matchday <- match(samples,samplesthres$samplecode)
    days <- samplesthres$days[matchday]
    match1 <- match(samples,datalist[[vec[1]]][,1])
    match2 <- match(samples,datalist[[vec[2]]][,1])
    data <- data.frame(samples,datalist[[vec[1]]][match1,-1],datalist[[vec[2]]][match2,-1])
  }
  else {
    inisamples <- intersect(datalist[[vec[1]]][,1],datalist[[vec[2]]][,1])
    for(i in 3:length(vec)){
      temp <- intersect(datalist[[vec[i]]][,1],inisamples)
      inisamples <- temp
    }
    samples <- inisamples
    matchday <- match(samples,samplesthres$samplecode)
    days <- samplesthres$days[matchday]
    matchlist <- matrix(0,length(samples),length(vec))
    for(i in 1:length(vec))
      matchlist[,i] <- match(samples,datalist[[vec[i]]][,1])
    data <- datalist[[vec[1]]][matchlist[,1],]
    for(i in 2:length(vec)){
      data <- data.frame(data,datalist[[vec[i]]][matchlist[,i],-1])
    }
  }
  combdata <- data.frame(data,days)
  combdata <- combdata[match(commonsamples,combdata[,1]),]
  combdata
}

datacombiner2 <- function(vec,type,datatypes,featureslist) {
  if(type == 1) 
    datalist <- emtdatalist(datatypes,featureslist)
  if(type == 2)
    datalist <- randomdatalist()
  if(type == 3)
    datalist <- simudatalist()
  data <- NULL
  days <- NULL
  #return data and days
  # replace samplesthres$samplecode, with common samples code
  days <- samplesthres$days[match(commonsamples,samplesthres$samplecode)]
  count <- 0
  for(i in vec){
    count <- count + 1
    if(count > 1){
      data <- data.frame(data,datalist[[i]][match(commonsamples,datalist[[i]][,1]),-1])
    }else{
      data <- data.frame(datalist[[i]][match(commonsamples,datalist[[i]][,1]),-1])
    }
  }
  data.frame(commonsamples,data,days)
}

rulegenerator <- function(t) {
  # ordert <- t[order(t$days),]
  # n=floor(nrow(t)/3)
  # sampleidx <- c(1:n,(nrow(ordert)-(n-1)):nrow(ordert))
  classvar <- factor(t$days >= thresupper,labels = c("long","short"))
  newdata <- data.frame(sapply(t[,-c(1,ncol(t))],function(x) factor(x,labels=c("high","low"))),classvar)
  transctns <- as(newdata,"transactions")
  rules <- apriori(transctns,appearance = list(rhs = c("classvar=long", "classvar=short"),
                                               default="lhs"),parameter = list(confidence=0.8,support = 0.1,target="rules",maxtime=10))
  # rules.pruned <- NULL
  # if(nrow(rules@quality)==0)
  #   nrules <- 0
  # else{
  #   rules.sorted <- sort(rules, by="support")
  #   #inspect(rules.sorted)
  #   subset.matrix <- is.subset(rules.sorted, rules.sorted)
  #   subset.matrix[lower.tri(subset.matrix, diag=T)] <- NA
  #   redundant <- colSums(subset.matrix, na.rm=T) >= 1
  #   rules.pruned <- rules.sorted[!redundant]
  #   # inspect(rules.pruned)
  #   nrules <- nrow(rules.pruned@quality)
  # }
  # rules.pruned
  rules <- rules[!is.redundant(rules)]
  rules <- sort(rules,by="confidence")
  list(rules,nrow(rules@quality),nrow(t))
}

ruletesterPCL <- function(inrules,test) {
  accuracy <- 0
  coverage <- 0
  if(length(inrules)==0){
    accuracy <- 0
    coverage <- 0
  }else{
    lhsruletable <- as.matrix(inrules@lhs@data)
    rhsruletable <- as.matrix(inrules@rhs@data)
    nr <- nrow(lhsruletable)
    nc <- ncol(lhsruletable)
    # lhsruletable[c(nr-1,nr),] <- rhsruletable[c(nr-1,nr),]
    # TRUE: long # high: odd number TRUE, low: even number TRUE
    # read in dataset
    # transform to tranction and then matrices, if a line match, count the class label.
    test$days <- factor(test$days >= thresupper,labels = c("long","short"))
    newdata <- data.frame(sapply(test[,-c(1,ncol(test))],function(x) factor(x,labels=c("high","low"))),test$days)
    transctns <- as(newdata,"transactions")
    testdata <- as.matrix(transctns@data)
    #find the number of rule matches for each sample
    truelabel <- testdata[nr,]
    prediction <- rep(0,ncol(testdata))
    for(i in 1:ncol(testdata)){
      predictionbysinglerules <- rep(-1,ncol(lhsruletable))
      for(j in 1:ncol(lhsruletable)){
        idx <- which(lhsruletable[,j]==TRUE)
        temp <- sum(testdata[idx,i])
        # print(c(sum(lhsruletable[,j])-1,temp))
        if(temp==length(idx)){
          predictionbysinglerules[j] <- rhsruletable[nr,j]
        }
      }
      # numtrue <- sum(predictionbysinglerules==TRUE)
      # numFALSE <- sum(predictionbysinglerules==FALSE)
      Tidx <- predictionbysinglerules==TRUE
      Fidx <- predictionbysinglerules==FALSE
      
      Trules <- which(rhsruletable[nr,]==TRUE)
      Frules <- which(rhsruletable[nr,]==FALSE)
      
      sumTidx <- if(sum(Tidx)>10) 10 else sum(Tidx)
      sumFidx <- if(sum(Fidx)>10) 10 else sum(Fidx)
      
      Tvec <- inrules@quality$support[Trules][order(inrules@quality$support[Trules],decreasing = TRUE)][1:sumTidx]
      Fvec <- inrules@quality$support[Frules][order(inrules@quality$support[Frules],decreasing = TRUE)][1:sumFidx]
      
      if(sumTidx<=0){
        pcltrue <- 0
      }else{
        temp <- inrules@quality$support[Tidx]
        Tvec1 <- temp[order(temp,decreasing = TRUE)][1:sumTidx]
        pcltrue <- sum(Tvec1/Tvec)/sumTidx
      }
      
      if(sumFidx<=0){
        pclfalse <- 0
      }else{
        temp <- inrules@quality$support[Fidx]
        Fvec1 <- temp[order(temp,decreasing = TRUE)][1:sumFidx]
        pclfalse <- sum(Fvec1/Fvec)/sumFidx
      }
      #decide the label
      label <- ifelse(pcltrue>=pclfalse, TRUE, FALSE)
      #store the prediction
      prediction[i] <- label
      if(sum(Tidx)==0 && sum(Fidx)==0)
        prediction[i] <- -1
      print(c(pcltrue,pclfalse,sum(Tidx),sum(Fidx),truelabel[i],prediction[i]==truelabel[i]))
    }
    accuracy <- sum(prediction==truelabel)/sum(prediction>=0)
    coverage <- sum(prediction>=0)/ncol(testdata)
  }
  c(accuracy,coverage)
}

ruletester <- function(inrules,data) {
  accuracy <- 0
  coverage <- 0
  if(length(inrules)==0){
    accuracy <- 0
    coverage <- 0
  }else{
    lhsruletable <- as.matrix(inrules@lhs@data)
    rhsruletable <- as.matrix(inrules@rhs@data)
    nr <- nrow(lhsruletable)
    nc <- ncol(lhsruletable)
    lhsruletable[c(nr-1,nr),] <- rhsruletable[c(nr-1,nr),]
    test$days <- factor(test$days >= thresupper,labels = c("long","short"))
    newdata <- data.frame(sapply(test[,-c(1,ncol(test))],function(x) factor(x,labels=c("high","low"))),test$days)
    transctns <- as(newdata,"transactions")
    testdata <- as.matrix(transctns@data)
    #find the number of rule matches for each sample
    truelabel <- testdata[nr,]
    prediction <- rep(0,ncol(testdata))
    for(i in 1:ncol(testdata)){
      predictionbysinglerules <- rep(-1,ncol(lhsruletable))
      for(j in 1:ncol(lhsruletable)){
        idx <- which(lhsruletable[,j]==TRUE)
        idx <- idx[-length(idx)]
        temp <- sum(testdata[idx,i])
        # print(c(sum(lhsruletable[,j])-1,temp))
        if(temp==length(idx)){
          predictionbysinglerules[j] <- lhsruletable[nr,j]
        }
      }
      numtrue <- sum(predictionbysinglerules==TRUE)
      numFALSE <- sum(predictionbysinglerules==FALSE)
      #decide the label
      label <- -5
      if(numtrue==0 && numFALSE==0){
        label <- -5
      }else{
        label <- ifelse(numtrue>=numFALSE, TRUE, FALSE)
      }
      #store the prediction
      prediction[i] <- label
      print(c(numtrue,numFALSE,truelabel[i],prediction[i]==truelabel[i]))
    }
    accuracy <- sum(prediction==truelabel)/sum(prediction>=0)
    coverage <- sum(prediction>=0)/ncol(testdata)
  }
  c(accuracy,coverage)
}

rulepredictor <- function(vec,type,rulelength) {
  data <- datacombiner(vec,type)
  y <- data$days >= thresupper
  accurule <- rep(0,50)
  coverage <- rep(0,50)
  for(rep in 1:10){
    flds <- createFolds(y, k = 5, list = TRUE, returnTrain = TRUE)
    i=(rep-1)*5
    for(f in 1:length(flds)){
      i <- i+1
      trainindex <- flds[[f]]
      train <- data[trainindex,]
      test <- data[-trainindex,]
      inrules <- rulegenerator(train,rulelength)[[1]]
      result <- ruletesterPCL(inrules,test)
      accurule[i] <- result[1]
      coverage[i] <- result[2]
    }
  }
  list(accurule,coverage)
}

numofrules <- function (type,datatypes, featureslist){
  resultrules <- list()
  dt <- character()
  # l <- length(rulelengthvec)
  nrule <- rep(0,7)
  nsample <- rep(0,7)
  checknumfeatures <- matrix(0,7)
  # for(p in 1:l){
  #   rulelength <- rulelengthvec[p]
    featureslist2 <- featureslist
    count <- 0
    for(i in 1:3){
      comb <- combn(3,i)
      for(j in 1:ncol(comb)){
        count <- count + 1
        print(count)
        dt[count] <- paste(comb[,j],collapse = "")
        if(length(comb[,j])==2){
          featureslist2 <- list(featureslist[[1]][1:10],featureslist[[2]][1:10],featureslist[[3]][1:10])
        }
        if(length(comb[,j])==3){
          featureslist2 <- list(featureslist[[1]][1:7],featureslist[[2]][1:7],featureslist[[3]][1:7])
        }
        # temp <- rulegenerator(datacombiner(comb[,j],type, datatypes, featureslist2))
        temp <- rulegenerator(datacombiner(comb[,j],type, datatypes, featureslist2))
        resultrules <- c(resultrules,temp[[1]])
        nrule[count] <- temp[[2]]
        nsample[count] <- temp[[3]]
        checknumfeatures[count] <- sum(sapply(featureslist2, length))
      # }
    }
  }
  rulecount <- data.frame(dt,nrule,nsample)
  names(rulecount) <- c("DataComb","numrules","NSample")
  list(rulecount, resultrules)
}

rulestatistics <- function(ruleset){
  #number of rules
  totalnumofrules <- length(ruleset)
  #number of rules for long and short
  numofruleslong <- unname(table(ruleset@rhs@data@i)[1])
  numofrulesshort <- unname(table(ruleset@rhs@data@i)[2])
  #number of rules of different length, for long survival group
  longindex <- ((ruleset@rhs@data@i) == (nrow(ruleset@rhs@data)-2))
  shortindex <- ((ruleset@rhs@data@i) == (nrow(ruleset@rhs@data)-1))
  
  if(sum(longindex) < 2)
    rulelengthlong <- sum(ruleset@lhs@data[,longindex]) else
      rulelengthlong <- table(colSums((ruleset@lhs@data[,longindex])))
  
  # rulelengthlong <- table(colSums((ruleset@lhs@data[,longindex])))
  #for short survival group
  if(sum(shortindex) < 2)
    rulelengthshort <- sum(ruleset@lhs@data[,shortindex]) else
      rulelengthshort <- table(colSums((ruleset@lhs@data[,shortindex])))
  
  if(totalnumofrules < 50)
    top <- totalnumofrules else
      top <- 50
  avgvalues <- colMeans(ruleset@quality[1:top,])
  avgrulelength<- mean(colSums((ruleset@lhs@data)))
  
  list(numofruleslong,numofrulesshort,rulelengthlong,rulelengthshort,avgvalues,avgrulelength)
}

datatypes <- datatypes74
results_rules74 <- foreach (i = 1:11) %dopar% {
  featureslist <- list(level1_74_features[[i]],level2_74_features[[i]],level3_74_features[[i]])
  numofrules(1,datatypes,featureslist)
}
numrules74 <- as.data.frame(sapply(results_rules74, function(x) x[[1]][,2]),row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(numrules74) <- c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble")

datatypes <- datatypes123
results_rules123 <- foreach (i = 1:11) %dopar% {
  featureslist <- list(level1_123_features[[i]],level2_123_features[[i]],level3_123_features[[i]])
  numofrules(1,datatypes,featureslist)
}
numrules123 <- as.data.frame(sapply(results_rules123, function(x) x[[1]][,2]),row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(numrules123) <- c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble")

datatypes <- datatypes455
results_rules455 <- foreach (i = 1:11) %dopar% {
  featureslist <- list(level1_455_features[[i]],level2_455_features[[i]],level3_455_features[[i]])
  numofrules(1,datatypes,featureslist)
}
numrules455 <- as.data.frame(sapply(results_rules455, function(x) x[[1]][,2]),row.names = c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA"))
names(numrules455) <- c("t-test","Lasso", "NetLasso", "addDA2","Netrank","stSVM","Cox","RegCox","MSS","Survnet","Ensemble")

# save.image('LUAD_with_rules.RData')

#find representative rules

# inspect(head(sort(results_rules74[[1]][[2]][[7]],by="confidence"),50))
# inspect(sort(results_rules74[[11]][[2]][[7]],by="confidence"))
# 
# toprulegenes <- unique(c("DNMT1","PAK1","MDM2","CDH3","miR-34a","HMGA2","LOXL2","ETS1","LOXL2","KRAS","AKT1","KRAS","KRT18","miR-101-1","AKT1","MAP2K1",
#                          "EZH2","FBXL14","BTRC"))
# temp <- subgraph(LEnet74,c("SNAI1","RAF1","ZEB2","ZEB1","CDH1","GSK3B","NFKB1","TWIST1","BMI1",toprulegenes))
# plot(temp)

# distribution of rules
ratiosmatrix <- matrix(0,11,7)
supportmatrix <- matrix(0,11,7)
confidencematrix <- matrix(0,11,7)
liftmatrix <- matrix(0,11,7)
rulelengthmatrix <- matrix(0,11,7)
#one graph for the number of rules
for(i in 1:11){
  rules7comb <- results_rules74[[i]][[2]]
  temp <- lapply(rules7comb, function(x) rulestatistics(x))
  # for(j in 1:7){
  #   temp2 <- bind_rows(temp[[j]][[3]],temp[[j]][[4]])
  #   for(m in 1:nrow(temp2)){
  #     for(n in 1:ncol(temp2)){
  #       if(is.na(temp2[m,n])==TRUE)
  #         temp2[m,n] <- 0
  #     }
  #   }
  #   barplot(as.matrix(as.data.frame(temp2)))
  # }
  ratiovec <- rep(0,7)
  for(j in 1:7){
    ratiosmatrix[i,j] <- temp[[j]][[1]]/temp[[j]][[2]]
    supportmatrix[i,j] <- unname(temp[[j]][[5]][1])
    confidencematrix[i,j] <- unname(temp[[j]][[5]][2])
    liftmatrix[i,j] <- unname(temp[[j]][[5]][3])
    rulelengthmatrix[i,j] <- temp[[j]][[6]]
  }
  # plot(1:7,ratiovec, type='l')
}
#replace NA with max(row)
for(i in 1:11){
  idx <- which(is.na(ratiosmatrix[i,]))
  ratiosmatrix[i,idx] <- max(ratiosmatrix[i,],na.rm = TRUE)
}
#scatter plot of the ratios
# ratiospoint <- NULL
# for(i in 1:11){
#   for(j in 1:7){
#     ratiospoint <- rbind(ratiospoint, c(ratiosmatrix[i,j],j))
#   }
# }

# plot of the number of rules
ggdata <- data.frame(t(numrules74)[-c(4,6,10),])
# ggdata <- data.frame(t(numrules74))
names(ggdata) <- c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA")
ggdata <- ggdata %>% gather(datacomb,value,1:7,factor_key=TRUE)
p1 <- ggplot(ggdata, aes(y=value,x=datacomb)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Number of rules") +
  scale_x_discrete(name = "Data Combinations") +
  theme_bw()+
  theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),
        text = element_text(size = 14, family = "Times"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12),
        legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) 

# plot of the ratio of rules
ggdata <- data.frame(ratiosmatrix[-c(4,6,10),])
# ggdata <- data.frame(ratiosmatrix)
names(ggdata) <- c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA")
ggdata <- ggdata %>% gather(datacomb,value,1:7,factor_key=TRUE)
p2 <- ggplot(ggdata, aes(y=value,x=datacomb)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Ratios between two classes") +
  scale_x_discrete(name = "Data Combinations") +
  theme_bw()+
  theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),
        text = element_text(size = 14, family = "Times"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12),
        legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) 

#plot of average confidence of top 50 rules

ggdata <- data.frame(confidencematrix[-c(4,6,10),])
# ggdata <- data.frame(confidencematrix)
names(ggdata) <- c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA")
ggdata <- ggdata %>% gather(datacomb,value,1:7,factor_key=TRUE)
p3 <- ggplot(ggdata, aes(y=value,x=datacomb)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Average confidence") +
  scale_x_discrete(name = "Data Combinations") +
  theme_bw()+
  theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),
        text = element_text(size = 14, family = "Times"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12),
        legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) 

# ggdata <- data.frame(supportmatrix)
# names(ggdata) <- c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA")
# ggdata <- ggdata %>% gather(datacomb,value,1:7,factor_key=TRUE)
# p4 <- ggplot(ggdata, aes(y=value,x=datacomb)) +
#   geom_boxplot(alpha=0.7) +
#   scale_y_continuous(name = "Average Support of rules") +
#   scale_x_discrete(name = "Data Combinations") +
#   theme_bw()+
#   theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),
#         text = element_text(size = 14, family = "Times"),
#         axis.title = element_text(face="bold"),
#         axis.text.x=element_text(size = 12),
#         legend.position = "top") + 
#   theme(axis.text.x = element_text(angle = 30, hjust = 1)) 

ggdata <- data.frame(liftmatrix[-c(4,6,10),])
# ggdata <- data.frame(liftmatrix)
names(ggdata) <- c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA")
ggdata <- ggdata %>% gather(datacomb,value,1:7,factor_key=TRUE)
p4 <- ggplot(ggdata, aes(y=value,x=datacomb)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Average lift") +
  scale_x_discrete(name = "Data Combinations") +
  theme_bw()+
  theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),
        text = element_text(size = 14, family = "Times"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12),
        legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) 

#rulelength
ggdata <- data.frame(rulelengthmatrix)
ggdata <- data.frame(rulelengthmatrix[-c(4,6,10),])
names(ggdata) <- c("GE","DM","CNA","GE+DM","GE+CNA","DM+CNA","GE+DM+CNA")
ggdata <- ggdata %>% gather(datacomb,value,1:7,factor_key=TRUE)
p5 <- ggplot(ggdata, aes(y=value,x=datacomb)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Average length") +
  scale_x_discrete(name = "Data Combinations") +
  theme_bw()+
  theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),
        text = element_text(size = 14, family = "Times"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12),
        legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) 


grid.arrange(p1,p5,p3,p4,ncol=2)

#the next good example:  inspect(head(sort(results_rules123[[11]][[2]][[7]],by="confidence"),60))

# {GATA6=high,CDC20=low,FOXA3=high}                          => {classvar=short} 0.1503759 1.0000000  2.046154
# [6]  {GATA6=high,CDH3=high,FOXA3=high}                          => {classvar=short} 0.1353383 1.0000000  2.046154
# [7]  {GATA6=high,miR.34a.1=low,FOXA3=high}                      => {classvar=short} 0.1428571 1.0000000  2.046154
# [8]  {miR.34a=low,CCND1=high,FOXA3=low}                         => {classvar=long}  0.1127820 1.0000000  1.955882
# [9]  {GATA6=low,miR.34a=low,CCND1=high}                         => {classvar=long}  0.1052632 1.0000000  1.955882
# [10] {BIRC3=low,CCND1=high,FOXA3=low}                           => {classvar=long}  0.1127820 1.0000000  1.955882

#get names from the rules
templeft <- as(lhs(head(sort(results_rules123[[11]][[2]][[7]],by="confidence"),20)),'list')
tempright <-  as(rhs(head(sort(results_rules123[[11]][[2]][[7]],by="confidence"),20)),'list')
# gsub("\\{|\\}", "", temp$lhs)
mollist <- c()
for(i in 1:length(templeft)){
  mollist <- c(mollist, gsub("\\=.*","",templeft[[i]]))
}
mollist <- sort(unique(mollist))[1:16]
mollist[12] <- 'HMGA2'
mollist <- gsub('\\.','-',mollist)
plot(subgraph(LEnet455,c(mollist,"SNAI1","CDH1","VIM","OCLN","miR-200a","miR-130b","ZEB2")))



