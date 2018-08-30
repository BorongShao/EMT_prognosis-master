cpy <- LUAD455_methy_700_1400
lengthcpy <- length(cpy)

# all features
auc_allemt455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[1]][[1]])
aupr_allemt455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[1]][[2]])
accu_allemt455 <- rowSums(sapply(1:length(cpy), function(x) cpy[[x]][[1]][[3]]))/lengthcpy
pvalue_allemt455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[1]][[4]])
timeauc_allemt455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[1]][[5]]),na.rm=TRUE)/lengthcpy

# clinical features 
# auc_clinical <- sapply(1:lengthcpy, function(x) LUAD_folds_clinical_700_1400[[x]][[1]][[1]])
# aupr_clinical <- sapply(1:lengthcpy, function(x) LUAD_folds_clinical_700_1400[[x]][[1]][[2]])
# accu_clinical <- rowSums(sapply(1:lengthcpy, function(x) LUAD_folds_clinical_700_1400[[x]][[1]][[3]]))/lengthcpy
# 
# auc_clinical_rnf <- sapply(1:lengthcpy, function(x) LUAD_folds_clinical_700_1400[[x]][[2]][[1]])
# aupr_clinical_rnf <- sapply(1:lengthcpy, function(x) LUAD_folds_clinical_700_1400[[x]][[2]][[2]])
# accu_clinical_rnf <- rowSums(sapply(1:lengthcpy, function(x) LUAD_folds_clinical_700_1400[[x]][[2]][[3]]))/lengthcpy

# ttest
auc_ttest455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[3]][[1]])
aupr_ttest455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[3]][[2]])
accu_ttest455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[3]][[3]]))/lengthcpy
pvalue_ttest455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[3]][[4]])
timeauc_ttest455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[3]][[5]]),na.rm=TRUE)/lengthcpy

#lasso
auc_lasso455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[5]][[1]])
aupr_lasso455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[5]][[2]])
accu_lasso455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[5]][[3]]))/lengthcpy
pvalue_lasso455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[5]][[4]])
timeauc_lasso455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[5]][[5]]),na.rm=TRUE)/lengthcpy

#netlasso
auc_netlasso455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[7]][[1]])
aupr_netlasso455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[7]][[2]])
accu_netlasso455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[7]][[3]]))/lengthcpy
pvalue_netlasso455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[7]][[4]])
timeauc_netlasso455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[7]][[5]]),na.rm=TRUE)/lengthcpy

#add
auc_stSVM455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[9]][[1]])
aupr_stSVM455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[9]][[2]])
accu_stSVM455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[9]][[3]]))/lengthcpy
pvalue_stSVM455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[9]][[4]])
timeauc_stSVM455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[9]][[5]]),na.rm=TRUE)/lengthcpy

#netrank
auc_add455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[11]][[1]])
aupr_add455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[11]][[2]])
accu_add455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[11]][[3]]))/lengthcpy
pvalue_add455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[11]][[4]])
timeauc_add455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[11]][[5]]),na.rm=TRUE)/lengthcpy

#stSVM
auc_netrank455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[13]][[1]])
aupr_netrank455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[13]][[2]])
accu_netrank455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[13]][[3]]))/lengthcpy
pvalue_netrank455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[13]][[4]])
timeauc_netrank455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[13]][[5]]),na.rm=TRUE)/lengthcpy

#cox
auc_cox455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[15]][[1]])
aupr_cox455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[15]][[2]])
accu_cox455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[15]][[3]]))/lengthcpy
pvalue_cox455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[15]][[4]])
timeauc_cox455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[15]][[5]]),na.rm=TRUE)/lengthcpy

#cox_lasso
auc_coxreg455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[17]][[1]])
aupr_coxreg455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[17]][[2]])
accu_coxreg455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[17]][[3]]))/lengthcpy
pvalue_coxreg455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[17]][[4]])
timeauc_coxreg455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[17]][[5]]),na.rm=TRUE)/lengthcpy

#RDS
auc_rds455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[19]][[1]])
aupr_rds455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[19]][[2]])
accu_rds455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[19]][[3]]))/lengthcpy
pvalue_rds455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[19]][[4]])
timeauc_rds455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[19]][[5]]),na.rm=TRUE)/lengthcpy

#survnet
auc_survnet455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[21]][[1]])
aupr_survnet455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[21]][[2]])
accu_survnet455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[21]][[3]]))/lengthcpy
pvalue_survnet455 <- sapply(1:lengthcpy, function(x) cpy[[x]][[21]][[4]])
timeauc_survnet455 <- rowSums(sapply(1:lengthcpy, function(x) cpy[[x]][[21]][[5]]),na.rm=TRUE)/lengthcpy



list(threstrainidx,svm_clittest2,svm_clilasso2,svm_clinetlasso2,svm_clistsvmfs2,svm_cliaddg2,svm_clinetrank2,svm_clicoxfs2,svm_clicoxregfs2,
     svm_clirdsfs2,svm_clisurvnetfs2,rnf_clittest2,rnf_clilasso2,rnf_clinetlasso2,rnf_clistsvmfs2,rnf_cliaddg2,rnf_clinetrank2,rnf_clicoxfs2,rnf_clicoxregfs2,
     rnf_clirdsfs2,rnf_clisurvnetfs2)

cpy <- LUAD455_FSFmethy_700_1400
lengthcpy <- length(cpy)
# auc_random_allfeatures <- sapply(1:lengthcpy, function(x) cpy[[x]][[1]][[1]])
# auc_random_allfeatures_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[2]][[1]])
auc_random_ttest <- sapply(1:lengthcpy, function(x) cpy[[x]][[2]][[1]])
auc_random_ttest_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[12]][[1]])
auc_random_lasso <- sapply(1:lengthcpy, function(x) cpy[[x]][[3]][[1]])
auc_random_lasso_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[13]][[1]])
auc_random_netlasso <- sapply(1:lengthcpy, function(x) cpy[[x]][[4]][[1]])
auc_random_netlasso_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[14]][[1]])
auc_random_add <- sapply(1:lengthcpy, function(x) cpy[[x]][[6]][[1]])
auc_random_add_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[16]][[1]])
auc_random_netrank <- sapply(1:lengthcpy, function(x) cpy[[x]][[7]][[1]])
auc_random_netrank_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[17]][[1]])

# aupr_random_allfeatures <- sapply(1:lengthcpy, function(x) cpy[[x]][[26]][[2]])
# aupr_random_allfeatures_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[36]][[2]])
aupr_random_ttest <- sapply(1:lengthcpy, function(x) cpy[[x]][[2]][[2]])
aupr_random_ttest_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[12]][[2]])
aupr_random_lasso <- sapply(1:lengthcpy, function(x) cpy[[x]][[3]][[2]])
aupr_random_lasso_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[13]][[2]])
aupr_random_netlasso <- sapply(1:lengthcpy, function(x) cpy[[x]][[4]][[2]])
aupr_random_netlasso_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[14]][[2]])
aupr_random_add <- sapply(1:lengthcpy, function(x) cpy[[x]][[6]][[2]])
aupr_random_add_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[16]][[2]])
aupr_random_netrank <- sapply(1:lengthcpy, function(x) cpy[[x]][[7]][[2]])
aupr_random_netrank_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[17]][[2]])

# accu_random_allfeatures <- sapply(1:lengthcpy, function(x) cpy[[x]][[1]][[3]][51])
# accu_random_allfeatures_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[2]][[3]][51])
accu_random_ttest <- sapply(1:lengthcpy, function(x) cpy[[x]][[2]][[3]][51])
accu_random_ttest_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[12]][[3]][51])
accu_random_lasso <- sapply(1:lengthcpy, function(x) cpy[[x]][[3]][[3]][51])
accu_random_lasso_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[13]][[3]][51])
accu_random_netlasso <- sapply(1:lengthcpy, function(x) cpy[[x]][[4]][[3]][51])
accu_random_netlasso_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[14]][[3]][51])
accu_random_add <- sapply(1:lengthcpy, function(x) cpy[[x]][[6]][[3]][51])
accu_random_add_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[16]][[3]][51])
accu_random_netrank <- sapply(1:lengthcpy, function(x) cpy[[x]][[7]][[3]][51])
accu_random_netrank_rnf <- sapply(1:lengthcpy, function(x) cpy[[x]][[17]][[3]][51])

#calculate p values
# cpy <- LUAD455_methy_700_1400
# rcpy <- LUAD455_climethy_700_1400
# matchedindex <- rep(0,length(cpy))
# for(i in 1:length(cpy)){
#   indexi <- cpy[[i]][[33]]
#   for(j in 1:length(rcpy)){
#     if(all.equal(indexi,rcpy[[j]][[23]])==TRUE){
#       matchedindex[i] <- j
#       break
#     }
#   }
# }

matchedindex <- sample(300,length(LUAD455_methy_700_1400))

# pvalue_allemt_auc <- t.test(auc_allemt455,auc_random_allfeatures[matchedindex],paired = TRUE)$p.value
# pvalue_allemt_aupr <- t.test(auc_allemt455,aupr_random_allfeatures[matchedindex],paired = TRUE)$p.value
# pvalue_allemt_accu <- t.test(dt2_accu_allemt455,accu_random_allfeatures[matchedindex],paired = TRUE)$p.value

pvalue_ttest_auc <- t.test(auc_ttest455,auc_random_ttest[matchedindex],paired = TRUE)$p.value
pvalue_ttest_aupr <- t.test(auc_ttest455,aupr_random_ttest[matchedindex],paired = TRUE)$p.value
pvalue_ttest_accu <- t.test(dt2_accu_ttest455,accu_random_ttest[matchedindex],paired = TRUE)$p.value

pvalue_lasso_auc <- t.test(auc_lasso455,auc_random_lasso[matchedindex],paired = TRUE)$p.value
pvalue_lasso_aupr <- t.test(auc_lasso455,aupr_random_lasso[matchedindex],paired = TRUE)$p.value
pvalue_lasso_accu <- t.test(dt2_accu_lasso455,accu_random_lasso[matchedindex],paired = TRUE)$p.value

pvalue_netlasso_auc <- t.test(auc_netlasso455,auc_random_netlasso[matchedindex],paired = TRUE)$p.value
pvalue_netlasso_aupr <- t.test(auc_netlasso455,aupr_random_netlasso[matchedindex],paired = TRUE)$p.value
pvalue_netlasso_accu <- t.test(dt2_accu_netlasso455,accu_random_netlasso[matchedindex],paired = TRUE)$p.value

pvalue_add_auc <- t.test(auc_add455,auc_random_add[matchedindex],paired = TRUE)$p.value
pvalue_add_aupr <- t.test(auc_add455,aupr_random_add[matchedindex],paired = TRUE)$p.value
pvalue_add_accu <- t.test(dt2_accu_add455,accu_random_add[matchedindex],paired = TRUE)$p.value

pvalue_netrank_auc <- t.test(auc_netrank455,auc_random_netrank[matchedindex],paired = TRUE)$p.value
pvalue_netrank_aupr <- t.test(auc_netrank455,aupr_random_netrank[matchedindex],paired = TRUE)$p.value
pvalue_netrank_accu <- t.test(dt2_accu_netrank455,accu_random_netrank[matchedindex],paired = TRUE)$p.value

pvaluedata <- c(pvalue_ttest_auc,pvalue_ttest_aupr,pvalue_ttest_accu,
                pvalue_lasso_auc,pvalue_lasso_aupr,pvalue_lasso_accu,pvalue_netlasso_auc,pvalue_netlasso_aupr,pvalue_netlasso_accu,
                pvalue_add_auc,pvalue_add_aupr,pvalue_add_accu,pvalue_netrank_auc,pvalue_netrank_aupr,pvalue_netrank_accu)
pvaluecharac <- data.frame(sapply(pvaluedata, function(x) paste("P-value =\n", format(x,digits=3))))
temp1 <- factor(rep(c(1,2,3),5),labels=c("AUC Value","AUPR Value","Accuracy at 0.5 Cutoff"))
temp2 <- factor(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3)),
                labels = c("t-test","Lasso", "NetLasso", "addDA2","Netrank"))
pvaluedataframe <- data.frame(temp1,temp2,pvaluecharac)
names(pvaluedataframe) <- c("metrics","Algorithms","pvalues")


nrun1 <- length(auc_ttest455)
nrun2 <- length(auc_random_ttest)
ngroup=2
Measurements <- c(auc_ttest455,auc_random_ttest,auc_lasso455,auc_random_lasso,auc_netlasso455,auc_random_netlasso,
                  auc_add455,auc_random_add,auc_netrank455,auc_random_netrank,aupr_ttest455,aupr_random_ttest,
                  aupr_lasso455,aupr_random_lasso,aupr_netlasso455,aupr_random_netlasso,aupr_add455,aupr_random_add,aupr_netrank455,aupr_random_netrank,
                  dt2_accu_ttest455,accu_random_ttest,dt2_accu_lasso455,accu_random_lasso,dt2_accu_netlasso455,
                  accu_random_netlasso,dt2_accu_add455,accu_random_add,dt2_accu_netrank455,accu_random_netrank)
Features <- factor(rep(c(rep(1,nrun1),rep(2,nrun2)),15),
                   labels = c("Individually selected features","Frequently selected features"))
Algorithms <- factor(rep(c(rep(1,nrun1+nrun2),rep(2,nrun1+nrun2),rep(3,nrun1+nrun2),rep(4,nrun1+nrun2),rep(5,nrun1+nrun2)),3),
                     labels = c("t-test","Lasso", "NetLasso", "addDA2","Netrank"))
metrics <- factor(c(rep(1,(nrun1+nrun2)*5),rep(2,(nrun1+nrun2)*5),rep(3,(nrun1+nrun2)*5)),labels = c("AUC Value","AUPR Value","Accuracy at 0.5 Cutoff"))
ggdata <- data.frame(Measurements,Features,Algorithms,metrics)
# ann_text <- data.frame(mpg = 15,wt = 5,lab = "Text",
#                        cyl = factor(8,levels = c("4","6","8")))
# p + geom_text(data = ann_text,label = "Text")
pdf('../Thesis_Results/alllevels/random_methy455_svm.pdf', width = 8.5,height = 11)
ggplot(ggdata) + 
  geom_density(aes(x=Measurements,y=..scaled..,fill=Features),alpha=0.15) + 
  facet_grid(Algorithms ~ metrics) +
  scale_y_continuous(name = "Density") +
  theme_bw()+
  theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),
        text = element_text(size = 14, family = "Times"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12),
        legend.position = "top") + 
  geom_label(aes(x=0.33,y=0.88,label=pvalues), size=3,
             data=pvaluedataframe)
# theme(axis.title.x=element_blank())
dev.off()

#AUC boxplot
nrun1 <- length(auc_ttest455)
nrun2 <- length(auc_random_ttest)
nmethod <- 2
ngroup <- 5
results <- c(auc_ttest455,auc_lasso455,auc_netlasso455,auc_add455,auc_netrank455,
             auc_random_ttest,auc_random_lasso,auc_random_netlasso,auc_random_add,auc_random_netrank)
results2 <- c(aupr_ttest455,aupr_lasso455,aupr_netlasso455,aupr_add455,aupr_netrank455,
             aupr_random_ttest,aupr_random_lasso,aupr_random_netlasso,aupr_random_add,aupr_random_netrank)

Features <- factor(c(rep(1,nrun1*ngroup),rep(2,nrun2*ngroup),rep(1,nrun1*ngroup),rep(2,nrun2*ngroup)),
                labels = c("Individually selected features","Frequently selected features"))
`FS Algorithms` <- rep(factor(c(as.vector(sapply(1:5,function(x) rep(x,nrun1))),as.vector(sapply(1:5,function(x) rep(x,nrun2)))),
                          labels = c("t-test","Lasso", "NetLasso", "addDA2","Netrank")),2)
metric <- factor(c(rep(1,length(results)),rep(2,length(results))),labels=c("AUC", "AUPR"))

ggdata <- as.data.frame(cbind(c(results,results2),Features,`FS Algorithms`,metric))
ggdata$Features <- factor(group,labels =c("Individually selected features","Frequently selected features"))
ggdata$`FS Algorithms` <- factor(`FS Algorithms`,labels = c("t-test","Lasso", "NetLasso", "addDA2","Netrank"))
names(ggdata)[1] <- "results"
ggdata$metric <- factor(ggdata$metric, labels=c("AUC","AUPR"))
variancevalues <- sapply(list(auc_ttest455,auc_lasso455,auc_netlasso455,auc_add455,auc_netrank455,
                              auc_random_ttest,auc_random_lasso,auc_random_netlasso,auc_random_add,auc_random_netrank,
                              aupr_ttest455,aupr_lasso455,aupr_netlasso455,aupr_add455,aupr_netrank455,
                              aupr_random_ttest,aupr_random_lasso,aupr_random_netlasso,aupr_random_add,aupr_random_netrank),
                         function(x) var(x))

ggplot(ggdata, aes(x = `FS Algorithms`,
                   y = results, fill=Features)) +
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "AUC Values") +
  scale_x_discrete(name = "Feature Selection Algorithms") +
  theme_bw()+
  facet_grid(. ~ metric) + 
  theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),
        text = element_text(size = 14, family = "Times"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13),
        legend.position = "top") + 
   theme(axis.text.x = element_text(angle = 30, hjust = 1)) 





# feature plots

# concatenation, t-test 2, lasso 4, cox 14, mss 18, netlasso 6, netrank12, regcox 16.
FSFstatisticcont2dt <- function(algidx, numfeatures,objlist){
  cpy <- objlist[[1]]
  lengthcpy <- length(cpy)
  numfeatures2 <- numfeatures*2
  cox_sf <- sapply(1:lengthcpy, function(x) { 
    sf_74 <- rep(0,numfeatures2)
    sf_74[cpy[[x]][[algidx]]] <- 1
    sf_74
  })
  temp1 <-  order(rowSums(cox_sf),decreasing = TRUE)[1:20]
  
  cpy <- objlist[[2]]
  lengthcpy <- length(cpy)
  numfeatures <- numfeatures
  cox_sf <- sapply(1:lengthcpy, function(x) { 
    sf_74 <- rep(0,numfeatures)
    sf_74[cpy[[x]][[algidx]]] <- 1
    sf_74
  })
  temp2 <-  order(rowSums(cox_sf),decreasing = TRUE)[1:20]
  
  cpy <- objlist[[3]]
  lengthcpy <- length(cpy)
  numfeatures <- numfeatures
  cox_sf <- sapply(1:lengthcpy, function(x) { 
    sf_74 <- rep(0,numfeatures)
    sf_74[cpy[[x]][[algidx]]] <- 1
    sf_74
  })
  temp3 <- order(rowSums(cox_sf),decreasing = TRUE)[1:20]+numfeatures
  # ratio, and overlap
  c(sum(temp1>numfeatures)/20, sum(!is.na(match(temp1, temp2)))/20, sum(!is.na(match(temp1, temp3)))/20)
  # data.frame(rbind(temp2,temp3,temp1))
}

FSFstatisticcont3dt <- function(algidx, numfeatures,objlist){
  cpy <- objlist[[1]]
  lengthcpy <- length(cpy)
  numfeatures2 <- numfeatures*3
  cox_sf <- sapply(1:lengthcpy, function(x) { 
    sf_74 <- rep(0,numfeatures2)
    sf_74[cpy[[x]][[algidx]]] <- 1
    sf_74
  })
  temp1 <-  order(rowSums(cox_sf),decreasing = TRUE)[1:20]
  
  cpy <- objlist[[2]]
  lengthcpy <- length(cpy)
  numfeatures <- numfeatures
  cox_sf <- sapply(1:lengthcpy, function(x) { 
    sf_74 <- rep(0,numfeatures)
    sf_74[cpy[[x]][[algidx]]] <- 1
    sf_74
  })
  temp2 <-  order(rowSums(cox_sf),decreasing = TRUE)[1:20] 
  
  cpy <- objlist[[3]]
  lengthcpy <- length(cpy)
  numfeatures <- numfeatures
  cox_sf <- sapply(1:lengthcpy, function(x) { 
    sf_74 <- rep(0,numfeatures)
    sf_74[cpy[[x]][[algidx]]] <- 1
    sf_74
  })
  temp3 <-  order(rowSums(cox_sf),decreasing = TRUE)[1:20] + numfeatures
  
  cpy <- objlist[[4]]
  lengthcpy <- length(cpy)
  numfeatures <- numfeatures
  cox_sf <- sapply(1:lengthcpy, function(x) { 
    sf_74 <- rep(0,numfeatures)
    sf_74[cpy[[x]][[algidx]]] <- 1
    sf_74
  })
  temp4 <-  order(rowSums(cox_sf),decreasing = TRUE)[1:20] + numfeatures*2
  
  c(sum(temp1>numfeatures)/20, sum(!is.na(match(temp1, temp2)))/20, sum(!is.na(match(temp1, temp3)))/20,sum(!is.na(match(temp1, temp4)))/20)
}


conct74stats2dt <- sapply(c(2,4,14,18), function(x) FSFstatisticcont2dt(x,74,list(LUAD74_inte_level12_1,LUAD74ge,LUAD74dm)))
conct123stats2dt <- sapply(c(2,4,14,18), function(x) FSFstatisticcont2dt(x,123,list(LUAD123_inte_level12_1,LUAD123ge,LUAD123dm)))
conct455stats2dt <- sapply(c(2,4,14,18), function(x) FSFstatisticcont2dt(x,455,list(LUAD455_inte_level12_1,LUAD455ge,LUAD455dm)))

conct74stats3dt <- sapply(c(2,4,14,18), function(x) FSFstatisticcont3dt(x,74,list(LUAD74_inte_level123_1,LUAD74ge_123,LUAD74dm_123,LUAD74cna_123)))
conct123stats3dt <- sapply(c(2,4,14,18), function(x) FSFstatisticcont3dt(x,123,list(LUAD123_inte_level123_1,LUAD123ge_123,LUAD123dm_123,LUAD123cna_123)))
conct455stats3dt <- sapply(c(2,4,14,18), function(x) FSFstatisticcont3dt(x,455,list(LUAD455_inte_level123_1,LUAD455ge_123,LUAD455dm_123,LUAD455cna_123)))

plot(1:4,conct74stats2dt[1,], main="Milage vs. Car Weight", 
     xlab="Algorithms", ylab="Feature ratios (GE/non-GE)", pch=18, col="blue")
text(wt, mpg, row.names(mtcars), cex=0.6, pos=4, col="red")



candidates <- list(LUAD74_inte_level12_1,LUAD74_inte_level12_all,LUAD74_inte_level12_scoring1,LUAD74_inte_level12_scoring5)

# idxnetlasso <- findthebest(candidates,7)[[2]]
cpyfeatures <- LUAD74_gedm_cluster3_all
source('retrievefeatures.R')
netlasso74gedm <- netlassofs

# idxadd <- findthebest(candidates,11)[[2]]
cpyfeatures <- LUAD74_gedm_cluster3_1
source('retrievefeatures.R')
add74gedm <- addg_sf[1:10]

# idxnetrank <- findthebest(candidates,13)[[2]]
cpyfeatures <- LUAD74_gedm_cluster3_all
source('retrievefeatures.R')
netrank74gedm <- netrankfs

# idxcoxreg <- findthebest(candidates,17)[[2]]
cpyfeatures <- LUAD74_gedm_cluster3_scoring1
source('retrievefeatures.R')
coxreg74gedm <- coxregfs

# idxsurvnet <-findthebest(candidates,21)[[2]]
cpyfeatures <- LUAD74_gedm_cluster3_scoring5
source('retrievefeatures.R')
survnet74gedm <- survnet_sf[1:10]


cpyfeatures <- LUAD74ge_cluster3
source('retrievefeatures.R')
netlasso74ge <-  netlassofs
netrank74ge <- netrankfs
add74ge <- addg_sf[1:10]
coxreg74ge <- coxregfs
survnet74ge <- survnet_sf[1:10]

cpyfeatures <- LUAD74dm_cluster3
source('retrievefeatures.R')
netlasso74dm <-  netlassofs+74
netrank74dm <- netrankfs+74
add74dm <- sapply(addg_sf[1:10], function(x) x+74)
coxreg74dm <- coxregfs+74
survnet74dm <- sapply(survnet_sf[1:10], function(x) x+74)










candidates <-  list(LUAD123_inte_level12_1,LUAD123_inte_level12_2,LUAD123_inte_level12_all,LUAD123_inte_level12_scoring1,LUAD123_inte_level12_scoring5)

# idxnetlasso <- findthebest(candidates,7)[[2]]
cpyfeatures <- LUAD123_gedm_cluster3_all
source('retrievefeatures.R')
netlasso123gedm <- netlassofs

# idxadd <- findthebest(candidates,11)[[2]]
cpyfeatures <- LUAD123_gedm_cluster3_all
source('retrievefeatures.R')
add123gedm <- addg_sf[1:10]

# idxnetrank <- findthebest(candidates,13)[[2]]
cpyfeatures <- LUAD123_gedm_cluster3_1
source('retrievefeatures.R')
netrank123gedm <- netrankfs

# idxcoxreg <- findthebest(candidates,17)[[2]]
cpyfeatures <- LUAD123_gedm_cluster3_2
source('retrievefeatures.R')
coxreg123gedm <- coxregfs

# idxsurvnet <-findthebest(candidates,21)[[2]]
cpyfeatures <- LUAD123_gedm_cluster3_1
source('retrievefeatures.R')
survnet123gedm <- survnet_sf[1:10]


cpyfeatures <- LUAD123ge_cluster3
source('retrievefeatures.R')
netlasso123ge <-  netlassofs
netrank123ge <- netrankfs
add123ge <- addg_sf[1:10]
coxreg123ge <- coxregfs
survnet123ge <- survnet_sf[1:10]

cpyfeatures <- LUAD123dm_cluster3
source('retrievefeatures.R')
netlasso123dm <-  netlassofs+123
netrank123dm <- netrankfs+123
add123dm <- sapply(addg_sf[1:10], function(x) x+123)
coxreg123dm <- coxregfs+123
survnet123dm <- sapply(survnet_sf[1:10], function(x) x+123)







candidates <-  list(LUAD455_inte_level12_1,LUAD455_inte_level12_2,LUAD455_inte_level12_all,LUAD455_inte_level12_scoring1,LUAD455_inte_level12_scoring5)
# idxnetlasso <- findthebest(candidates,7)[[2]]
cpyfeatures <- LUAD455_gedm_cluster3_scoring1
source('retrievefeatures.R')
netlasso455gedm <- netlassofs

# idxadd <- findthebest(candidates,11)[[2]]
cpyfeatures <- LUAD455_gedm_cluster3_scoring1
source('retrievefeatures.R')
add455gedm <- addg_sf[1:10]

# idxnetrank <- findthebest(candidates,13)[[2]]
cpyfeatures <- LUAD455_gedm_cluster3_2
source('retrievefeatures.R')
netrank455gedm <- netrankfs

# idxcoxreg <- findthebest(candidates,17)[[2]]
cpyfeatures <- LUAD455_gedm_cluster3_all
source('retrievefeatures.R')
coxreg455gedm <- coxregfs

# idxsurvnet <-findthebest(candidates,21)[[2]]
cpyfeatures <- LUAD455_gedm_cluster3_1
source('retrievefeatures.R')
survnet455gedm <- survnet_sf[1:10]


cpyfeatures <- LUAD455ge_cluster3
source('retrievefeatures.R')
netlasso455ge <-  netlassofs
netrank455ge <- netrankfs
add455ge <- addg_sf[1:10]
coxreg455ge <- coxregfs
survnet455ge <- survnet_sf[1:10]

cpyfeatures <- LUAD455dm_cluster3
source('retrievefeatures.R')
netlasso455dm <-  netlassofs+455
netrank455dm <- netrankfs+455
add455dm <- sapply(addg_sf[1:10], function(x) x+455)
coxreg455dm <- coxregfs+455
survnet455dm <- sapply(survnet_sf[1:10], function(x) x+455)


FSFstatisticsmultivec2dt <- function(vectorlist,numfeatures) {
  temp1 <- vectorlist[[1]]
  temp2 <- vectorlist[[2]]
  temp3 <- vectorlist[[3]]
  c(sum(temp1>numfeatures)/length(temp1), sum(!is.na(match(temp1, temp2)))/length(temp1), sum(!is.na(match(temp1, temp3)))/length(temp1))
}
FSFstatisticsmultivec3dt <- function(vectorlist,numfeatures) {
  temp1 <- vectorlist[[1]]
  temp2 <- vectorlist[[2]]
  temp3 <- vectorlist[[3]]
  temp4 <- vectorlist[[4]]
  c(sum(temp1>numfeatures)/length(temp1), sum(!is.na(match(temp1, temp2)))/length(temp1), sum(!is.na(match(temp1, temp3)))/length(temp1),
    sum(!is.na(match(temp1, temp4)))/length(temp1))
}
FSFstatisticsmultinet2dt <- function(vectorlist) {
  multinet <- vectorlist[[1]]
  genenet <- vectorlist[[2]]
  methynet <- vectorlist[[3]]
  matchge <- rep(0,10)
  matchdm <- rep(0,10)
  for(i in 1:length(genenet)){
    for(j in 1:length(multinet)){
      if(sum(!is.na(match(genenet[[i]],multinet[[j]])))==length(genenet[[i]])){
        matchge[i] <- matchge[i]+1
      }
    }
  }
  for(i in 1:length(methynet)){
    for(j in 1:length(multinet)){
      if(sum(!is.na(match(methynet[[i]],multinet[[j]])))==length(methynet[[i]])){
        matchdm[i] <- matchdm[i] + 1
      }
    }
  }
  list(c(matchge,matchdm),genenet[which(matchge>0)],methynet[which(matchdm>0)])
}
FSFstatisticsmultinet3dt <- function(vectorlist) {
  multinet <- vectorlist[[1]]
  genenet <- vectorlist[[2]]
  methynet <- vectorlist[[3]]
  cnanet <- vectorlist[[4]]
  matchge <- rep(0,10)
  matchdm <- rep(0,10)
  matchcna <- rep(0,10)
  for(i in 1:length(genenet)){
    for(j in 1:length(multinet)){
      if(sum(!is.na(match(genenet[[i]],multinet[[j]])))==length(genenet[[i]])){
        matchge[i] <- matchge[i]+1
      }
    }
  }
  for(i in 1:length(methynet)){
    for(j in 1:length(multinet)){
      if(sum(!is.na(match(methynet[[i]],multinet[[j]])))==length(methynet[[i]])){
        matchdm[i] <- matchdm[i] + 1
      }
    }
  }
  for(i in 1:length(cnanet)){
    for(j in 1:length(multinet)){
      if(sum(!is.na(match(cnanet[[i]],multinet[[j]])))==length(cnanet[[i]])){
        matchcna[i] <- matchcna[i] + 1
      }
    }
  }
  list(c(matchge,matchdm,matchcna),genenet[which(matchge>0)],methynet[which(matchdm>0)],cnanet[which(matchcna>0)])
}


FSFstatisticsmultivec2dt(list(netlasso74gedm,netlasso74ge,netlasso74dm),74)
FSFstatisticsmultivec2dt(list(unique(unlist(add74gedm)),unique(unlist(add74ge)),unique(unlist(add74dm))),74)
FSFstatisticsmultivec2dt(list(netrank74gedm,netrank74ge,netrank74dm),74)
FSFstatisticsmultivec2dt(list(coxreg74gedm,coxreg74ge,coxreg74dm),74)
FSFstatisticsmultivec2dt(list(unique(unlist(survnet74gedm)),unique(unlist(survnet74ge)),unique(unlist(survnet74dm))),74)
FSFstatisticsmultinet2dt(list(add74gedm,add74ge,add74dm))
FSFstatisticsmultinet2dt(list(survnet74gedm,survnet74ge,survnet74dm))

FSFstatisticsmultivec2dt(list(netlasso123gedm,netlasso123ge,netlasso123dm),123)
FSFstatisticsmultivec2dt(list(unique(unlist(add123gedm)),unique(unlist(add123ge)),unique(unlist(add123dm))),123)
FSFstatisticsmultivec2dt(list(netrank123gedm,netrank123ge,netrank123dm),123)
FSFstatisticsmultivec2dt(list(coxreg123gedm,coxreg123ge,coxreg123dm),123)
FSFstatisticsmultivec2dt(list(unique(unlist(survnet123gedm)),unique(unlist(survnet123ge)),unique(unlist(survnet123dm))),123)
FSFstatisticsmultinet2dt(list(add123gedm,add123ge,add123dm))
FSFstatisticsmultinet2dt(list(survnet123gedm,survnet123ge,survnet123dm))

FSFstatisticsmultivec2dt(list(netlasso455gedm,netlasso455ge,netlasso455dm),455)
FSFstatisticsmultivec2dt(list(unique(unlist(add455gedm)),unique(unlist(add455ge)),unique(unlist(add455dm))),455)
FSFstatisticsmultivec2dt(list(netrank455gedm,netrank455ge,netrank455dm),455)
FSFstatisticsmultivec2dt(list(coxreg455gedm,coxreg455ge,coxreg455dm),455)
FSFstatisticsmultivec2dt(list(unique(unlist(survnet455gedm)),unique(unlist(survnet455ge)),unique(unlist(survnet455dm))),455)
FSFstatisticsmultinet2dt(list(add455gedm,add455ge,add455dm))
FSFstatisticsmultinet2dt(list(survnet455gedm,survnet455ge,survnet455dm))




candidates <- list(LUAD74_inte_level123_1,LUAD74_inte_level123_2,LUAD74_inte_level123_all,LUAD74_inte_level123_scoring1,LUAD74_inte_level123_scoring5)

# idxnetlasso <- findthebest(candidates,7)[[2]]
cpyfeatures <- LUAD74_gedmcna_cluster3_all
source('retrievefeatures.R')
netlasso74gedm <- netlassofs

# idxadd <- findthebest(candidates,11)[[2]]
cpyfeatures <- LUAD74_gedmcna_cluster3_2
source('retrievefeatures.R')
add74gedm <- addg_sf[1:10]

# idxnetrank <- findthebest(candidates,13)[[2]]
cpyfeatures <- LUAD74_gedmcna_cluster3_all
source('retrievefeatures.R')
netrank74gedm <- netrankfs

# idxcoxreg <- findthebest(candidates,17)[[2]]
cpyfeatures <- LUAD74_gedmcna_cluster3_scoring5
source('retrievefeatures.R')
coxreg74gedm <- coxregfs

# idxsurvnet <-findthebest(candidates,21)[[2]]
cpyfeatures <- LUAD74_gedmcna_cluster3_all
source('retrievefeatures.R')
survnet74gedm <- survnet_sf[1:10]


cpyfeatures <- LUAD74ge_123_cluster3
source('retrievefeatures.R')
netlasso74ge <-  netlassofs
netrank74ge <- netrankfs
add74ge <- addg_sf[1:10]
coxreg74ge <- coxregfs
survnet74ge <- survnet_sf[1:10]

cpyfeatures <- LUAD74dm_123_cluster3
source('retrievefeatures.R')
netlasso74dm <-  netlassofs+74
netrank74dm <- netrankfs+74
add74dm <- sapply(addg_sf[1:10], function(x) x+74)
coxreg74dm <- coxregfs+74
survnet74dm <- sapply(survnet_sf[1:10], function(x) x+74)

cpyfeatures <- LUAD74cna_123_cluster3
source('retrievefeatures.R')
netlasso74cna <-  netlassofs+148
netrank74cna <- netrankfs+148
add74cna <- sapply(addg_sf[1:10], function(x) x+148)
coxreg74cna <- coxregfs+148
survnet74cna <- sapply(survnet_sf[1:10], function(x) x+148)






candidates <- list(LUAD123_inte_level123_1,LUAD123_inte_level123_2,LUAD123_inte_level123_all,LUAD123_inte_level123_scoring1,LUAD123_inte_level123_scoring5)

# idxnetlasso <- findthebest(candidates,7)[[2]]
cpyfeatures <- LUAD123_gedmcna_cluster3_all
source('retrievefeatures.R')
netlasso123gedm <- netlassofs

# idxadd <- findthebest(candidates,11)[[2]]
cpyfeatures <- LUAD123_gedmcna_cluster3_2
source('retrievefeatures.R')
add123gedm <- addg_sf[1:10]

# idxnetrank <- findthebest(candidates,13)[[2]]
cpyfeatures <- LUAD123_gedmcna_cluster3_all
source('retrievefeatures.R')
netrank123gedm <- netrankfs

# idxcoxreg <- findthebest(candidates,17)[[2]]
cpyfeatures <- LUAD123_gedmcna_cluster3_2
source('retrievefeatures.R')
coxreg123gedm <- coxregfs

# idxsurvnet <-findthebest(candidates,21)[[2]]
cpyfeatures <- LUAD123_gedmcna_cluster3_2
source('retrievefeatures.R')
survnet123gedm <- survnet_sf[1:10]


cpyfeatures <- LUAD123ge_123_cluster3
source('retrievefeatures.R')
netlasso123ge <-  netlassofs
netrank123ge <- netrankfs
add123ge <- addg_sf[1:10]
coxreg123ge <- coxregfs
survnet123ge <- survnet_sf[1:10]

cpyfeatures <- LUAD123dm_123_cluster3
source('retrievefeatures.R')
netlasso123dm <-  netlassofs+123
netrank123dm <- netrankfs+123
add123dm <- sapply(addg_sf[1:10], function(x) x+123)
coxreg123dm <- coxregfs+123
survnet123dm <- sapply(survnet_sf[1:10], function(x) x+123)

cpyfeatures <- LUAD123cna_123_cluster3
source('retrievefeatures.R')
netlasso123cna <-  netlassofs+246
netrank123cna <- netrankfs+246
add123cna <- sapply(addg_sf[1:10], function(x) x+246)
coxreg123cna <- coxregfs+246
survnet123cna <- sapply(survnet_sf[1:10], function(x) x+246)




candidates <- list(LUAD455_inte_level123_1,LUAD455_inte_level123_2,LUAD455_inte_level123_all,LUAD455_inte_level123_scoring1,LUAD455_inte_level123_scoring5)

# idxnetlasso <- findthebest(candidates,7)[[2]]
cpyfeatures <- LUAD455_gedmcna_cluster3_1
source('retrievefeatures.R')
netlasso455gedm <- netlassofs

# idxadd <- findthebest(candidates,11)[[2]]
cpyfeatures <- LUAD455_gedmcna_cluster3_1
source('retrievefeatures.R')
add455gedm <- addg_sf[1:10]

# idxnetrank <- findthebest(candidates,13)[[2]]
cpyfeatures <- LUAD455_gedmcna_cluster3_scoring1
source('retrievefeatures.R')
netrank455gedm <- netrankfs

# idxcoxreg <- findthebest(candidates,17)[[2]]
cpyfeatures <- LUAD455_gedmcna_cluster3_scoring5
source('retrievefeatures.R')
coxreg455gedm <- coxregfs

# idxsurvnet <-findthebest(candidates,21)[[2]]
cpyfeatures <- LUAD455_gedmcna_cluster3_1
source('retrievefeatures.R')
survnet455gedm <- survnet_sf[1:10]


cpyfeatures <- LUAD455ge_123_cluster3
source('retrievefeatures.R')
netlasso455ge <-  netlassofs
netrank455ge <- netrankfs
add455ge <- addg_sf[1:10]
coxreg455ge <- coxregfs
survnet455ge <- survnet_sf[1:10]

cpyfeatures <- LUAD455dm_123_cluster3
source('retrievefeatures.R')
netlasso455dm <-  netlassofs+455
netrank455dm <- netrankfs+455
add455dm <- sapply(addg_sf[1:10], function(x) x+455)
coxreg455dm <- coxregfs+455
survnet455dm <- sapply(survnet_sf[1:10], function(x) x+455)

cpyfeatures <- LUAD455cna_123_cluster3
source('retrievefeatures.R')
netlasso455cna <-  netlassofs+910
netrank455cna <- netrankfs+910
add455cna <- sapply(addg_sf[1:10], function(x) x+910)
coxreg455cna <- coxregfs+910
survnet455cna <- sapply(survnet_sf[1:10], function(x) x+910)




FSFstatisticsmultivec3dt(list(netlasso74gedm,netlasso74ge,netlasso74dm,netlasso74cna),74)
FSFstatisticsmultivec3dt(list(unique(unlist(add74gedm)),unique(unlist(add74ge)),unique(unlist(add74dm)),unique(unlist(add74cna))),74)
FSFstatisticsmultivec3dt(list(netrank74gedm,netrank74ge,netrank74dm,netrank74cna),74)
FSFstatisticsmultivec3dt(list(coxreg74gedm,coxreg74ge,coxreg74dm,coxreg74cna),74)
FSFstatisticsmultivec3dt(list(unique(unlist(survnet74gedm)),unique(unlist(survnet74ge)),unique(unlist(survnet74dm)),unique(unlist(survnet74cna))),74)
FSFstatisticsmultinet3dt(list(add74gedm,add74ge,add74dm,add74cna))
FSFstatisticsmultinet3dt(list(survnet74gedm,survnet74ge,survnet74dm,survnet74cna))

FSFstatisticsmultivec3dt(list(netlasso123gedm,netlasso123ge,netlasso123dm,netlasso123cna),123)
FSFstatisticsmultivec3dt(list(unique(unlist(add123gedm)),unique(unlist(add123ge)),unique(unlist(add123dm)),unique(unlist(add123cna))),123)
FSFstatisticsmultivec3dt(list(netrank123gedm,netrank123ge,netrank123dm,netrank123cna),123)
FSFstatisticsmultivec3dt(list(coxreg123gedm,coxreg123ge,coxreg123dm,coxreg123cna),123)
FSFstatisticsmultivec3dt(list(unique(unlist(survnet123gedm)),unique(unlist(survnet123ge)),unique(unlist(survnet123dm)),unique(unlist(survnet123cna))),123)
FSFstatisticsmultinet3dt(list(add123gedm,add123ge,add123dm,add123cna))
FSFstatisticsmultinet3dt(list(survnet123gedm,survnet123ge,survnet123dm,survnet123cna))

FSFstatisticsmultivec3dt(list(netlasso455gedm,netlasso455ge,netlasso455dm,netlasso455cna),455)
FSFstatisticsmultivec3dt(list(unique(unlist(add455gedm)),unique(unlist(add455ge)),unique(unlist(add455dm)),unique(unlist(add455cna))),455)
FSFstatisticsmultivec3dt(list(netrank455gedm,netrank455ge,netrank455dm,netrank455cna),455)
FSFstatisticsmultivec3dt(list(coxreg455gedm,coxreg455ge,coxreg455dm,coxreg455cna),455)
FSFstatisticsmultivec3dt(list(unique(unlist(survnet455gedm)),unique(unlist(survnet455ge)),unique(unlist(survnet455dm)),unique(unlist(survnet455cna))),455)
FSFstatisticsmultinet3dt(list(add455gedm,add455ge,add455dm,add455cna))
FSFstatisticsmultinet3dt(list(survnet455gedm,survnet455ge,survnet455dm,survnet455cna))


#netlasso74gedm
FSFrules <- function(features){
      temp <- threslevel12[,c(1,features+1,ncol(threslevel12))]
      t <- temp
      classvar <- factor(t$days > 0,labels = c("long","short"))
      newdata <- data.frame(sapply(t[,-c(1,ncol(t))],function(x) factor(x > mean(x),labels=c("high","low"))),classvar)
      transctns <- as(newdata,"transactions")
      rules <- apriori(transctns,appearance = list(rhs = c("classvar=long", "classvar=short"),
                                                   default="lhs"),parameter = list(confidence=0.8,support = 0.1,target="rules",maxtime=10))
      
      rules <- rules[!is.redundant(rules)]
      rules <- sort(rules,by="confidence")
      list(rules,nrow(rules@quality),nrow(t))
}


threslevel1 <- read.csv('../Thesis_Results/74network/emt74_geneexp_700_1400.csv')
threslevel2 <-  read.csv('../Thesis_Results/74network/emt74_methy_700_1400.csv')
nothreslevel1 <- read.csv('../Thesis_Results/74network/emt74_geneexp_700_1400_nothres.csv')
nothreslevel2 <- read.csv('../Thesis_Results/74network/emt74_methy_700_1400_nothres.csv')

thres12samples <- intersect(threslevel1[,1],threslevel2[,1])
nothres12samples <- intersect(nothreslevel1[,1],nothreslevel2[,1])
threslevel12 <- as.data.frame(cbind(threslevel1[match(thres12samples,threslevel1[,1]),-ncol(threslevel1)], threslevel2[match(thres12samples,threslevel2[,1]),-1]))
nothreslevel12 <- as.data.frame(cbind(nothreslevel1[match(nothres12samples,nothreslevel1[,1]),-ncol(nothreslevel1)], 
                                      nothreslevel2[match(nothres12samples,nothreslevel2[,1]),-1]))
methyfeaturenames <-sapply(names(threslevel1)[-c(1,ncol(threslevel1))], function(x) paste("mty_", x, sep=""))
featurenames <- c(names(threslevel1)[-c(ncol(threslevel1))],methyfeaturenames,"days")
names(threslevel12) <- featurenames
names(nothreslevel12) <- featurenames

threslevel12[,-c(1,ncol(threslevel12))] <- scale(threslevel12[,-c(1,ncol(threslevel12))])
nothreslevel12[,-c(1,ncol(nothreslevel12))] <-  scale(nothreslevel12[,-c(1,ncol(nothreslevel12))])

featurelist <- list(netlasso74ge,netlasso74dm,netlasso74gedm)
netlassorules74 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(unique(unlist(add74ge)),unique(unlist(add74dm)),unique(unlist(add74gedm)))
netaddrules74 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(netrank74ge,netrank74dm,netrank74gedm)
netrankrules74 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(coxreg74ge,coxreg74dm,coxreg74gedm)
coxregrules74 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(unique(unlist(survnet74ge)),unique(unlist(survnet74dm)),unique(unlist(survnet74gedm)))
netsurvnetrules74 <- lapply(featurelist, function(x) FSFrules(x))



threslevel1 <- read.csv('../Thesis_Results/123network/emt123_geneexp_700_1400.csv')
threslevel2 <-  read.csv('../Thesis_Results/123network/emt123_methy_700_1400.csv')
nothreslevel1 <- read.csv('../Thesis_Results/123network/emt123_geneexp_700_1400_nothres.csv')
nothreslevel2 <- read.csv('../Thesis_Results/123network/emt123_methy_700_1400_nothres.csv')

thres12samples <- intersect(threslevel1[,1],threslevel2[,1])
nothres12samples <- intersect(nothreslevel1[,1],nothreslevel2[,1])
threslevel12 <- as.data.frame(cbind(threslevel1[match(thres12samples,threslevel1[,1]),-ncol(threslevel1)], threslevel2[match(thres12samples,threslevel2[,1]),-1]))
nothreslevel12 <- as.data.frame(cbind(nothreslevel1[match(nothres12samples,nothreslevel1[,1]),-ncol(nothreslevel1)], 
                                      nothreslevel2[match(nothres12samples,nothreslevel2[,1]),-1]))
methyfeaturenames <-sapply(names(threslevel1)[-c(1,ncol(threslevel1))], function(x) paste("mty_", x, sep=""))
featurenames <- c(names(threslevel1)[-c(ncol(threslevel1))],methyfeaturenames,"days")
names(threslevel12) <- featurenames
names(nothreslevel12) <- featurenames


threslevel12[,-c(1,ncol(threslevel12))] <- scale(threslevel12[,-c(1,ncol(threslevel12))])
nothreslevel12[,-c(1,ncol(nothreslevel12))] <-  scale(nothreslevel12[,-c(1,ncol(nothreslevel12))])

featurelist <- list(netlasso123ge,netlasso123dm,netlasso123gedm)
netlassorules123 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(unique(unlist(add123ge)),unique(unlist(add123dm)),unique(unlist(add123gedm)))
netaddrules123 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(netrank123ge,netrank123dm,netrank123gedm)
netrankrules123 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(coxreg123ge,coxreg123dm,coxreg123gedm)
coxregrules123 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(unique(unlist(survnet123ge)),unique(unlist(survnet123dm)),unique(unlist(survnet123gedm)))
netsurvnetrules123 <- lapply(featurelist, function(x) FSFrules(x))


threslevel1 <- read.csv('../Thesis_Results/455network/emt455_geneexp_700_1400.csv')
threslevel2 <-  read.csv('../Thesis_Results/455network/emt455_methy_700_1400.csv')
# threslevel3 <- read.csv('../Thesis_Results/455network/emt445_cna_700_1400.csv')
nothreslevel1 <- read.csv('../Thesis_Results/455network/emt455_geneexp_700_1400_nothres.csv')
nothreslevel2 <- read.csv('../Thesis_Results/455network/emt455_methy_700_1400_nothres.csv')
# nothreslevel3 <- read.csv('../Thesis_Results/455network/emt445_cna_700_1400_nothres.csv')

thres12samples <- intersect(threslevel1[,1],threslevel2[,1])
nothres12samples <- intersect(nothreslevel1[,1],nothreslevel2[,1])
threslevel12 <- as.data.frame(cbind(threslevel1[match(thres12samples,threslevel1[,1]),-ncol(threslevel1)], threslevel2[match(thres12samples,threslevel2[,1]),-1]))
nothreslevel12 <- as.data.frame(cbind(nothreslevel1[match(nothres12samples,nothreslevel1[,1]),-ncol(nothreslevel1)], 
                                      nothreslevel2[match(nothres12samples,nothreslevel2[,1]),-1]))
methyfeaturenames <-sapply(names(threslevel2)[-c(1,ncol(threslevel2))], function(x) paste("mty_", x, sep=""))
featurenames <- c(names(threslevel1)[-c(ncol(threslevel1))],methyfeaturenames,"days")
names(threslevel12) <- featurenames
names(nothreslevel12) <- featurenames

threslevel12[,-c(1,ncol(threslevel12))] <- scale(threslevel12[,-c(1,ncol(threslevel12))])
nothreslevel12[,-c(1,ncol(nothreslevel12))] <-  scale(nothreslevel12[,-c(1,ncol(nothreslevel12))])

featurelist <- list(netlasso455ge,netlasso455dm,netlasso455gedm)
netlassorules455 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(unique(unlist(add455ge)),unique(unlist(add455dm)),unique(unlist(add455gedm)))
netaddrules455 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(netrank455ge,netrank455dm,netrank455gedm)
netrankrules455 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(coxreg455ge,coxreg455dm,coxreg455gedm)
coxregrules455 <- lapply(featurelist, function(x) FSFrules(x))

featurelist <- list(unique(unlist(survnet455ge)),unique(unlist(survnet455dm)),unique(unlist(survnet455gedm)))
netsurvnetrules455 <- lapply(featurelist, function(x) FSFrules(x))

rulelists <- list(netlassorules74,netaddrules74,netrankrules74,coxregrules74,netsurvnetrules74,
             netlassorules123,netaddrules123,netrankrules123,coxregrules123,netsurvnetrules123,
             netlassorules455,netaddrules455,netrankrules455,coxregrules455,netsurvnetrules455)

ratiosmatrix2dt <- matrix(0,15,3)
supportmatrix2dt <- matrix(0,15,3)
confidencematrix2dt <- matrix(0,15,3)
liftmatrix2dt <- matrix(0,15,3)
rulelengthmatrix2dt <- matrix(0,15,3)

rulestatistics2 <- function(ruleset,ntop){
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
  
  if(totalnumofrules < ntop)
    top <- totalnumofrules else
      top <- ntop
  avgvalues <- colMeans(ruleset@quality[1:top,])
  avgrulelength<- mean(colSums((ruleset@lhs@data)))
  
  list(numofruleslong,numofrulesshort,rulelengthlong,rulelengthshort,avgvalues,avgrulelength)
}

for(j in 1:15){
  for (i in 1:3){
    temp <- rulestatistics2(rulelists[[j]][[i]][[1]],20)
    ratiosmatrix2dt[j,i] <- temp[[1]]/temp[[2]]
    supportmatrix2dt[j,i] <- unname(temp[[5]][1])
    confidencematrix2dt[j,i] <- unname(temp[[5]][2])
    liftmatrix2dt[j,i] <- unname(temp[[5]][3])
    rulelengthmatrix2dt[j,i] <- temp[[6]]
  }
}

meanval <- as.vector(confidencematrix2dt)
Algorithms <- factor(rep(c(1,2,3,4,5),9),labels = c("NetLasso","addDA2","Netrank","RegCox","Survnet"))
Networks <- factor(rep(c(rep(1,5),rep(2,5),rep(3,5)),3),labels = c("Core Network","Filtered Network","Extended Network"))
`Data Levels` <-  factor(c(rep(1,15),rep(2,15),rep(3,15)), 
                         labels = c("GE", "DM", "Multiplex"))
meanaccu <- data.frame(meanval,Algorithms,Networks,`Data Levels`)
names(meanaccu)[4] <- "Data Levels"
ggplot(meanaccu, aes(fill=`Data Levels`,y=meanval,x=Algorithms)) +
  geom_bar(position = "dodge",stat="identity") +
  scale_y_continuous(name = "Average Confidence of Rules") +
  scale_x_discrete(name = "Feature Selection Algorithms") +
  theme_bw()+
  theme(plot.title = element_text(size = 15, family = "Times", face = "bold"),
        text = element_text(size = 14, family = "Times"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12),
        legend.position = "top") + 
  coord_cartesian(ylim=c(0.83,1)) + 
  facet_grid(Networks ~ .) + 
  scale_fill_brewer(palette="Set3") 



