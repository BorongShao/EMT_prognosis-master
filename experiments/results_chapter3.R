#--------------------74 nodes network-------------------------------#

cpy <- LUAD74_inte_level13_1
lengthcpy <- length(cpy)
cor1net74 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cor1net74) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74_inte_level13_2
lengthcpy <- length(cpy)
cor2net74 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cor2net74) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74_inte_level13_all
lengthcpy <- length(cpy)
corallnet74 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(corallnet74) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74_inte_level13_scoring1
lengthcpy <- length(cpy)
score1net74 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(score1net74) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74_inte_level13_scoring5
lengthcpy <- length(cpy)
score5net74 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(score5net74) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74ge
lengthcpy <- length(cpy)
ge74 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(ge74) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74cna
lengthcpy <- length(cpy)
cna74 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cna74) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

# 74 nodes network _ clinical features

cpy <- LUAD74_intecli_level13_1
lengthcpy <- length(cpy)
clicor1net74 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicor1net74) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74_intecli_level13_2
lengthcpy <- length(cpy)
clicor2net74 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicor2net74) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74_intecli_level13_all
lengthcpy <- length(cpy)
clicorallnet74 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicorallnet74) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74_intecli_level13_scoring1
lengthcpy <- length(cpy)
cliscore1net74 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cliscore1net74) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74_intecli_level13_scoring5
lengthcpy <- length(cpy)
cliscore5net74 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cliscore5net74) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74ge_cli
lengthcpy <- length(cpy)
clige74 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clige74) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD74cna_cli
lengthcpy <- length(cpy)
clicna74 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicna74) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

p74cor1cluster2 <- sapply(1:12, function(x) LUAD74_gecna_cluster2_1[[1]][[x]][[1]])
names(p74cor1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74cor1cluster3 <- sapply(1:12, function(x) LUAD74_gecna_cluster3_1[[1]][[x]][[1]])
names(p74cor1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")
temp

p74cor2cluster2 <- sapply(1:12, function(x) LUAD74_gecna_cluster2_2[[1]][[x]][[1]])
names(p74cor2cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74cor2cluster3 <- sapply(1:12, function(x) LUAD74_gecna_cluster3_2[[1]][[x]][[1]])
names(p74cor2cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74allcluster2 <- sapply(1:12, function(x) LUAD74_gecna_cluster2_all[[1]][[x]][[1]])
names(p74allcluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74allcluster3 <- sapply(1:12, function(x) LUAD74_gecna_cluster3_all[[1]][[x]][[1]])
names(p74allcluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74score1cluster2 <- sapply(1:12, function(x) LUAD74_gecna_cluster2_scoring1[[1]][[x]][[1]])
names(p74score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74score1cluster3 <- sapply(1:12, function(x) LUAD74_gecna_cluster3_scoring1[[1]][[x]][[1]])
names(p74score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74score5cluster2 <- sapply(1:12, function(x) LUAD74_gecna_cluster2_scoring5[[1]][[x]][[1]])
names(p74score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74score5cluster3 <- sapply(1:12, function(x) LUAD74_gecna_cluster3_scoring5[[1]][[x]][[1]])
names(p74score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74gecluster2 <- sapply(1:12, function(x) LUAD74ge_cluster2[[1]][[x]][[1]])
names(p74score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74gecluster3 <- sapply(1:12, function(x) LUAD74ge_cluster3[[1]][[x]][[1]])
names(p74score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74cnacluster2 <- sapply(1:12, function(x) LUAD74cna_cluster2[[1]][[x]][[1]])
names(p74score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p74cnacluster3 <- sapply(1:12, function(x) LUAD74cna_cluster3[[1]][[x]][[1]])
names(p74score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

auc74 <- rbind(ge74,cna74,cor2net74,cor1net74,score1net74,score5net74,corallnet74)
cliauc74 <- rbind(clige74,clicna74,clicor2net74,clicor1net74,cliscore1net74,cliscore5net74,clicorallnet74)
pvalue74cluster2 <- rbind(p74gecluster2,p74cnacluster2,p74cor2cluster2,p74cor1cluster2,p74score1cluster2,p74score5cluster2,p74allcluster2)
pvalue74cluster3 <- rbind(p74gecluster3,p74cnacluster3,p74cor2cluster3,p74cor1cluster3,p74score1cluster3,p74score5cluster3,p74allcluster3)


#--------------------123 nodes network-------------------------------#

cpy <- LUAD123_inte_level13_1
lengthcpy <- length(cpy)
cor1net123 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cor1net123) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123_inte_level13_2
lengthcpy <- length(cpy)
cor2net123 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cor2net123) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123_inte_level13_all
lengthcpy <- length(cpy)
corallnet123 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(corallnet123) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123_inte_level13_scoring1
lengthcpy <- length(cpy)
score1net123 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(score1net123) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123_inte_level13_scoring5
lengthcpy <- length(cpy)
score5net123 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(score5net123) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123ge
lengthcpy <- length(cpy)
ge123 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(ge123) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123cna
lengthcpy <- length(cpy)
cna123 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cna123) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

# 123 nodes network _ clinical features

cpy <- LUAD123_intecli_level13_1
lengthcpy <- length(cpy)
clicor1net123 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicor1net123) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123_intecli_level13_2
lengthcpy <- length(cpy)
clicor2net123 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicor2net123) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123_intecli_level13_all
lengthcpy <- length(cpy)
clicorallnet123 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicorallnet123) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123_intecli_level13_scoring1
lengthcpy <- length(cpy)
cliscore1net123 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cliscore1net123) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123_intecli_level13_scoring5
lengthcpy <- length(cpy)
cliscore5net123 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cliscore5net123) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123ge_cli
lengthcpy <- length(cpy)
clige123 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clige123) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD123cna_cli
lengthcpy <- length(cpy)
clicna123 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicna123) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")


p123cor1cluster2 <- sapply(1:12, function(x) LUAD123_gecna_cluster2_1[[1]][[x]][[1]])
names(p123cor1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123cor1cluster3 <- sapply(1:12, function(x) LUAD123_gecna_cluster3_1[[1]][[x]][[1]])
names(p123cor1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")
temp

p123cor2cluster2 <- sapply(1:12, function(x) LUAD123_gecna_cluster2_2[[1]][[x]][[1]])
names(p123cor2cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123cor2cluster3 <- sapply(1:12, function(x) LUAD123_gecna_cluster3_2[[1]][[x]][[1]])
names(p123cor2cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123allcluster2 <- sapply(1:12, function(x) LUAD123_gecna_cluster2_all[[1]][[x]][[1]])
names(p123allcluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123allcluster3 <- sapply(1:12, function(x) LUAD123_gecna_cluster3_all[[1]][[x]][[1]])
names(p123allcluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123score1cluster2 <- sapply(1:12, function(x) LUAD123_gecna_cluster2_scoring1[[1]][[x]][[1]])
names(p123score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123score1cluster3 <- sapply(1:12, function(x) LUAD123_gecna_cluster3_scoring1[[1]][[x]][[1]])
names(p123score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123score5cluster2 <- sapply(1:12, function(x) LUAD123_gecna_cluster2_scoring5[[1]][[x]][[1]])
names(p123score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123score5cluster3 <- sapply(1:12, function(x) LUAD123_gecna_cluster3_scoring5[[1]][[x]][[1]])
names(p123score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123gecluster2 <- sapply(1:12, function(x) LUAD123ge_cluster2[[1]][[x]][[1]])
names(p123score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123gecluster3 <- sapply(1:12, function(x) LUAD123ge_cluster3[[1]][[x]][[1]])
names(p123score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123cnacluster2 <- sapply(1:12, function(x) LUAD123cna_cluster2[[1]][[x]][[1]])
names(p123score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p123cnacluster3 <- sapply(1:12, function(x) LUAD123cna_cluster3[[1]][[x]][[1]])
names(p123score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

auc123 <- rbind(ge123,cna123,cor2net123,cor1net123,score1net123,score5net123,corallnet123)
cliauc123 <- rbind(clige123,clicna123,clicor2net123,clicor1net123,cliscore1net123,cliscore5net123,clicorallnet123)
pvalue123cluster2 <- rbind(p123gecluster2,p123cnacluster2,p123cor2cluster2,p123cor1cluster2,p123score1cluster2,p123score5cluster2,p123allcluster2)
pvalue123cluster3 <- rbind(p123gecluster3,p123cnacluster3,p123cor2cluster3,p123cor1cluster3,p123score1cluster3,p123score5cluster3,p123allcluster3)











#
#--------------------455 nodes network-------------------------------#

cpy <- LUAD455_inte_level13_1
lengthcpy <- length(cpy)
cor1net455 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cor1net455) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455_inte_level13_2
lengthcpy <- length(cpy)
cor2net455 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cor2net455) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455_inte_level13_all
lengthcpy <- length(cpy)
corallnet455 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(corallnet455) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455_inte_level13_scoring1
lengthcpy <- length(cpy)
score1net455 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(score1net455) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455_inte_level13_scoring5
lengthcpy <- length(cpy)
score5net455 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(score5net455) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455ge
lengthcpy <- length(cpy)
ge455 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(ge455) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455cna
lengthcpy <- length(cpy)
cna455 <- sapply(seq(1,21,2), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cna455) <- c("EMT","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

# 455 nodes network _ clinical features

cpy <- LUAD455_intecli_level13_1
lengthcpy <- length(cpy)
clicor1net455 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicor1net455) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455_intecli_level13_2
lengthcpy <- length(cpy)
clicor2net455 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicor2net455) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455_intecli_level13_all
lengthcpy <- length(cpy)
clicorallnet455 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicorallnet455) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455_intecli_level13_scoring1
lengthcpy <- length(cpy)
cliscore1net455 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cliscore1net455) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455_intecli_level13_scoring5
lengthcpy <- length(cpy)
cliscore5net455 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(cliscore5net455) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455ge_cli
lengthcpy <- length(cpy)
clige455 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clige455) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

cpy <- LUAD455cna_cli
lengthcpy <- length(cpy)
clicna455 <- sapply(c(24,2:11), function(y) mean(sapply(1:lengthcpy, function(x) cpy[[x]][[y]][[1]])))           
names(clicna455) <- c("clinical","t-test","Lasso", "NetLasso", "stSVM","addDA2","netrank","Cox","RegCox","MSS","Survnet")

p455cor1cluster2 <- sapply(1:12, function(x) LUAD455_gecna_cluster2_1[[1]][[x]][[1]])
names(p455cor1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455cor1cluster3 <- sapply(1:12, function(x) LUAD455_gecna_cluster3_1[[1]][[x]][[1]])
names(p455cor1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")
temp

p455cor2cluster2 <- sapply(1:12, function(x) LUAD455_gecna_cluster2_2[[1]][[x]][[1]])
names(p455cor2cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455cor2cluster3 <- sapply(1:12, function(x) LUAD455_gecna_cluster3_2[[1]][[x]][[1]])
names(p455cor2cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455allcluster2 <- sapply(1:12, function(x) LUAD455_gecna_cluster2_all[[1]][[x]][[1]])
names(p455allcluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455allcluster3 <- sapply(1:12, function(x) LUAD455_gecna_cluster3_all[[1]][[x]][[1]])
names(p455allcluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455score1cluster2 <- sapply(1:12, function(x) LUAD455_gecna_cluster2_scoring1[[1]][[x]][[1]])
names(p455score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455score1cluster3 <- sapply(1:12, function(x) LUAD455_gecna_cluster3_scoring1[[1]][[x]][[1]])
names(p455score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455score5cluster2 <- sapply(1:12, function(x) LUAD455_gecna_cluster2_scoring5[[1]][[x]][[1]])
names(p455score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455score5cluster3 <- sapply(1:12, function(x) LUAD455_gecna_cluster3_scoring5[[1]][[x]][[1]])
names(p455score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455gecluster2 <- sapply(1:12, function(x) LUAD455ge_cluster2[[1]][[x]][[1]])
names(p455score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455gecluster3 <- sapply(1:12, function(x) LUAD455ge_cluster3[[1]][[x]][[1]])
names(p455score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455cnacluster2 <- sapply(1:12, function(x) LUAD455cna_cluster2[[1]][[x]][[1]])
names(p455score1cluster2) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

p455cnacluster3 <- sapply(1:12, function(x) LUAD455cna_cluster3[[1]][[x]][[1]])
names(p455score1cluster3) <- c("EMT","t-test","Lasso", "NetLasso","addDA2","netrank","stSVM","Cox","RegCox","MSS","Survnet","ensemble")

auc455 <- rbind(ge455,cna455,cor2net455,cor1net455,score1net455,score5net455,corallnet455)
cliauc455 <- rbind(clige455,clicna455,clicor2net455,clicor1net455,cliscore1net455,cliscore5net455,clicorallnet455)
# pvalue455cluster2 <- rbind(p455gecluster2,p455cnacluster2,p455cor2cluster2,p455cor1cluster2,p455score1cluster2,p455score5cluster2,p455allcluster2)
# pvalue455cluster3 <- rbind(p455gecluster3,p455cnacluster3,p455cor2cluster3,p455cor1cluster3,p455score1cluster3,p455score5cluster3,p455allcluster3)
pvalue455cluster2 <- rbind(p455gecluster2,p455cnacluster2,p455cor2cluster2,p455cor1cluster2,p455score1cluster2,p455score5cluster2)
pvalue455cluster3 <- rbind(p455gecluster3,p455cnacluster3,p455cor2cluster3,p455cor1cluster3,p455score1cluster3,p455score5cluster3)

