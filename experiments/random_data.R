thresupper = 1200
threslower = 900
alivesamples <- read.csv(file='../Experiments/LUAD/Data/clinical/alivesamples.txt',header=TRUE, sep=' ')
deadsamples <- read.csv(file='../Experiments/LUAD/Data/clinical/deadsamples.txt',header=TRUE, sep=' ')
allsamples <- rbind(deadsamples,alivesamples)
alivesensor <- alivesamples[alivesamples$days >= thresupper,]
#alivesensor <- alivesamples
samples <- rbind(deadsamples,alivesensor)
samplesthres <- samples[samples$days >= thresupper | samples$days <threslower, ]

# all RNA-Seq features
toremove <- which(sapply(LErnaseq_all[,-1], var) < 0.5)+1
allrnaseqfeatures <- data.frame(LErnaseq_all[match(rmoutlier[,1], LErnaseq_all[,1]),-toremove], rmoutlier$days)
names(allrnaseqfeatures)[ncol(allrnaseqfeatures)] <- "days"
allrnaseqfeatures[,-c(1,ncol(allrnaseqfeatures))] <- log2(allrnaseqfeatures[,-c(1,ncol(allrnaseqfeatures))]+1)

random430auc <- rep(1000,0)
random430aupr <- rep(1000,0)

for(i in 1:1000){
  # take the same samples
  randidx <- sample(ncol(allrnaseqfeatures)-2,414)+1
  LErnaseq_rand <- allrnaseqfeatures[match(rmoutlier[,1], allrnaseqfeatures[,1]),c(1,randidx)]
  # names(LErnaseq_rand)[ncol(LErnaseq_rand)] <- "days"
  randidx <- sample(ncol(LEmirna_all)-1,16)+1
  LEmirna_rand <- LEmirna_all[match(rmoutlier[,1], LEmirna_all[,1]),c(randidx)]
  rnaseq_mirna_rand <- data.frame(LErnaseq_rand,LEmirna_rand, rmoutlier$days)
  rnaseq_mirna_rand[,-c(1,ncol(rnaseq_mirna_rand))] <- log2(rnaseq_mirna_rand[,-c(1,ncol(rnaseq_mirna_rand))]+1)
  LEthresgeneexp_rand <- rnaseq_mirna_rand
  names(LEthresgeneexp_rand)[length(LEthresgeneexp_rand)] <- "days"
  
  
  MCCV_random <- foreach(x = 1:100,.errorhandling='remove')  %dopar% {
    threstrainidx <-sampling(LEthresgeneexp_rand,0.9)
    trainset <- LEthresgeneexp_rand[threstrainidx,]
    testset <- LEthresgeneexp_rand[-threstrainidx,]
    svm1(trainset[,-c(1,ncol(trainset))],trainset$days, testset[,-c(1,ncol(testset))], testset$days)
  }
  
  random430auc[i] <- mean(sapply(MCCV_random, function(x) x[1]))
  random430aupr[i] <- mean(sapply(MCCV_random, function(x) x[2]))
}

#kernel density plot
d <- density(random430auc)
plot(d, main="Kernel Density plot of the AUC values of 100 groups of random features", xlab="AUC values")
polygon(d,col="red",border="blue")




