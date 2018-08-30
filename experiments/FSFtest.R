
LEnet <- LEnet455
numfeatures <- 455
thresdata <- read.csv('../Thesis_Results/455network/emt455_geneexp_700_1400.csv')
nothresdata <- read.csv('../Thesis_Results/455network/emt455_geneexp_700_1400_nothres.csv')
LUAD455_geneexp_700_1400_6algs <- foreach(x = 1:50,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_geneexp_flds_700_1400_thres[[x]]
  nothrestrainidx <- LUAD_geneexp_flds_700_1400_nothres[[x]]
  evaluation_algorithms_6algs(thresdata,threstrainidx,nothresdata,nothrestrainidx)
}
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
LUAD455_methy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_flds_700_1400_thres[[x]]
  nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
  evaluation_algorithms(thresdata,threstrainidx,nothresdata,nothrestrainidx)
}
thresdata2 <- getclinicaldata(1400,700,thresdata)
cpy <- LUAD455_methy_700_1400
source('selectedfeatures_clustering.R')
LUAD455_FSFmethy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_flds_700_1400_thres[[x]]
  nothrestrainidx <- LUAD_methy_flds_700_1400_nothres[[x]]
  evaluation_algorithms_FSF(thresdata,threstrainidx)
}
LUAD455_methy_cluster2 <- cluster_measure(nothresdata,2,add_dire455[[2]])
LUAD455_methy_cluster3 <- cluster_measure(nothresdata,3,add_dire455[[2]])
LUAD455_climethy_700_1400 <- foreach(x = 1:repetition,.errorhandling='remove')  %dopar% {
  threstrainidx <-LUAD_methy_cliflds700_1400[[x]]
  single_level_methods1_cliplus_2sets(thresdata2[[1]],thresdata2[[2]],threstrainidx,thresdata2[[3]],20,10)
}




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

