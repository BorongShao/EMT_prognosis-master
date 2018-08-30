#degree centrality
emtnetwork <- LEnet455

for(j in 1:10){
  if(j==4 | j==10){
    level1nodes <-  unique(unlist(LUAD74_geneexp_cluster3[[2]][[j]]))
    level2nodes <- unique(unlist(LUAD74_methy_cluster3[[2]][[j]]))
    level3nodes <- match(V(LEnet74_cna)$name[unique(unlist(LUAD74_cna_cluster3[[2]][[j]]))], V(emtnetwork)$name) 
  } else{
    level1nodes <-  LUAD74_geneexp_cluster3[[2]][[j]]
    level2nodes <- LUAD74_methy_cluster3[[2]][[j]]
    level3nodes <-  match(V(LEnet74_cna)$name[LUAD74_cna_cluster3[[2]][[j]]], V(emtnetwork)$name) 
  }

}

j=4
level1nodes <-  LUAD455_geneexp_cluster3[[2]][[j]]
level2nodes <- LUAD455_methy_cluster3[[2]][[j]]
level3nodes <-  match(V(LEnet455_cna)$name[LUAD455_cna_cluster3[[2]][[j]]], V(emtnetwork)$name)

level1nodes <-  unique(unlist(LUAD455_geneexp_cluster3[[2]][[j]]))
level2nodes <- unique(unlist(LUAD455_methy_cluster3[[2]][[j]]))
level3nodes <- match(V(LEnet455_cna)$name[unique(unlist(LUAD455_cna_cluster3[[2]][[j]]))], V(emtnetwork)$name) 

#5 tables for 5 measures, each with 3 levels and 11 algorithms. 
#centrality
centrality_measure <- function(FUN){
  table_centrality <- matrix(0,11,9)
  for(j in 1:3){
    if(j==1){
      emtnetwork <- LEnet74
      temp <- LUAD74_geneexp_cluster3
    }else if(j==2){
      emtnetwork <- LEnet123
      temp <- LUAD123_geneexp_cluster3
    }else if(j==3){
      emtnetwork <- LEnet455
      temp <- LUAD455_geneexp_cluster3
    }
    for(i in 1:11){
      if(i !=4 & i!=10){
        table_centrality[i,j] <- mean(FUN(emtnetwork)$vector[temp[[2]][[i]]])
      }else{
        table_centrality[i,j] <- mean(FUN(emtnetwork)$vector[unique(unlist(temp[[2]][[i]]))])
      }
    }
  }
  
  for(j in 4:6){
    if(j==4){
      emtnetwork <- LEnet74
      temp <- LUAD74_methy_cluster3
    }else if(j==5){
      emtnetwork <- LEnet123
      temp <- LUAD123_methy_cluster3
    }else if(j==6){
      emtnetwork <- LEnet455
      temp <- LUAD455_methy_cluster3
    }
    for(i in 1:11){
      if(i !=4 & i!=10){
        table_centrality[i,j] <- mean(FUN(emtnetwork)$vector[temp[[2]][[i]]])
      }else{
        table_centrality[i,j] <- mean(FUN(emtnetwork)$vector[unique(unlist(temp[[2]][[i]]))])
      }
    }
  }
  
  for(j in 7:9){
    if(j==7){
      emtnetwork <- LEnet74
      temp <- LUAD74_cna_cluster3
    }else if(j==8){
      emtnetwork <- LEnet123
      temp <- LUAD123_cna_cluster3
    }else if(j==9){
      emtnetwork <- LEnet455
      temp <- LUAD455_cna_cluster3
    }
    for(i in 1:11){
      if(i !=4 & i!=10){
        table_centrality[i,j] <- mean(FUN(emtnetwork)$vector[temp[[2]][[i]]])
      }else{
        table_centrality[i,j] <- mean(FUN(emtnetwork)$vector[unique(unlist(temp[[2]][[i]]))])
      }
    }
  }
  table_centrality
}



# centr_degree(emtnetwork)[[1]][level1nodes]
# betweenness(emtnetwork)[level1nodes]
# closeness(emtnetwork)[level1nodes]
# eccentricity(emtnetwork)[level1nodes]
# eigen_centrality(emtnetwork)$vector[level1nodes]

table_degree <- centrality_measure(centr_degree)
table_betweenness <- centrality_measure(betweenness)
table_closeness <- centrality_measure(closeness)
table_eccentricity <- centrality_measure(eccentricity)
table_eigen <- centrality_measure(eigen_centrality)

table_degree <- as.data.frame(table_degree)
names(table_degree) <- c("GE74","GE123","GE455","DM74","DM123","DM455","CA70","CA117","CA445")

table_betweenness <- as.data.frame(table_betweenness)
names(table_betweenness) <- c("GE74","GE123","GE455","DM74","DM123","DM455","CA70","CA117","CA445")

# table_closeness <- as.data.frame(table_closeness)
# names(table_closeness) <- c("GE74","GE123","GE455","DM74","DM123","DM455","CA70","CA117","CA445")

# table_eccentricity <- as.data.frame(table_eccentricity)
# names(table_eccentricity) <- c("GE74","GE123","GE455","DM74","DM123","DM455","CA70","CA117","CA445")

# table_eigen <- as.data.frame(table_eigen)
# names(table_eigen) <- c("GE74","GE123","GE455","DM74","DM123","DM455","CA70","CA117","CA445")


# > write.csv(table_degree[,c(1,4,7,2,5,8,3,6,9)], file='../../Desktop/degree.csv', row.names = algorithmstrings)
# > write.csv(table_betweenness[,c(1,4,7,2,5,8,3,6,9)], file='../../Desktop/betweenness.csv', row.names = algorithmstrings)


#use Kruskal-Wallis test and Mann-Whitney U test to determine the difference
# algorithms overall or for each algorithm?

#Kruskal-Wallis test
kruskal.test(c(table_betweenness[,1], table_betweenness[,4],table_betweenness[,7]), c(rep(1,11),rep(2,11),rep(3,11))) 
#Kruskal-Wallis chi-squared = 2.8721, df = 2, p-value = 0.2379
kruskal.test(c(table_betweenness[,2], table_betweenness[,5],table_betweenness[,8]), c(rep(1,11),rep(2,11),rep(3,11))) 
#Kruskal-Wallis chi-squared = 8.8751, df = 2, p-value = 0.01183
kruskal.test(c(table_betweenness[,3], table_betweenness[,6],table_betweenness[,9]), c(rep(1,11),rep(2,11),rep(3,11))) 
#Kruskal-Wallis chi-squared = 7.2591, df = 2, p-value = 0.02653

#Mann-Whitney U test

wilcox.test(table_betweenness[,1],table_betweenness[,4],paired = TRUE) #V = 37, p-value = 0.7646
wilcox.test(table_betweenness[,1],table_betweenness[,7],paired = TRUE) # V = 56, p-value = 0.04199
wilcox.test(table_betweenness[,4],table_betweenness[,7],paired = TRUE) #V = 46, p-value = 0.2783


wilcox.test(table_betweenness[,2],table_betweenness[,5],paired = TRUE) #V = 1, p-value = 0.001953
wilcox.test(table_betweenness[,2],table_betweenness[,8],paired = TRUE) # V = 10, p-value = 0.04199
wilcox.test(table_betweenness[,5],table_betweenness[,8],paired = TRUE) # V = 66, p-value = 0.0009766

wilcox.test(table_betweenness[,3],table_betweenness[,6],paired = TRUE) # V = 6, p-value = 0.01367
wilcox.test(table_betweenness[,3],table_betweenness[,9],paired = TRUE) # V = 24, p-value = 0.4648
wilcox.test(table_betweenness[,6],table_betweenness[,9],paired = TRUE) # V = 56, p-value = 0.04199

#Kruskal-Wallis test
kruskal.test(c(table_degree[,1], table_degree[,4],table_degree[,7]), c(rep(1,11),rep(2,11),rep(3,11))) # 0.26
kruskal.test(c(table_degree[,2], table_degree[,5],table_degree[,8]), c(rep(1,11),rep(2,11),rep(3,11))) 
#Kruskal-Wallis chi-squared = 7.6637, df = 2, p-value = 0.02167
kruskal.test(c(table_degree[,3], table_degree[,6],table_degree[,9]), c(rep(1,11),rep(2,11),rep(3,11))) 
# Kruskal-Wallis chi-squared = 6.6038, df = 2, p-value = 0.03681

#Mann-Whitney U test

wilcox.test(table_degree[,1],table_degree[,4],paired = TRUE) # V = 23, p-value = 0.5286
wilcox.test(table_degree[,1],table_degree[,7],paired = TRUE) # V = 55, p-value = 0.05545
wilcox.test(table_degree[,4],table_degree[,7],paired = TRUE) # V = 47.5, p-value = 0.213

wilcox.test(table_degree[,2],table_degree[,5],paired = TRUE) # V = 2, p-value = 0.00293
wilcox.test(table_degree[,2],table_degree[,8],paired = TRUE) # V = 3, p-value = 0.008686
wilcox.test(table_degree[,5],table_degree[,8],paired = TRUE) # V = 61, p-value = 0.009766

wilcox.test(table_degree[,3],table_degree[,6],paired = TRUE) # V = 5, p-value = 0.009766
wilcox.test(table_degree[,3],table_degree[,9],paired = TRUE) # V = 26, p-value = 0.5771
wilcox.test(table_degree[,6],table_degree[,9],paired = TRUE) # V = 59, p-value = 0.01855

# closeness (modest)
wilcox.test(table_closeness[,1],table_closeness[,4],paired = TRUE) # 
wilcox.test(table_closeness[,1],table_closeness[,7],paired = TRUE) # 
wilcox.test(table_closeness[,4],table_closeness[,7],paired = TRUE) # 

wilcox.test(table_closeness[,2],table_closeness[,5],paired = TRUE) # 
wilcox.test(table_closeness[,2],table_closeness[,8],paired = TRUE) #
wilcox.test(table_closeness[,5],table_closeness[,8],paired = TRUE) # 

wilcox.test(table_closeness[,3],table_closeness[,6],paired = TRUE) # 
wilcox.test(table_closeness[,3],table_closeness[,9],paired = TRUE) # 
wilcox.test(table_closeness[,6],table_closeness[,9],paired = TRUE) # 

# eccentricity (modest)
wilcox.test(table_eccentricity[,1],table_eccentricity[,4],paired = TRUE) # 
wilcox.test(table_eccentricity[,1],table_eccentricity[,7],paired = TRUE) #
wilcox.test(table_eccentricity[,4],table_eccentricity[,7],paired = TRUE) # 

wilcox.test(table_eccentricity[,2],table_eccentricity[,5],paired = TRUE) # 
wilcox.test(table_eccentricity[,2],table_eccentricity[,8],paired = TRUE) # 
wilcox.test(table_eccentricity[,5],table_eccentricity[,8],paired = TRUE) # 

wilcox.test(table_eccentricity[,3],table_eccentricity[,6],paired = TRUE) #
wilcox.test(table_eccentricity[,3],table_eccentricity[,9],paired = TRUE) # 
wilcox.test(table_eccentricity[,6],table_eccentricity[,9],paired = TRUE) # 

# eigenvector centrality (modest)
# wilcox.test(table_eigen[,1],table_eigen[,4],paired = TRUE) # V = 23, p-value = 0.5286
# wilcox.test(table_eigen[,1],table_eigen[,7],paired = TRUE) # V = 55, p-value = 0.05545
# wilcox.test(table_eigen[,4],table_eigen[,7],paired = TRUE) # V = 47.5, p-value = 0.213
# 
# wilcox.test(table_eigen[,2],table_eigen[,5],paired = TRUE) # V = 2, p-value = 0.00293
# wilcox.test(table_eigen[,2],table_eigen[,8],paired = TRUE) # V = 3, p-value = 0.008686
# wilcox.test(table_eigen[,5],table_eigen[,8],paired = TRUE) # V = 61, p-value = 0.009766
# 
# wilcox.test(table_eigen[,3],table_eigen[,6],paired = TRUE) # V = 5, p-value = 0.009766
# wilcox.test(table_eigen[,3],table_eigen[,9],paired = TRUE) # V = 26, p-value = 0.5771
# wilcox.test(table_eigen[,6],table_eigen[,9],paired = TRUE) # V = 59, p-value = 0.01855



