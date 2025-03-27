# Make a random map for SLiM
map1 <- matrix(runif(100, 0, 1), nrow=10)
write.table(map1, file = "/home/tianlin/Documents/github/data/slim_data/map/random_map1.csv",
          quote = F, row.names = F, col.names = F, sep=",")
          
#Make ten other maps for recurrent changes in M3a (Not used)
for (i in 1:10){
  map1 <- matrix(runif(100, 0, 1), nrow=10)
  outName <- paste(paste("/home/anadem/github/arg-for-gea/map/random_historical_map", 
                         i, sep = ""), ".csv", sep = "")
  write.table(map1, file = outName,
              quote = F, row.names = F, col.names = F, sep=",")
}
