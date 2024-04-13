# Make a random map for SLiM
map1 <- matrix(runif(100, 0, 1), nrow=10)
write.table(map1, file = "/home/tianlin/Documents/github/data/slim_data/map/random_map1.csv",
          quote = F, row.names = F, col.names = F, sep=",")
