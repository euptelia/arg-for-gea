# Mean of fst of each combination
df <- read.table(file = "/home/tianlin/Documents/github/data/tskit_data/figure/20240426/fst20240426.txt",
                 sep = ",",
                 header = F)
df2 <- data.frame(model=apply(df[df$V12 == "tick110000", 3:10], 1, paste, collapse="_"),
                  fst=df[df$V12 == "tick110000",14])
tapply(df2$fst, INDEX=as.factor(df2$model), FUN = mean)
 