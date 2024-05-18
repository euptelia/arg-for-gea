#Try lfmm2 for GEA, using package LEA
# Following instructions from
# http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf
# install.packages("devtools")
# devtools::install_github("bcm-uga/LEA")
library(LEA)
setwd("/home/tianlin/Documents/github/data/R/")

#### PCA ####
# # Memory error for the whole dataset: Use subsampling instead
# ftemp <- read.table("Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1464742874.lfmm")
# pc = pca("Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1464742874.lfmm", 
#          scale = T)

# Sample 10000 SNPs 
genoFileName = "Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1707743957.geno"
ftemp <- read.table(genoFileName,
                    colClasses = "character")
lf <- read.table("Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1707743957_LFsite.txt")

sampleSNP <- sort(sample(seq(nrow(ftemp)), 10000))
fout <- as.character(ftemp[sampleSNP,])
write.table(fout, 
            file = paste(substr(genoFileName, 1, (nchar(genoFileName)-5)),
                         "10000SNP_sample2.geno",
                         sep = "_"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

pc = pca("Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1707743957_10000SNP.geno", 
         scale = T)
tw = tracy.widom(pc)
# plot the percentage of variance explained by each component
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)



#### STRUCTURE ####
project = NULL
project = snmf("Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1707743957.geno",
               K = 1:20,
               entropy = TRUE,
               repetitions = 5,
               project = "new",
               CPU = 8)

#### lfmm2 ####
# data("offset_example")
# 
# dim(offset_example$geno)
# dim(offset_example$env)
# dim(offset_example$env.pred)

# 4604 diploid individuals with 305435 biallelic SNPs
Y <- read.geno("Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1707743957_10000SNP.geno")
# 1 environmental variables
X <- read.env("Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1707743957.env")

# Took 40%-60% memory for the complete dataset (25-38 Gb)
mod.lfmm2 <- lfmm2(input = Y, env = X, K = 12)

# showing the estimated factors
par(mfrow = c(1,1))
plot(mod.lfmm2@U, col = "grey", pch = 19,
     xlab = "Factor 1",
     ylab = "Factor 2")
pv <- lfmm2.test(object = mod.lfmm2,
                 input = Y,
                 env = X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

