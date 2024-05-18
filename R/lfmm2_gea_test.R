#Try lfmm2 for GEA, using package LEA
# Following instructions from
# http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf
# install.packages("devtools")
# devtools::install_github("bcm-uga/LEA")
library(LEA)
setwd("/home/tianlin/Documents/github/data/R/")
data("tutorial")
# The data include 400 SNPs for 50 individuals.
dim(tutorial.R)
dim(tutorial.C)
write.lfmm(tutorial.R, "genotypes.lfmm")
write.geno(tutorial.R, "genotypes.geno")
write.env(tutorial.C, "gradients.env")

#### PCA ####
pc = pca("genotypes.lfmm", scale = T)
# test if the eigenvectors provide evidence for a population structure
tw = tracy.widom(pc) 
tw


#### lfmm ####
project = NULL
project = lfmm("genotypes.lfmm",
               "gradients.env",
               K = 6,
               repetitions = 5,
               project = "new")
# Combineing z-scores obtained from multile runs and adjusted p-values
p = lfmm.pvalues(project, K = 6)
pvalues = p$pvalues
pvalues
# Plot pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "seagreen3")
plot(-log10(pvalues), 
     pch = 19, col = "seagreen", 
     cex = .7)


#### lfmm2 ####
data("offset_example")

dim(offset_example$geno)
dim(offset_example$env)
dim(offset_example$env.pred)

# 200 diploid individuals genotyped at 510 SNP
Y <- offset_example$geno
# 4 environmental variables
X <- offset_example$env

mod.lfmm2 <- lfmm2(input = Y, env = X, K = 2)

# showing the K = 2 estimated factors
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
