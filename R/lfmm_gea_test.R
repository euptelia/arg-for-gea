#Try lfmm for GEA
# install.packages("foreach")
install.packages("RSpectra")
# devtools::install_github("bcm-uga/lfmm")
library(lfmm)
data("example.data")
data("skin.exposure")
# pca
Y <- example.data$genotype
pc <- prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = "pc", ylab = "Variance explained")

X <- example.data$phenotype
# 170 individuals, 1 phenotype
Y <- example.data$genotype
# 170 individuals, 26943 loci

## Fit an LFMM, i.e. compute B, U, V values
m1_lfmm <- lfmm_ridge(Y = Y,
                      X = X,
                      K = 6)

# Perform association test using the fitted model
pv <- lfmm_test(Y = Y,
                X = X,
                lfmm = m1_lfmm,
                calibrate = "gif") # 
pvalues <- pv$calibrated.pvalue
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues),
       xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)


# Manhattan plot
plot(-log10(pvalues),
     pch = 19,
     cex = .2,
     xlab = "SNP",
     ylab = "-Log P",
     col = "grey")
points(example.data$causal.set,
       -log10(pvalues)[example.data$causal.set],
       type = "h",
       col = "blue")
