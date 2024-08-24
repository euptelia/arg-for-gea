#Try lfmm2 for GEA, using package LEA
# Following instructions from
# http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf
# install.packages("devtools")
# devtools::install_github("bcm-uga/LEA")
library(LEA)
library(ggplot2)
setwd("/home/tianlin/Documents/github/data/R/20240507/")
# sample_size = 500
#### PCA ####
# # Memory error for the whole dataset
# # Use 10000 randomly chosen SNPs
# ftemp <- read.table("Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1464742874.lfmm")
# pc = pca("Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed1464742874.lfmm", 
#          scale = T)

# runName = "Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.06_mateD0.15_seed101863569341019597_tick110000_mutSeed1191756589"
# shortName="LowPoly_HighMig_ClineMap"

# runName = "Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.03_mateD0.12_seed166332697017507196_tick110000_mutSeed1297446955"
# shortName="HighPoly_LowMig_ClineMap"

# runName = "Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed147836278104916407_tick110000_mutSeed1594694943"
# shortName="HighPoly_HighMig_ClineMap"

# runName = "Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.03_mateD0.12_seed110027493345619720_tick110000_mutSeed1253419139"
# shortName="LowPoly_LowMig_ClineMap"

# runName = "Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.03_mateD0.12_seed1464898069755111446_tick110000_mutSeed1426871279"
# shortName = "LowPoly_LowMig_patchyMap"

# runName = "Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.03_mateD0.12_seed2050996937094533611_tick110000_mutSeed2065548675"
# shortName = "HighPoly_LowMig_patchyMap"

# runName = "Continuous_nonWF_M2b_glacialHistory_patchyMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed15101314668201809_tick110000_mutSeed166593282"
# shortName = "HighPoly_highMig_patchyMap"

# runName = "Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.06_mateD0.15_seed101863569341019597_tick110000_mutSeed1191756589"
# shortName = "LowPoly_highMig_patchyMap"

# runName = "Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.03_mateD0.12_seed166332697017507196_tick101000_mutSeed964677274"
# shortName = "HighPoly_LowMig_ClineMap_101000"

# runName = "Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-08_sigmaM0.01_sigmaW0.4_sigmaD0.06_mateD0.15_seed147836278104916407_tick101000_mutSeed915207499"
# shortName = "HighPoly_HighMig_ClineMap_101000"

# runName = "Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.03_mateD0.12_seed110027493345619720_tick101000_mutSeed495925480"
# shortName = "LowPoly_LowMig_ClineMap_101000"

runName = "Continuous_nonWF_M2b_glacialHistoryOptimum0_clineMap_mu1.0e-09_sigmaM0.1_sigmaW0.4_sigmaD0.06_mateD0.15_seed101863569341019597_tick101000_mutSeed116934221"
shortName = "LowPoly_HighMig_ClineMap_101000"

genoFileName = paste(runName, "_500Inds.geno", sep="")
envFileName =  paste(runName, "_500Inds.env", sep="")
lfFileName =  paste(runName, "_LFsite_500Inds.txt", sep="")
# ftemp <- read.table(genoFileName,
#                     colClasses = "character")
# # Sample 10000 SNPs for PCA
# sampleSNP <- sort(sample(seq(nrow(ftemp)), 10000))
# fout <- as.character(ftemp[sampleSNP,])

# # Sample 500 diploid individuals
# N <- length(ftemp[[1]])
# fout <- data.frame(matrix(NA, nrow = N, ncol = 1))
# sampleInd <- sort(sample(seq(nchar(ftemp[[1]][1])), sample_size))
# for (i in seq(N)){
#   fout[i, 1] <- paste(substring(ftemp[[1]][i], sampleInd, sampleInd), 
#                          collapse = '')
# }
# # Remove non-polymorphic sites
# fout <- fout[fout[1]]
# sample_name <- paste(substr(genoFileName, 1, (nchar(genoFileName)-5)),
#                      "500Inds.geno",
#                      sep = "_")
# write.table(fout, 
#             file = sample_name,
#             row.names = FALSE,
#             col.names = FALSE,
#             quote = FALSE)

pc = pca(genoFileName, 
         scale = T)
tw = tracy.widom(pc) 
tw
# plot the percentage of variance explained by each component
setwd("/home/tianlin/Documents/github/data/R/202408/figures")
png(file = paste(runName, "_pca_percentage.png", sep = ""),
    height = 800,
    width = 1200,
    res = 150)
plot(tw$percentage, pch = 19, col = "grey42", cex = .8)
dev.off()

#### STRUCTURE ####
# # Takes a long time
# project = NULL
# project = snmf(genoFileName,
#                K = 1:20,
#                entropy = TRUE,
#                repetitions = 5,
#                project = "new",
#                CPU = 8)
# snmfName = paste(runName, 
#                  "_500Inds.snmfProject",
#                  sep = "")
# project = load.snmfProject(file=snmfName)
# 
# png(file = paste(runName, "_snmf.png", sep = ""),
#     height = 800,
#     width = 1200,
#     res = 150)
# plot(project, col = "grey42", pch = 19, cex = 1.2,
#      main = shortName)
# dev.off()

#### lfmm2 ####
# data("offset_example")
# 
# dim(offset_example$geno)
# dim(offset_example$env)
# dim(offset_example$env.pred)

setwd("/home/tianlin/Documents/github/data/R/20240507/")
# diploid individuals with biallelic SNPs
Y <- read.geno(genoFileName)
# 1 environmental variables
X <- read.env(envFileName)
#lf_site
lf <- read.table(lfFileName)
# real_positive <- as.factor(relative_lf > 0.01)


# Took 40%-60% memory for the complete dataset (25-38 Gb)
# 7% memory for 500-individual sample
K <- 4
# K <- 8
mod.lfmm2 <- lfmm2(input = Y, env = X, K = K)

# showing the estimated factors
par(mfrow = c(1,1))
plot(mod.lfmm2@U, col = "grey", pch = 19,
     xlab = "Factor 1",
     ylab = "Factor 2")
pv <- lfmm2.test(object = mod.lfmm2,
                 input = Y,
                 env = X,
                 full = TRUE)
# plot(-log10(pv$pvalues), 
#      col = c("grey", "firebrick")[real_positive], 
#      cex = c(0.5, 1)[real_positive], 
#      pch = c(19, 17)[real_positive])
# abline(h = -log10(0.1/510), lty = 2, col = "orange")

df <- data.frame(pvalues=-log10(pv$pvalues),
                 lf=lf[1],
                 relative_lf=(lf/sum(lf))[1],
                 real_positive=as.factor((lf/sum(lf))[1] > 0.01),
                 temp_id = seq(length(pv$pvalues)))
colnames(df) <- c("pvalues", "lf", "relative_lf", "real_positive", "temp_id")
df_realPositive = df[df$real_positive==TRUE,]


# # Method1
# color2 <- c("grey", "firebrick")      
# shape2 <- c(19,17)
# ggplot(df, aes(temp_id, pvalues, 
#                color=real_positive,
#                shape=real_positive,
#                size=real_positive,
#                #alpha=c(0.5,1)[real_positive]
#                ))+
#   geom_point()+
#   scale_shape_manual(values=c(19,17))+
#   scale_size_manual(values=c(1,4))+
#   scale_color_manual(values=color2)+
#   theme_bw()
setwd("/home/tianlin/Documents/github/data/R/20240507/figures")
png(file = paste(runName, "_lfmm2_K", K, "_pValues_and_lfsite.png", sep = ""),
    height = 800,
    width = 1600,
    res = 150)
color2 <- c("grey", "firebrick")      
shape2 <- c(19,17)
ggplot(df, aes(temp_id, pvalues))+
  geom_point(aes(), 
             colour="grey", size=1, shape=19)+
  geom_point(data=df_realPositive,
             aes(x=temp_id,
                 y=pvalues),
             colour="firebrick", size=4, shape=17, alpha=0.4)+
  geom_text(data=df_realPositive,
            aes(x=temp_id,
                y=(pvalues+max(df$pvalues)/20),
                label=round(relative_lf,2)),
            colour="firebrick",
            alpha=0.8)+
  labs(title=shortName,
        x ="SNP ID", y = "-log10(p-value)")+
# scale_shape_manual(values=c(19,17))+
# scale_size_manual(values=c(0.5,1))+
# scale_color_manual(values=color2)+
theme_bw()
dev.off()
