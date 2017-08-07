# RNA-seq data from NAPLS - miRNAs
# correlations with gray matter longitudinal phenotype
# plasticity-related miRNA

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/RNAseq/"
setwd(workdir)

# source("http://bioconductor.org/biocLite.R")
libs <- c("dplyr", "psych", "ggplot2", "CALF", "glmnet", "caret", 
          "nlme", "car", "mediation")
invisible(lapply(libs, require, character.only = TRUE))

#### load in data
# 136 miRNAs; 3 groups (NP, P, UC)
data  <- read.table("miRNA_norm.txt",header=T)
gm    <- read.table("GMpheno.txt",header=T)
covar <- read.table("age-sex.txt",header=T)

# data$cc <- data$ProdromeStatus
# levels(data$cc) <- c(0,1)
gm$SiteSubjID    <- paste(gm$SiteNumber,gm$SubjectNumber,sep="") %>% as.integer()
covar$SiteSubjID <- paste(covar$SiteNumber,covar$SubjectNumber,sep="") %>% as.integer()
data$gm          <- gm[match(data$SiteSubjID,gm$SiteSubjID),3]
data$age         <- covar[match(data$SiteSubjID,covar$SiteSubjID),1]
data$sex         <- covar[match(data$SiteSubjID,covar$SiteSubjID),2]

# review: Mellios & Sur, 2012, Frontiers in Psychiatry
# miR-21 (2), miR-30a (1), miR-30d (1), miR-34a (0), miR-128 (1), 
# miR-132 (1), miR-134 (0), miR-137 (0), miR-138 (0), miR-181b (1), 
# miR-195 (0), miR-212 (0), miR-219 (0), miR-346 (0)
miRNA <- c("miR.21.5p", "miR.21.3p", "miR.30a.5p", "miR.30d.5p", 
           "miR.128", "miR.132.3p", "miR.181b.5p")

# multiple regression
m1 <- lm(gm ~ miR.21.5p + miR.21.3p + miR.30a.5p + miR.30d.5p + miR.128 + miR.132.3p + miR.181b.5p,
         data = data, na.action = na.omit)

# univariate analyses
models <- lapply(miRNA, function(x) {
  lm(eval(substitute(gm ~ i, list(i = as.name(x)))),data = data, na.action = na.omit)
})

# stats
model_stats <- lapply(models, function(x) summary(x))

# list of overall effects
gm = NULL
for (i in 1:7) {
  temp    <- model_stats[[i]]$coefficients %>% as.data.frame()
  temp$R2 <- model_stats[[i]]$r.squared
  gm      <- rbind(gm,temp[2,3:5])
}

colnames(gm)[1:2] <- c("tvalue","pvalue")
gm$marker         <- rownames(gm)
gm                <- gm[order(gm$pvalue),]
