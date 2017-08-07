# PRS x miRNA classifier

#### housekeeping
RNAdir <- "/data/swe_gwas/ABZ/NAPLS/RNAseq/"
PRSdir <- "/data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/"
setwd("/data/swe_gwas/ABZ/NAPLS/")

libs <- c("dplyr", "ggplot2", "gdata")
invisible(lapply(libs, require, character.only = TRUE))

#### load in miRNA data
# 136 miRNAs; 3 groups (NP, P, UC)
data  <- read.table(paste0(RNAdir, "miRNA_norm.txt"), header=T)
gm    <- read.table(paste0(RNAdir, "GMpheno.txt"), header=T)
covar <- read.table(paste0(RNAdir, "age-sex.txt"), header=T)

gm$SiteSubjID    <- paste(gm$SiteNumber,gm$SubjectNumber,sep="") %>% as.integer()
covar$SiteSubjID <- paste(covar$SiteNumber,covar$SubjectNumber,sep="") %>% as.integer()
data$gm          <- gm[match(data$SiteSubjID,gm$SiteSubjID),3]
data$age         <- covar[match(data$SiteSubjID,covar$SiteSubjID),1]
data$sex         <- covar[match(data$SiteSubjID,covar$SiteSubjID),2]

# miRNA classifier based on 
# past subsetting without rerunning
# read in; don't rerun CALF
hits <- read.table(paste0(RNAdir, "final_hits_miRNA_9.txt"), stringsAsFactors = F)
hits <- hits[,1]
hits <- gsub("hsa-","", hits)
hits <- gsub("-",".", hits)
weights  <- c(-1,-1,1,-1,1,-1,1,1,-1)
hits9    <- hits
weights9 <- weights

data$group2         <- data$group
levels(data$group2) <- c("not","P","not") 
levels(data$group2) <- c("Nonconverter or Control", "Converter")
levels(data$group)  <- c("Nonconverter","Converter", "Control")

data$gm.score.9 <- dplyr::select(data,one_of(hits9)) %>% mapply("*",.,weights9) %>% rowSums()
data$group      <- factor(data$group,levels(data$group)[c(2,1,3)])

data.inc <- data[complete.cases(data$gm),]

######## load in PRS data
# genetics PCs from MDS in plink (see plink_PCA_cauc.sh)
PC.df <- read.table(paste0(PRSdir,"pca.pruned.data/napls_cauc_PCs.eigenvec"),header=T,stringsAsFactors = F)

# SZ risk scores
# read in all the score files (8 pT)
# score input lists from preprocessing.R
# scores generated using prs/SZrisk-prs.sh by plink
score_files   <- list.files(path="PSYCH-CHIP/prs",pattern="napls_SZrisk_.*\\.profile")
score_files   <- lapply(score_files, function(x) gsub("napls","PSYCH-CHIP/prs/napls",x)) %>% unlist
scores        <- lapply(score_files, function(x) read.table(x, header=T, stringsAsFactors=F))
names(scores) <- c("pT001","pT01","pT05","pT1","pT2","pT5","pTall","sigsnps")

# just select the score columns
prs.df <- lapply(scores, "[[", 4) %>% as.data.frame
prs.df <- cbind(scores$pT001$FID, scores$pT001$IID, scores$pT001$PHENO, prs.df)
names(prs.df)[1:3] <- c("FID","IID","PHENO")


######## combine to test PRS x miRNA classifier
# convert between IDs
# gene IDs --> imaging / demo IDs
# EDITED XLS: for UNC subjects, Collaborator.Sample.ID_2 seems like 
#             it is site and subject number so manually edited the
#             corresponding columns -- NAPLS2.Site.Number + NAPLS2.Site.ID
# PCs are from caucasian only, so subj loss here
id.conv <- read.xls(paste0(PRSdir,"NAPLS2_IDconversion_edited.xls"), stringsAsFactors = F)
prs.PC  <- merge(PC.df,prs.df,by="IID")

prs.PC$site <- id.conv[match(prs.PC$IID, id.conv$Collaborator.Sample.ID), "NAPLS2.Site.Number"]
prs.PC$subj <- id.conv[match(prs.PC$IID, id.conv$Collaborator.Sample.ID), "NAPLS2.Site.ID"]
prs.PC$SiteSubjID <- paste0(prs.PC$site,prs.PC$subj)

data.all <- merge(data.inc, prs.PC[,c("SiteSubjID", "pTall", "PC1", "PC2", "PC3")], by="SiteSubjID")

# stats (N = 53)
m1 <- lm(pTall ~ gm.score.9 + group + age + sex + PC1 + PC2 + PC3, data = data.all)
summary(m1)

