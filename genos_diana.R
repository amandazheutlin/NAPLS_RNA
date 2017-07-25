# extracting genotypes for diana

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/"
setwd(workdir)

libs <- c("dplyr", "ggplot2", "gdata")
invisible(lapply(libs, require, character.only = TRUE))

#### data
# read in genotypes
genos <- read.csv("genos_diana.ped", sep="\t", header=F) # genotypes
snps  <- read.table("genos_diana.map", header=F) # snp IDs

# remove spaces + name cols
genos1 <- lapply(genos, function(x) gsub(" ", "", x, fixed = T)) %>% as.data.frame
genos1 <- genos1[,c(1:2,7:31)] # remove excess cols

snp.names     <- snps[,2] %>% as.character 
names(genos1) <- c("FID", "IID", snp.names)

# convert geno IDs to site/subj IDs
id.conv <- read.xls("napls2_IDconversion.xls", stringsAsFactors = F)

genos1$siteID <- id.conv[match(genos1$IID,id.conv$Collaborator.Participant.ID),9]
genos1$subjID <- id.conv[match(genos1$IID,id.conv$Collaborator.Participant.ID),10]

# write genos
write.table(genos1[,c(28,29,3:27)],"genos_diana_cleaned.txt",
            sep="\t", col.names=T, row.names=F, quote=F)

