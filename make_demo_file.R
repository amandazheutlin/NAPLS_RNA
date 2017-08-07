# make demo text file

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2")
invisible(lapply(libs, require, character.only = TRUE))

#### load in data
# df = everyone from NAPLS with imaging
# conv = list of converters (N = 94)
df   <- read.table("napls_pheno-demo.txt", header=T, stringsAsFactors = F)
conv <- read.table("ConversionAndPopsCriteria_13Jan2016.csv", header=T, sep=",", stringsAsFactors = F)

# make unique subject IDs
df$SiteSubjID   <- paste0(df$SiteNumber,df$SubjectNumber) %>% as.integer()
conv$SiteSubjID <- paste0(conv$SiteNumber,conv$SubjectNumber) %>% as.integer()

# conversion variable
df$dx <- df$SubjectType
df$dx <- ifelse(df$SiteSubjID %in% conv$SiteSubjID, "Converter", df$dx)
df$dx <- as.factor(df$dx)
levels(df$dx) <- c("Control", "Converter", "NonConverter")

write.table(df[,c(1:12,23)],"napls_demo.txt",sep="\t",
            quote=F,row.names=F,col.names=T)


