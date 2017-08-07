# miRNA classifier extension (built in NAPLS originally)
# correlations between miRNA classifier mRNA targets + cortical thickness
# current data = swedish twins

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/mRNA_GM/"
setwd(workdir)

source("https://bioconductor.org/biocLite.R")
biocLite("illuminaHumanv4.db")
libs <- c("dplyr", "psych", "ggplot2", "nlme", "illuminaHumanv4.db")
invisible(lapply(libs, require, character.only = TRUE))

#### load data
# add updated dx info
load("/data/swe_gwas/ABZ/RNA_GWAS/swedenclean.rdata")
dx = read.table("/data/swe_gwas/ABZ/RNA_GWAS/swedishNP_PCA.txt",header=T)
swedenclean$dx_new <- dx[match(swedenclean$StudyID,dx$IID),7] %>% as.factor()

# zygosity + tvab
zyg              <- read.table("/data/swe_gwas/ABZ/immunegenes/zygosity_alltwins.txt",header=T)
swedenclean$tvab <- zyg[match(swedenclean$StudyID,zyg$IID),4]
swedenclean$zyg  <- zyg[match(swedenclean$StudyID,zyg$IID),3]

# phenotype (superior frontal cortical thickness)
pheno <- read.table("/data/swe_gwas/ABZ/immunegenes/supfrontal_thickness.txt",header=T)
pheno$IID <- substr(pheno$ID,1,10)
swedenclean$SFthick <- pheno[match(swedenclean$StudyID,pheno$IID),"supfront"]

# gene targets of miRNA classifier
erk5 <- read.table("ERK5-genes.txt",header=T,sep="\t")
pka  <- read.table("PKA-genes.txt",header=T,sep="\t")
fpcv <- read.table("FPCV-genes.txt",header=T,sep="\t")
ppar <- read.table("PPAR-genes.txt",header=T,sep="\t")
hgf  <- read.table("HGF-genes.txt",header=T,sep="\t")
reel <- read.table("reelin-genes.txt",header=T,sep="\t")

#### correlations between mRNA targets and cortical thickness
xx        <- as.data.frame(illuminaHumanv4ALIAS2PROBE) # all the genes

## ERK5 signaling
erk5.list <- erk5$EntrezGeneName %>% as.character # 28 genes
erk5.xx   <- xx[xx$alias_symbol %in% erk5.list,] # genes + probe IDs
probes    <- colnames(swedenclean)[colnames(swedenclean) %in% erk5.xx$probe_id]
vars      <- c("StudyID","Family","Sex","Age","Diagnosis","SFthick",probes)

df.erk5 <- swedenclean[,names(swedenclean) %in% vars]

# univariate
RNA.erk5 <- matrix(nrow=length(probes),ncol=4)
for (i in 1:length(probes)){
  f <- formula(paste("SFthick ~", probes[i],"+ Age + Sex + Diagnosis", sep=""))
  RNA.m2 <- lme(fixed=f, random = ~1|Family, data = df.erk5, na.action = na.omit)
  tval   <- summary(RNA.m2)$tTable[2,4]
  RNA.erk5[i,1] <- probes[i]
  RNA.erk5[i,2] <- tval
  RNA.erk5[i,3] <- sqrt(tval^2/(tval^2+summary(RNA.m2)$tTable[2,3]))
  RNA.erk5[i,4] <- summary(RNA.m2)$tTable[2,5]
}

RNA.erk5        <- as.data.frame(RNA.erk5)
names(RNA.erk5) <- c("item","t","corr","p")
RNA.erk5$gene   <- erk5.xx[match(RNA.erk5$item,erk5.xx$probe_id),2]
RNA.erk5        <- RNA.erk5[order(RNA.erk5$p),]

# sum score
df.erk5$erk5 <- dplyr::select(df.erk5,one_of(probes)) %>% rowSums()
m1           <- lme(SFthick ~ erk5 + Age + Sex + Diagnosis, random = ~1|Family, data = df.erk5, na.action = na.omit)
summary(m1) # t = .614, p = .542 (huge effect of age)

## PKA signaling
pka.list  <- pka$Symbol %>% as.character # 98 genes
pka.xx    <- xx[xx$alias_symbol %in% pka.list,] # genes + probe IDs
probes    <- colnames(swedenclean)[colnames(swedenclean) %in% pka.xx$probe_id]
vars      <- c("StudyID","Family","Sex","Age","Diagnosis","SFthick",probes)

df.pka <- swedenclean[,names(swedenclean) %in% vars]

# sum score
df.pka$pka <- dplyr::select(df.pka,one_of(probes)) %>% rowSums()
m2         <- lme(SFthick ~ pka + Age + Sex + Diagnosis, random = ~1|Family, data = df.pka, na.action = na.omit)
summary(m2) # t = .921, p = .362 (huge effect of age)

## factors promoting cardiogenesis in vertebrates
fpcv.list <- fpcv$Symbol %>% as.character # 34 genes
fpcv.xx   <- xx[xx$alias_symbol %in% fpcv.list,] # genes + probe IDs
probes    <- colnames(swedenclean)[colnames(swedenclean) %in% fpcv.xx$probe_id]
vars      <- c("StudyID","Family","Sex","Age","Diagnosis","SFthick",probes)

df.fpcv <- swedenclean[,names(swedenclean) %in% vars]

# sum score
df.fpcv$fpcv    <- dplyr::select(df.fpcv,one_of(probes)) %>% rowSums()
m3  <- lme(SFthick ~ fpcv + Age + Sex + Diagnosis, random = ~1|Family, data = df.fpcv, na.action = na.omit)
summary(m3) # t = .123, p = .903 (huge effect of age)

## PPARa/RXRa activation
ppar.list <- ppar$Symbol %>% as.character # 50 genes
ppar.xx   <- xx[xx$alias_symbol %in% ppar.list,] # genes + probe IDs
probes    <- colnames(swedenclean)[colnames(swedenclean) %in% ppar.xx$probe_id]
vars      <- c("StudyID","Family","Sex","Age","Diagnosis","SFthick",probes)

df.ppar <- swedenclean[,names(swedenclean) %in% vars]

# sum score
df.ppar$ppar    <- dplyr::select(df.ppar,one_of(probes)) %>% rowSums()
m4  <- lme(SFthick ~ ppar + Age + Sex + Diagnosis, random = ~1|Family, data = df.ppar, na.action = na.omit)
summary(m4) # t = .629, p = .533 (huge effect of age)

## HGF signaling
hgf.list <- hgf$Symbol %>% as.character # 43 genes
hgf.xx   <- xx[xx$alias_symbol %in% hgf.list,] # genes + probe IDs
probes    <- colnames(swedenclean)[colnames(swedenclean) %in% hgf.xx$probe_id]
vars      <- c("StudyID","Family","Sex","Age","Diagnosis","SFthick",probes)

df.hgf <- swedenclean[,names(swedenclean) %in% vars]

# sum score
df.hgf$hgf    <- dplyr::select(df.hgf,one_of(probes)) %>% rowSums()
m5  <- lme(SFthick ~ hgf + Age + Sex + Diagnosis, random = ~1|Family, data = df.hgf, na.action = na.omit)
summary(m5) # t = .831, p = .410 (huge effect of age)

## Reeling signaling in neurons
reel.list <- reel$Symbol %>% as.character # 32 genes
reel.xx   <- xx[xx$alias_symbol %in% reel.list,] # genes + probe IDs
probes    <- colnames(swedenclean)[colnames(swedenclean) %in% reel.xx$probe_id]
vars      <- c("StudyID","Family","Sex","Age","Diagnosis","SFthick",probes)

df.reel <- swedenclean[,names(swedenclean) %in% vars]

# sum score
df.reel$reel    <- dplyr::select(df.reel,one_of(probes)) %>% rowSums()
m6  <- lme(SFthick ~ reel + Age + Sex + Diagnosis, random = ~1|Family, data = df.reel, na.action = na.omit)
summary(m6) # t = .526, p = .602 (huge effect of age)

##### sum across genes
allgenes    <- c(pka.list, erk5.list, fpcv.list, ppar.list, hgf.list, reel.list) # 285 genes
allgenes.xx <- xx[xx$alias_symbol %in% allgenes,] # genes + probe IDs
probes      <- colnames(swedenclean)[colnames(swedenclean) %in% allgenes.xx$probe_id]
vars        <- c("StudyID","Family","Sex","Age","Diagnosis","SFthick",probes)

df.allgenes <- swedenclean[,names(swedenclean) %in% vars]

# sum score
df.allgenes$allgenes  <- dplyr::select(df.allgenes,one_of(probes)) %>% rowSums()
m7  <- lme(SFthick ~ allgenes + Age + Sex + Diagnosis, random = ~1|Family, data = df.allgenes, na.action = na.omit)
summary(m7) # t = .485, p = .630 (huge effect of age)
