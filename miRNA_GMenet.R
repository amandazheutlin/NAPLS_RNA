# elastic net of miRNA and gray matter loss
# NAPLS pilot data (N = 74)
# functional annotation of selected miRNA

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/RNAseq/"
setwd(workdir)
set.seed(1)

# source("http://bioconductor.org/biocLite.R")
libs <- c("dplyr", "psych", "ggplot2", "CALF", "glmnet", "caret", 
          "nlme", "car", "illuminaHumanv4.db")
invisible(lapply(libs, require, character.only = TRUE))

#### load in data
# 136 miRNAs; 3 groups (NP, P, UC)
data  <- read.table("miRNA_norm.txt",header=T)
gm    <- read.table("GMpheno.txt",header=T)

gm$SiteSubjID    <- paste(gm$SiteNumber,gm$SubjectNumber,sep="") %>% as.integer()
data$gm          <- gm[match(data$SiteSubjID,gm$SiteSubjID),3]

#### LASSO against GM pheno
miRNA      <- names(data)[6:141]
df         <- data[complete.cases(data$gm),]
predictors <- df %>% dplyr::select(one_of(miRNA)) %>% scale() 
target     <- df[,142]

fitControl <- trainControl(method = "LOOCV", savePredictions = TRUE)
eGrid      <- expand.grid(.alpha = seq(1), .lambda = seq(0,1,by=0.05))
enet       <- train(x= as.matrix(predictors), y = as.numeric(target),
                    method = "glmnet", 
                    metric = "Rsquared",
                    trControl = fitControl)
enet

# variable importance
vars         <- varImp(enet)$importance %>% as.data.frame 
vars$markers <- rownames(vars)
vars         <- vars[order(-vars$Overall),]
vars$markers <- factor(vars$markers,
                       levels = vars$markers[order(vars$Overall)])

# permutations
perm.df <- NULL
for (i in 1:1000){
  target.rand <- sample(target)
  enet.rand   <- train(x= as.matrix(predictors), y = as.numeric(target.rand),
                      method = "glmnet", 
                      metric = "Rsquared",
                      trControl = fitControl)
  parameters  <- enet.rand$results   
  perm.df     <- rbind(perm.df,parameters[parameters$Rsquared==max(parameters$Rsquared),])
}

perm.df        <- as.data.frame(perm.df)
names(perm.df) <- names(parameters)
perm.df        <- perm.df[order(perm.df$Rsquared),]
perm.better    <- perm.df[perm.df$Rsquared > 0.096569564,] # (143 + 1)/(868 + 1) = .166

# graph it
ggplot(vars[1:20,], aes(x=markers,y=Overall)) + 
  geom_bar(stat="identity",width=.85, fill="#33CCFF") +
  geom_text(stat="identity",y=22, hjust=0, 
            aes(label=markers), colour="black", size=5) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13,color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept=60, colour="red", linetype="longdash") +
  ylab("Variable Importance") +
  xlab("miRNAs") +
  coord_flip(ylim = c(20,100))

###### functional annotations
miranda <- read.csv("human_predictions_S_C_aug2010.txt",
                    header=T, stringsAsFactors = F, sep="\t")

# lists of miRNAs to classify
# match miranda formatting
vars$markers.m <- gsub(".","-",vars$markers,fixed = T)
vars$markers.m <- paste0("hsa-",vars$markers.m)
vars$markers.c <- lapply(strsplit(vars$markers.m,"-"), function(x) paste(x[1:3],collapse="-")) 

top3    <- vars[1:3,4] %>% as.character()
above60 <- vars[vars$Overall>60,4] %>% as.character()

# remove miRNAs that are not predictors
# remove gene targets with mirsvr scores < -.5
top3.m <- miranda[miranda$mirna_name %in% top3,] # uses 2/3 of them (9529 rows)
top3.m <- top3.m[top3.m$mirsvr_score < -.1,] # 9527 rows

above60.m <- miranda[miranda$mirna_name %in% above60,] # uses 5/7 of them (21097 rows)
above60.m <- above60.m[above60.m$mirsvr_score < -.1,] # 21093 rows

convscore <- c("hsa-miR-941", "hsa-miR-103a", "hsa-miR-199a", "hsa-miR-92a", "hsa-miR-31")
convscore.m <- miranda[miranda$mirna_name %in% convscore,]

# use table to determine how many times each gene is listed
# discard either non-duplicates or all not listed X times
genefreq <- table(top3.m$gene_symbol) %>% as.data.frame()
genes    <- genefreq[genefreq$Freq > 1,]

genefreq2 <- table(above60.m$gene_symbol) %>% as.data.frame()
genes2a   <- genefreq2[genefreq2$Freq > 1,]
genes2b   <- genefreq2[genefreq2$Freq > 4,]

genefreq3 <- table(convscore.m$gene_symbol) %>% as.data.frame()
genes3    <- genefreq3[genefreq3$Freq > 1,]

# use intersection list for functional annotation
write.table(genes[,1],"genelists/top3.miranda.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(genes2a[,1],"genelists/top7more.miranda.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(genes2b[,1],"genelists/top7less.miranda.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(genes3[,1],"genelists/convscore.miranda.txt",col.names=F,row.names=F,quote=F,sep="\t")

# overlap in miRNA genes + exome regions
exome   <- read.table("../exome_gwas/150genes_GM.txt",header=F,stringsAsFactors = F)
overlap <- genes2b[genes2b$Var1 %in% exome$V1,] # 6 genes

#### check genes in swedish sample
load("/data/swe_gwas/ABZ/RNA_GWAS/swedenclean.rdata")
dx = read.table("/data/swe_gwas/ABZ/RNA_GWAS/swedishNP_PCA.txt",header=T)
swedenclean$dx_new <- dx[match(swedenclean$StudyID,dx$IID),7] %>% as.factor()

xx <- as.data.frame(illuminaHumanv4ALIAS2PROBE) # all the genes
zz <- xx[xx$probe_id %in% colnames(swedenclean),] # all the swedish probes
aa <- zz[zz$alias_symbol %in% overlap$Var1,]

columns <- c("StudyID","Family","dx_new", aa$probe_id)
df      <- swedenclean %>% dplyr::select(one_of(columns)) 

pheno     <- read.table("/data/swe_gwas/ABZ/MRI_ROI_swedes/FS_VolMeanThick_20151120.txt",header=T)
pheno$IID <- substr(pheno$scanid,1,10)
df$lhThickness <- pheno[match(df$StudyID,pheno$IID),20]
df$rhThickness <- pheno[match(df$StudyID,pheno$IID),21]

# linear mixed models 
genesmod.lme.l  <- lme(lhThickness ~ ILMN_1664750 + ILMN_1704369 + ILMN_1730223 + ILMN_1771689 + ILMN_1788149 + ILMN_2365223 + ILMN_2404746,
                       random = ~1|Family, data = df, na.action = na.omit)
genesmod.lme.r  <- lme(rhThickness ~ ILMN_1664750 + ILMN_1704369 + ILMN_1730223 + ILMN_1771689 + ILMN_1788149 + ILMN_2365223 + ILMN_2404746,
                       random = ~1|Family, data = df, na.action = na.omit)
genesmod.lme.l.sz  <- lme(lhThickness ~ ILMN_1664750 + ILMN_1704369 + ILMN_1730223 + ILMN_1771689 + ILMN_1788149 + ILMN_2365223 + ILMN_2404746,
                       random = ~1|Family, data = cc1, na.action = na.omit)
genesmod.lme.r.sz  <- lme(rhThickness ~ ILMN_1664750 + ILMN_1704369 + ILMN_1730223 + ILMN_1771689 + ILMN_1788149 + ILMN_2365223 + ILMN_2404746,
                       random = ~1|Family, data = cc1, na.action = na.omit)

# case-control differences in gene expression
# univariate models
overlap.genes  <- aa$probe_id
cc             <- subset(df,dx_new==1|dx_new==7)
cc1            <- subset(df,dx_new==1|dx_new==5|dx_new==6|dx_new==7)

cc$dx_new          <- droplevels(cc$dx_new)
cc1$dx_new         <- droplevels(cc1$dx_new)
levels(cc1$dx_new) <- c(1,7,7,7)

genemodels <- lapply(overlap.genes, function(x) {
  lme(eval(substitute(i ~ dx_new, list(i = as.name(x)))),
      random = ~1|Family, data = cc, na.action = na.omit)
})

# stats
genemodel_stats <- lapply(genemodels, function(x) summary(x)$tTable)

# list of overall effects
g.lthick = NULL
for (i in 1:7) {
  temp <- as.data.frame(genemodel_stats[[i]])
  g.lthick <- rbind(g.lthick,temp[2,])
}

colnames(g.lthick)[4:5] <- c("tvalue","pvalue")
g.lthick$marker         <- rownames(g.lthick)
g.lthick                <- g.lthick[order(g.lthick$pvalue),]

# fit a net
df2         <- df[complete.cases(df$lhThickness),]
predictors2 <- df2 %>% dplyr::select(one_of(overlap.genes))  
target2     <- df2$lhThickness

fitControl <- trainControl(method = "LOOCV", savePredictions = TRUE)
eGrid      <- expand.grid(.alpha = seq(1), .lambda = seq(0,1,by=0.05))
enet2      <- train(x= as.matrix(predictors2), y = as.numeric(target2),
                    method = "glmnet", 
                    metric = "Rsquared",
                    trControl = fitControl)
enet2

#### CALF for GM (binary)
# covariates
covar            <- read.table("age-sex.txt",header=T)
covar$SiteSubjID <- paste(covar$SiteNumber,covar$SubjectNumber,sep="") %>% as.integer()
data$age         <- covar[match(data$SiteSubjID,covar$SiteSubjID),1]
data$sex         <- covar[match(data$SiteSubjID,covar$SiteSubjID),2]

# convert gm to binary
data$gm.binary <- ifelse(data$gm < mean(data$gm,na.rm=T),0,1)
calf.df        <- data[complete.cases(data$gm.binary),c(145,6:141)]

# true data
m.true <- calf(data = calf.df, nMarkers = 136) # 9 miRNAs chosen
hits   <- m.true$selection$Marker %>% as.character()

# calf.df.hits <- calf.df %>% dplyr::select(one_of(hits))
# calf.df.hits <- cbind(calf.df$gm.binary,calf.df.hits)

# make a score
data$gm.score <- data$miR.29a.3p + data$miR.93.5p - data$miR.106b.5p + data$miR.193a.5p - data$miR.92b.3p - data$miR.181a.5p + data$miR.3615 - data$miR.181a.2.3p + data$let.7c 
gm.m1         <- lm(gm ~ gm.score + age + sex, data = data, na.action = na.omit)
gm.m2         <- lm(gm ~ gm.score, data = data, na.action = na.omit)
true.roc      <- roc(calf.df[,1],predict.lm(gm.m2)) # AUC = .851

# permutations
m.rand <- calf_randomize(data = calf.df, nMarkers = 136, randomize = T, times = 1000) # mean AUC = .860
m.sub  <- calf_subset(data = calf.df, nMarkers = 9, proportion = .8, times = 1000) # mean AUC = .884



###### empirical GO 
# mature_homosapiens.fa = mature miRNA in humans
# ensembl3utr.txt = sequences to scan
# hsa-ensembl-go.txt = annotations by GO terms
# sequences        <- read.csv("ensembl3utr.txt",header=F,stringsAsFactors = F)
# names(sequences) <- "col1"
# 
# missing <- subset(sequences,col1=="Sequence unavailable")
# dlt     <- row.names(missing) %>% as.numeric
# dlt2    <- dlt-1
# dltall  <- c(dlt,dlt2)
# 
# seqclean <- sequences[!rownames(sequences) %in% dltall,,drop=F]
# write.csv(seqclean,"ensembl3utr_noNA.txt",row.names=F,quote=F)
# 
# seqclean2 <- ifelse(!(substr(sequences$col1,1,1) %in% strings),NULL,sequences$col1)
# 
# # miranda
# mature <- read.csv("mature_homosapiens.fa", header=F, stringsAsFactors = F)
