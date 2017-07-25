# RNA-seq data from NAPLS - miRNAs
# correlations with gray matter longitudinal phenotype

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/RNAseq/"
setwd(workdir)

install.packages("~/Downloads/CALF_0.1.3(1).tar.gz", repos=NULL) # local version of CALF
libs <- c("dplyr", "psych", "ggplot2", "CALF", "nlme", "car", 
          "mediation", "RColorBrewer", "colorRamps", "MASS", "sfsmisc")
invisible(lapply(libs, require, character.only = TRUE))

#### load in data
# 136 miRNAs; 3 groups (NP, P, UC)
data  <- read.table("miRNA_norm.txt",header=T)
gm    <- read.table("GMpheno.txt",header=T)
covar <- read.table("age-sex.txt",header=T)

gm$SiteSubjID    <- paste(gm$SiteNumber,gm$SubjectNumber,sep="") %>% as.integer()
covar$SiteSubjID <- paste(covar$SiteNumber,covar$SubjectNumber,sep="") %>% as.integer()
data$gm          <- gm[match(data$SiteSubjID,gm$SiteSubjID),3]
data$age         <- covar[match(data$SiteSubjID,covar$SiteSubjID),1]
data$sex         <- covar[match(data$SiteSubjID,covar$SiteSubjID),2]

#### CALF of real vector
gm.comp <- data[complete.cases(data$gm),c(142,6:141)] # N = 74
gm.real <- calf(data = gm.comp, nMarkers = 136, targetVector="real")
hits    <- gm.real$selection$Marker %>% as.character()
weights <- gm.real$selection$Weight 

data$gm.score <- dplyr::select(data,one_of(hits)) %>% mapply("*",.,weights) %>% rowSums()
# gm.real.m1    <- lm(gm ~ gm.score + age + sex, data = data, na.action = na.omit)
# gm.real.m2    <- lm(gm ~ gm.score, data = data, na.action = na.omit)
# gm.real.m3    <- lm(gm ~ gm.score + age + sex + group, data = data, na.action = na.omit)

# remove 'outliers'
# hits10    <- hits[c(1:6,8,10)] # determined by subsetting below
# weights10 <- weights[c(1:6,8,10)] # ditto
# 
# data$gm.score.top10inc <- dplyr::select(data,one_of(hits10)) %>% mapply("*",.,weights10) %>% rowSums()
# 
# outs <- data[complete.cases(data$gm),]
# cor.test(outs$gm,outs$gm.score.top10inc) # r = .588 (8 miRNA score)
# 
# outs <- outs[order(outs$gm.score.top10inc),]
# outs <- outs[1:73,] # remove sub with highest miRNA score 
# outs <- outs[order(outs$gm),]
# outs <- outs[2:73,] # remove sub with greatest reduction
# cor.test(outs$gm,outs$gm.score.top10inc) # r = .460 ('outliers' removed)

#### permutations
perm.data <- gm.comp
gm.vector <- gm.comp$gm
perm.df   <- matrix(nrow=2000,ncol=5)
for (i in 1:2000){
  perm.data$gm    <- sample(gm.vector)
  perm.calf       <- calf(data = perm.data, nMarkers = 10, targetVector="real")
  perm.calf.df    <- perm.calf$selection[!duplicated(perm.calf$selection$Marker),]
  hits.perm       <- perm.calf.df$Marker %>% as.character()
  weights.perm    <- perm.calf.df$Weight 
  perm.data$score <- dplyr::select(perm.data,one_of(hits.perm)) %>% 
    mapply("*",.,weights.perm) %>% rowSums()
  perm.corr       <- cor.test(perm.data$gm,perm.data$score)
  perm.df[i,1]    <- unname(perm.corr$statistic)
  perm.df[i,2]    <- perm.corr$p.value
  perm.df[i,3]    <- unname(perm.corr$estimate)
  perm.df[i,4]    <- unlist(perm.corr$conf.int)[1]
  perm.df[i,5]    <- unlist(perm.corr$conf.int)[2]
}

perm.df        <- as.data.frame(perm.df)
names(perm.df) <- c("tval","pval","corr","CI.low","CI.high")
perm.df        <- perm.df[order(perm.df$corr),]
perm.better    <- perm.df[perm.df$corr <= -0.695 | perm.df$corr >= 0.695,] # (1+77)/(1+2000) = .039

perm.df$corr.p <- abs(perm.df$corr) # make all positive

# plot permutations
ggplot(perm.df,aes(x=corr.p)) + 
  geom_histogram(binwidth=.05,aes(fill=..count..)) + 
  scale_fill_continuous(high="#0000FF",low="#00CCFF") +
  geom_segment(x=.69,y=0,xend=.69,yend=600,linetype="dashed") +
  geom_hline(yintercept=0) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        legend.position = "none") +
  ylab("Frequency") +
  xlab("Pearson Correlation (r)")

#### subsetting
sub.df   <- matrix(nrow=2000,ncol=5)
sub.hits <- NULL
for (i in 1:2000){
  sub.data       <- dplyr::sample_frac(gm.comp,0.8,replace=F)
  sub.calf       <- calf(data = sub.data, nMarkers = 10, targetVector="real")
  hits.sub       <- sub.calf$selection$Marker %>% as.character()
  sub.hits       <- c(sub.hits,hits.sub)
}

sub.hits.tbl  <- table(sub.hits) %>% as.data.frame() %>% dplyr::arrange(-Freq)

sub.hits.tbl$in_score <- ifelse(sub.hits.tbl$sub.hits %in% hits,"selected","not selected")
sub.hits.tbl$perc     <- sub.hits.tbl$Freq / 2000

write.table(sub.hits.tbl,"subsetting.txt",col.names=T,row.names=F,quote=F,sep="\t")

top20          <- sub.hits.tbl[1:20,]
top20$sub.hits <- droplevels(top20$sub.hits)
top20$sub.hits <- factor(top20$sub.hits,levels=top20$sub.hits[order(-top20$perc)],ordered=T)

# check weights for R2
sub.df   <- matrix(nrow=100,ncol=5)
sub.hits <- NULL
for (i in 1:100){
  sub.data       <- dplyr::sample_frac(gm.comp,0.8,replace=F)
  sub.calf       <- calf(data = sub.data, nMarkers = 10, targetVector="real")
  hits.sub       <- sub.calf$selection$Marker %>% as.character()
  save.weight    <- data.frame(hits.sub,sub.calf$selection$Weight)
  sub.hits       <- rbind(sub.hits,save.weight)
}

weights.tbl <- table(sub.hits$hits.sub,sub.hits$sub.calf.selection.Weight) %>% 
               as.data.frame 

uni.miRNAs <- unique(weights.tbl$Var1) # check if duplicates

names(weights.tbl) <- c("miRNA","weight","freq")
write.table(weights.tbl,"subsetting_weights.txt",col.names=T,row.names=F,quote=F,sep="\t")


# graph subsets
# colourCount <- length(unique(top20$sub.hits))
# getPalette  <- colorRampPalette(brewer.pal(9, "Spectral"))
ggplot(top20,aes(x=sub.hits,y=perc)) + 
  geom_bar(stat="identity",aes(fill=in_score)) + 
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("#0000FF","#00CCFF")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank()) +
  ylab("Frequency (%)") +
  xlab("miRNA")
 
#### bootstrapping
boot.df   <- matrix(nrow=2000,ncol=5)
boot.hits <- NULL
for (i in 1:2000){
  boot.data       <- dplyr::sample_frac(gm.comp,1,replace=T)
  boot.calf       <- calf(data = boot.data, nMarkers = 10, targetVector="real")
  hits.boot       <- boot.calf$selection$Marker %>% as.character()
  boot.hits       <- c(boot.hits,hits.boot)
}

boot.hits.tbl  <- table(boot.hits) %>% as.data.frame() %>% dplyr::arrange(-Freq)

boot.hits.tbl$in_score <- ifelse(boot.hits.tbl$boot.hits %in% hits,"selected","not selected")
boot.hits.tbl$perc     <- boot.hits.tbl$Freq / 2000

write.table(boot.hits.tbl,"bootstrapping.txt",col.names=T,row.names=F,quote=F,sep="\t")

top25boot           <- boot.hits.tbl[1:25,]
top25boot$boot.hits <- droplevels(top25boot$boot.hits)
top25boot$boot.hits <- factor(top25boot$boot.hits,levels=top25boot$boot.hits[order(-top25boot$perc)],ordered=T)

# graph subsets
ggplot(top25boot,aes(x=boot.hits,y=perc)) + 
  geom_bar(stat="identity",aes(fill=in_score)) + 
  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("#0000FF","#00CCFF")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        axis.text.x = element_text(angle=60, hjust=1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank()) +
  ylab("Frequency (%)") +
  xlab("miRNA")



#### FINAL CLASSIFIER (9 miRNA)

# miRNA classifier based on 
# past subsetting without rerunning
hits9    <- hits[c(1:8,10)]
weights9 <- weights[c(1:8,10)]

hits <- read.table("final_hits_miRNA_9.txt", stringsAsFactors = F)
hits <- hits[,1]
hits <- gsub("hsa-","", hits)
hits <- gsub("-",".", hits)
weights <- c(-1,-1,1,-1,1,-1,1,1,-1)
hits9 <- hits
weights9 <- weights

data$group2         <- data$group
levels(data$group2) <- c("not","P","not") 
levels(data$group2) <- c("Nonconverter or Control", "Converter")
levels(data$group)  <- c("Nonconverter","Converter", "Control")

# probability of overlap between classifiers
# sum(dhyper(t, a, n-a, b))
# t = overlap:b
# a = largest list n; b = smaller list n
# n = total list
# overlap = 1, conversion classifier has 5 miRNAs, GM has 9 miRNAs
# 136 total miRNAs tested (136-9 = 127)
sum(dhyper(1:5, 9, 127, 5)) # p = .294

# matrix(c(n - union(A,B), setdiff(A,B), setdiff(B,A), intersect(A,B)), nrow=2)
# A and B are list Ns; n = total list
fisher.test(matrix(c(123, 4, 8, 1), nrow=2)) # p = .294

# correlation between miRNA classifier and GM
data$gm.score.9 <- dplyr::select(data,one_of(hits9)) %>% mapply("*",.,weights9) %>% rowSums()
gm.real.9.m1    <- lm(gm ~ gm.score.9 + age + sex + group2, data = data, na.action = na.omit)

# case/control effect on classifier
# miR.cc.m1 <- lm(gm.score.9 ~ group + age + sex, data = data, na.action = na.omit)
# summary(miR.cc.m1) # converters vs non-converters
# anova(miR.cc.m1) # overall group effect

# interaction between dx and miRNA classifier on GM
data$group <- factor(data$group,levels(data$group)[c(2,1,3)])

inter  <- lm(gm ~ gm.score.9*group, data=data, na.action=na.omit) 
inter2 <- lm(gm ~ gm.score.9*group + sex + age, data = data, na.action = na.omit)
anova(inter2) # con vs. non-con vs. controls

inter3 <- lm(gm ~ gm.score.9*group2 + sex + age, data = data, na.action = na.omit)
anova(inter3) # converters vs. non-converters/controls

#### graph GM miRNA score against GM "#9999FF" "#99CCFF"
print(ggplot(data, aes(x=gm.score.9,y=gm,colour=group2)) + 
  geom_point(aes(colour=group2, shape=group), size=3) +
  scale_colour_manual(values=c("#0000FF","#00CCFF")) +
  scale_fill_manual(values=c("#0000FF","#00CCFF")) +
  geom_smooth(method=lm,aes(fill=group2), alpha=.2) + 
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(.8,.3),
        legend.text = element_text(size=13, color="black"),
        legend.text.align = 1,
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  ylab("Gray Matter Change") +
  xlab("miRNA Score"))
pdf("miRNA-GM-dxslopes3.pdf", width =6, height =4)
dev.copy(pdf, "miRNA-GM-dxslopes3.pdf")
dev.off()

#### other ROIs 
# primary ROI is label 3 - right sup frontal
# label 1 is medial OFC; label 2 is middle frontal
ROIs            <- read.table("napls_ROIs.txt",header=T,stringsAsFactors = F)
ROIs$SiteSubjID <- paste(ROIs$SiteNumber,ROIs$SubjectNumber,sep="") %>% as.integer()
data$label1     <- ROIs[match(data$SiteSubjID,ROIs$SiteSubjID),4]
data$label2     <- ROIs[match(data$SiteSubjID,ROIs$SiteSubjID),5]
data$label4     <- ROIs[match(data$SiteSubjID,ROIs$SiteSubjID),8]
data$label5     <- ROIs[match(data$SiteSubjID,ROIs$SiteSubjID),9]

label1 <- lm(label1 ~ gm.score.9 + age + sex + group2, data = data, na.action = na.omit)
label2 <- lm(label2 ~ gm.score.9 + age + sex + group2, data = data, na.action = na.omit)

#### cytokines x miRNA classifier (9)
# build cytokine vars (sum of pro; anti)
# based on the neg values, i infer these have already been z-scored
# BP paper indicated z-scores were based on mean / SD of control group
cyto            <- read.table("cytokines.txt",header=T,stringsAsFactors = F)
cyto$SiteSubjID <- paste(cyto$sitenumber,cyto$subjectnumber,sep="") %>% as.integer()

# pro-inflammatory = TNF-alpha, IL-2, interferon-gamma
# anti-inflammatory = IL-10, granulocyte-macrophage colony-stimulating factor,
# IL-1 receptor antagonist, chemokine (C-C motif) ligand; 
# (Ty said he also included 'CC12' but not in the sheet or paper)
pro  <- c("TNF","IL2","IFNG")
anti <- c("IL10","CSF2", "IL1RN", "TGFB3")

cyto$pro.score  <- dplyr::select(cyto,one_of(pro)) %>% rowSums()
cyto$anti.score <- dplyr::select(cyto,one_of(anti)) %>% rowSums()
data$pro.score  <- cyto[match(data$SiteSubjID,cyto$SiteSubjID),188]
data$anti.score <- cyto[match(data$SiteSubjID,cyto$SiteSubjID),189]

# only subjects included
data.inc <- data[complete.cases(data$gm),]
data.inc$pro.score.win  <- psych::winsor(data.inc$pro.score, trim=0.02)
# data.inc$anti.score.win <- psych::winsor(data.inc$anti.score, trim=0.02)

cyto.pro.m1  <- lm(pro.score ~ gm.score.9 + age + sex + group, data = data.inc, na.action = na.omit)
cyto.anti.m1 <- lm(anti.score ~ gm.score.9 + age + sex + group, data = data.inc, na.action = na.omit)
cyto.pro.m2  <- lm(pro.score ~ gm.score.9 + age + sex + group2, data = data.inc, na.action = na.omit)
cyto.anti.m2 <- lm(anti.score ~ gm.score.9 + age + sex + group2, data = data.inc, na.action = na.omit)

# robust regression bc outlier
cyto.pro.m3 <- rlm(pro.score ~ gm.score.9 + age + sex + group2, method="M", data = data.inc, na.action = na.omit)
f.robftest(cyto.pro.m3, var="gm.score.9") # Wald test

# r2 <- function(x){
#   SSe      <- sum((x$resid)^2);  
#   observed <- x$resid + x$fitted;  
#   SSt      <- sum((observed-mean(observed))^2);  
#   value    <- 1-SSe/SSt;   
#   return(value);  
# } 
# r2(cyto.pro.m3)

cyto.pro.m4 <- lm(pro.score.win ~ gm.score.9 + age + sex + group2, data = data.inc, na.action = na.omit)
# cyto.pro.m5 <- lmrob(pro.score ~ gm.score.9 + age + sex + group2, method="M", init="S", data = data.inc, na.action = na.omit)
# cyto.pro.m6 <- lmRob(pro.score ~ gm.score.9 + age + sex + group2, data = data.inc, na.action = na.omit)


#### graph GM miRNA score against pro-inflammatory cytokines
levels(data.inc$group) <- c("Nonconverter", "Converter", "Control")
ggplot(data.inc, aes(x=gm.score.9,y=pro.score)) + 
  geom_point(aes(colour=group, shape=group), size=4) +
  scale_colour_manual(values=c("#0000FF","#00CCFF","#9999FF")) +
  geom_smooth(method=lm, fill="#99CCFF") + 
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  ylab("Pro-Inflammatory Cytokines") +
  xlab("miRNA Score")


### graph gene target pathways
# pathway file is sort of nonsense but the order is correct
pathway <- read.table("ipa_pathways.txt",header=T,sep="\t")
pathway <- cbind(rownames(pathway),pathway)
names(pathway) <- c("pathway","BH.pval","ratio","zscore","molecules")


#### CALF excluding 10 selected
# gm.void      <- gm.comp[,!(colnames(gm.comp) %in% hits)]
# gm.void.calf <- calf(data = gm.void, nMarkers = 10, targetVector="real")
# 
# hits.void    <- gm.void.calf$selection$Marker %>% as.character()
# weights.void <- gm.void.calf$selection$Weight 
# data$gm.void <- dplyr::select(data,one_of(hits.void)) %>% mapply("*",.,weights.void) %>% rowSums()
# cor.test(data$gm.void,data$gm) # r = .626
# 
# # permutations
# perm.void <- gm.void
# gm.vector <- gm.void$gm
# void.df   <- matrix(nrow=2000,ncol=5)
# for (i in 1:2000){
#   perm.void$gm    <- sample(gm.vector)
#   perm.calf       <- calf(data = perm.void, nMarkers = 10, targetVector="real")
#   perm.calf.df    <- perm.calf$selection[!duplicated(perm.calf$selection$Marker),]
#   hits.perm       <- perm.calf.df$Marker %>% as.character()
#   weights.perm    <- perm.calf.df$Weight 
#   perm.void$score <- dplyr::select(perm.void,one_of(hits.perm)) %>% 
#     mapply("*",.,weights.perm) %>% rowSums()
#   perm.corr       <- cor.test(perm.void$gm,perm.void$score)
#   void.df[i,1]    <- unname(perm.corr$statistic)
#   void.df[i,2]    <- perm.corr$p.value
#   void.df[i,3]    <- unname(perm.corr$estimate)
#   void.df[i,4]    <- unlist(perm.corr$conf.int)[1]
#   void.df[i,5]    <- unlist(perm.corr$conf.int)[2]
# }
# 
# void.df        <- as.data.frame(void.df)
# names(void.df) <- c("tval","pval","corr","CI.low","CI.high")
# void.df        <- void.df[order(void.df$corr),]
# void.better    <- void.df[void.df$corr <= -0.626 | void.df$corr >= 0.626,] # (1+435)/(1+2000) = .218

# test <- read.table("miRNA_GM_agesex.txt",stringsAsFactors = F,header=T)
# gm.test <- test[,c(4,7:142)]
# gm.calf.test <- calf(data = gm.test, nMarkers = 10, targetVector="real")
# clark.df <- read.table("GM_data.csv",stringsAsFactors = F,header=T,sep=",")
# 
# clark <- read.table("Clark_solution.txt",stringsAsFactors = F,header=F)
# clark.hits <- clark$V1
# clark.hits <- gsub("-",".",clark.hits,fixed=T)
# clark.weights <- clark$V2
# data$clark.score <- dplyr::select(data,one_of(clark.hits)) %>% mapply("*",.,clark.weights) %>% rowSums()
# cor.test(data$clark.score,data$gm)


###### functional annotations
# miranda <- read.csv("human_predictions_S_C_aug2010.txt",
#                     header=T, stringsAsFactors = F, sep="\t")
# 
# # lists of miRNAs to classify
# # match miranda formatting
# hits.m <- gsub(".","-",hits,fixed=T)
# hits.m <- paste0("hsa-",hits.m)
# hits.m <- lapply(strsplit(hits.m,"-"), function(x) paste(x[1:3],collapse="-")) %>% unlist
# 
# # manual edits bc I need fuzzy searches and don't know how
# hits.m[1]  <- "hsa-miR-103"
# hits.m[2]  <- "hsa-miR-140-5p" # ours is 3p
# hits.m[3]  <- "hsa-miR-142-3p" # ours is 5p
# hits.m[8]  <- "hsa-miR-339-5p" # ours is 3p
# hits.m[10] <- "hsa-miR-193a-3p" # ours is 5p
#   
# # remove miRNAs that are not predictors
# # remove gene targets with mirsvr scores < -.1
# all.m <- miranda[miranda$mirna_name %in% hits.m,] # uses 9 of 10 (35,922 rows)
#          # no miR-501 
# all.m <- all.m[all.m$mirsvr_score < -.1,] # 35,911 rows
# 
# # use table to determine how many times each gene is listed
# # discard either non-duplicates or all not listed X times
# genefreq <- table(all.m$gene_symbol) %>% as.data.frame()
# genes    <- genefreq[genefreq$Freq > 9,] # 348 genes overlapping all miRNAs *****
# genes2   <- genefreq[genefreq$Freq > 8,] # 455 genes overlapping all miRNAs
# 
# # use intersection list for functional annotation
# write.table(genes[,1],"genelists/gm.calf.miranda.txt",col.names=F,row.names=F,quote=F,sep="\t")
# write.table(genes2[,1],"genelists/gm.calf.miranda2.txt",col.names=F,row.names=F,quote=F,sep="\t")



