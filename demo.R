# demo stats for miRNA paper

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/RNAseq/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2","Hmisc")
invisible(lapply(libs, require, character.only = TRUE))

#### load in data
demo <- read.table("../napls_pheno-demo.txt",header=T,stringsAsFactors = F)
data <- read.table("miRNA_GM_agesex.txt",header=T)

demo$SiteSubjID <- paste(demo$SiteNumber,demo$SubjectNumber,sep="") %>% as.integer()
demo$group      <- data[match(demo$SiteSubjID,data$SiteSubjID),2]
df              <- demo[demo$SiteSubjID %in% data$SiteSubjID,]

#### demo stats
# age
describeBy(df$demo_age_ym,group=df$group)
# group.CI(demo_age_ym~group, data=df, ci=0.95)

m1 <- lm(demo_age_ym ~ group, data = df)
summary(m1)

# sex
table(df$demo_sex,df$group)
chisq.test(df$demo_sex, df$group)

# race
# european = 8; black = 5; mixed = 10
# asian = 2 + 3 + 4 + 7
# american = 1 + 6 + 9
df$race_group         <- df$demo_racial %>% as.factor
levels(df$race_group) <- list(american=c("1","6","9"), asian=c("2","3","4","7"),
                              european=8, black=5, mixed=10)

race     <- table(df$race_group,df$group) %>% as.data.frame
race$tot <- ifelse(race$Var2=="NP",34,ifelse(race$Var2=="P",13,27))
race$per <- race$Freq / race$tot

chisq.test(df$race_group,df$group)

# SES
# treated as continuous (??)
describeBy(df$demo_income,group=df$group)
m2 <- lm(demo_income ~ group, data = df)
summary(m2)

# education
describeBy(df$demo_education_years,group=df$group)
m3 <- lm(demo_education_years ~ group, data = df)
summary(m3)

# SOPS
sops            <- read.table("../SOPS.txt",header=T,stringsAsFactors = F)
sops$SiteSubjID <- paste(sops$SiteNumber,sops$SubjectNumber,sep="") %>% as.integer()
sops$group      <- data[match(sops$SiteSubjID,data$SiteSubjID),2]
sops.use        <- sops[sops$SiteSubjID %in% data$SiteSubjID,]

sops.vars           <- c("cSOPSPositive","cSOPSNegative","cSOPSDisorganization","cSOPSGeneral")
sops.use$cSOPSTotal <- rowSums(sops.use[,sops.vars])
describeBy(sops.use[,c(sops.vars,"cSOPSTotal")],sops.use$group)

s1 <- lm(cSOPSPositive ~ group, data = sops.use[sops.use$group!="UC",])
s2 <- lm(cSOPSNegative ~ group, data = sops.use[sops.use$group!="UC",])
s3 <- lm(cSOPSDisorganization ~ group, data = sops.use[sops.use$group!="UC",])
s4 <- lm(cSOPSGeneral ~ group, data = sops.use[sops.use$group!="UC",])
s5 <- lm(cSOPSTotal ~ group, data = sops.use[sops.use$group!="UC",])
s6 <- lm(cSOPSTotal ~ group, data = sops.use)

# medications
# multiple entries for each subject (1 drug per entry)
# use only baseline meds
# same med is entered multiple times for same subject + time point
# these entries have different "instance" numbers (???) -- ignored
meds <- read.csv("../meds.txt",sep="\t",header=T,stringsAsFactors = F)
meds$SiteSubjID <- paste(meds$SiteNumber,meds$SubjectNumber,sep="") %>% as.integer()
meds.use        <- meds[meds$SiteSubjID %in% data$SiteSubjID & meds$VisitLabel=="BL",]

meds.long       <- table(meds.use$SiteSubjID,meds.use$medlist_med_name) %>% as.data.frame()
meds.wide       <- reshape(meds.long[,c(1:3)],idvar="Var1",direction="wide",timevar="Var2")
meds.wide$group <- data[match(meds.wide$Var1,data$SiteSubjID),2]

meds.wide[meds.wide > 1] <- 1 # change meaningless high freq to 1 (result = only 1s & 0s)
names(meds.wide) <- gsub("Freq.","",names(meds.wide),fixed=T)

nomeds <- c(" No meds", "No Meds", "None", "No meds")
aps    <- c("Abilify","Haldol","Haloperidol","Risperdal","Risperidone","Seroquel","Zyprexa")
ads    <- c("Amitriptyline","Celexa","Effexor","Fluoxetine","Lexapro","Luvox",
            "Paxil","Prozac","Strattera","Trazadone","Trazodone","Wellbutrin","Zoloft")
stims  <- c("adderall","Adderall","Concerta","Dexadrine","dexedrine","Dexedrine",
            "Intuniv","Metadate","Ritalin","Vyvanse")
moodst <- c("Depakote")
benzo  <- c("Ativan","Clonazepam","Phenobarbital","Xanax")
sedat  <- c("Ambien","Clonidine","Oxycodone","Percocet","Tylenol Codeine")
NSAID  <- c("Aleve","Ibuprophen","Toradol")
otc    <- c("Midol","Motrin","Vitamin","Vitamins")
hist   <- c("Benadryl","Benedryl","Claritan","Claritin","Flonase","Loratadine","Meclizine","Zyrtec D")
asthma <- c("Advair","albuterol","Albuterol","Inhaler","Proair")
other  <- c("Accutane","acne meds","Amoxicillin","Antibiotic","Antibiotics",
            "Birth control","Birth Control","Clindomycin","Cogentin","Loestrin24Fe",
            "Penicillan","Penicillin","Trileptal","Z-Pac/Azithromyein")

meds.wide$none   <- rowSums(meds.wide[,nomeds])
meds.wide$aps    <- rowSums(meds.wide[,aps])
meds.wide$ads    <- rowSums(meds.wide[,ads])
meds.wide$stims  <- rowSums(meds.wide[,stims])
meds.wide$moodst <- meds.wide$Depakote
meds.wide$benzo  <- rowSums(meds.wide[,benzo])
meds.wide$sedat  <- rowSums(meds.wide[,sedat])
meds.wide$NSAID  <- rowSums(meds.wide[,NSAID])
meds.wide$otc    <- rowSums(meds.wide[,otc])
meds.wide$hist   <- rowSums(meds.wide[,hist])
meds.wide$asthma <- rowSums(meds.wide[,asthma])
meds.wide$other  <- rowSums(meds.wide[,other])

meds.wide[meds.wide > 1] <- 1 # again bc some ind on many of same category

meds.NP      <- meds.wide[meds.wide$group=="NP",85:96] 
meds.NP.sums <- colSums(meds.NP) %>% as.data.frame

meds.P      <- meds.wide[meds.wide$group=="P",85:96] 
meds.P.sums <- colSums(meds.P) %>% as.data.frame

meds.UC      <- meds.wide[meds.wide$group=="UC",85:96] 
meds.UC.sums <- colSums(meds.UC) %>% as.data.frame

meds.all        <- cbind(meds.NP.sums,meds.P.sums,meds.UC.sums)
meds.all$class  <- rownames(meds.all)
names(meds.all)[1:3] <- c("NP","P","UC")

meds.all$NPper <- meds.all$NP / 34
meds.all$Pper  <- meds.all$P / 13
meds.all$UCper <- meds.all$UC / 27

write.table(meds.all,"../meds_summary.txt",col.names=T,row.names=F,quote=F,sep="\t")

meds.CHR <- meds.wide[meds.wide$group!="UC",]
chisq.test(meds.CHR$aps,meds.CHR$group)
chisq.test(meds.CHR$ads,meds.CHR$group)
chisq.test(meds.CHR$stims,meds.CHR$group)
chisq.test(meds.CHR$moodst,meds.CHR$group)
chisq.test(meds.CHR$benzo,meds.CHR$group)
chisq.test(meds.CHR$NSAID,meds.CHR$group)
chisq.test(meds.CHR$hist,meds.CHR$group)
chisq.test(meds.CHR$sedat,meds.CHR$group)
chisq.test(meds.CHR$other,meds.CHR$group)
chisq.test(meds.CHR$none,meds.CHR$group)

# group differences in excluded subjects
# original sample: 764 CHR, 280 controls
# miRNA sample: 71 CHR, 30 controls
# excluded subs: 23 CHR, 3 controls
chisq  <- matrix(data=c(764,280,23,3),nrow=2) %>% data.frame
chisq.test(chisq)
chisq2 <- matrix(data=c(71,30,23,3),nrow=2) %>% data.frame
chisq.test(chisq2)
