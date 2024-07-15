#dx download PCH_allcovariates_sex_participant.tsv
library(dplyr)
df = read.delim("PCH_allcovariates_sex_participant.tsv")

str(df)
#begin with analysis of M/F discrepency in autoimmune diseases and cancers
lupus = df[grep("M32",df$X41270.0.0),]
df$lupus <- ifelse(grepl("M32", df$X41270.0.0), 1, 0)

Model = glm(lupus~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(lupus$X31.0.0)[1]
Males = table(lupus$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
lupusmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "lupus")

hashimoto = df[grep("E06",df$X41270.0.0),]
df$hashimoto <- ifelse(grepl("E06", df$X41270.0.0), 1, 0)

Model = glm(hashimoto~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(hashimoto$X31.0.0)[1]
Males = table(hashimoto$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
hashimotomerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "hashimoto")

graves = df[grep("E05",df$X41270.0.0),]
df$graves <- ifelse(grepl("E05", df$X41270.0.0), 1, 0)

Model = glm(graves~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(graves$X31.0.0)[1]
Males = table(graves$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
gravesmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "graves")

MS = df[grep("G35",df$X41270.0.0),]
df$MS <- ifelse(grepl("G35", df$X41270.0.0), 1, 0)

Model = glm(MS~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(MS$X31.0.0)[1]
Males = table(MS$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MSmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "MS")

RA = df[grep("M05",df$X41270.0.0),]
df$RA <- ifelse(grepl("M05", df$X41270.0.0), 1, 0)

Model = glm(RA~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(RA$X31.0.0)[1]
Males = table(RA$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
RAmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "RA")

ThyroidC = df[grep("C73",df$X41270.0.0),]
df$ThyroidC <- ifelse(grepl("C73", df$X41270.0.0), 1, 0)

Model = glm(ThyroidC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(ThyroidC$X31.0.0)[1]
Males = table(ThyroidC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
ThyroidCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "ThyroidC")

PBC = df[grep("K743",df$X41270.0.0),]
df$PBC <- ifelse(grepl("K743", df$X41270.0.0), 1, 0)

Model = glm(PBC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(PBC$X31.0.0)[1]
Males = table(PBC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PBCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "PBC")

Alzheimers = df[grep("G30",df$X41270.0.0),]
df$Alzheimers <- ifelse(grepl("G30", df$X41270.0.0), 1, 0)

Model = glm(Alzheimers~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(Alzheimers$X31.0.0)[1]
Males = table(Alzheimers$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Alzheimersmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "Alzheimers")

Celiac = df[grep("K900",df$X41270.0.0),]
df$Celiac <- ifelse(grepl("K900", df$X41270.0.0), 1, 0)

Model = glm(Celiac~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(Celiac$X31.0.0)[1]
Males = table(Celiac$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Celiacmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "Celiac")

PA = df[grep("L405",df$X41270.0.0),]
df$PA <- ifelse(grepl("L405", df$X41270.0.0), 1, 0)

Model = glm(PA~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(PA$X31.0.0)[1]
Males = table(PA$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PAmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "PA")

MG = df[grep("G700",df$X41270.0.0),]
df$MG <- ifelse(grepl("G700", df$X41270.0.0), 1, 0)

Model = glm(MG~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(MG$X31.0.0)[1]
Males = table(MG$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MGmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "MG")

CD = df[grep("K50",df$X41270.0.0),]
df$CD <- ifelse(grepl("K50", df$X41270.0.0), 1, 0)

Model = glm(CD~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(CD$X31.0.0)[1]
Males = table(CD$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
CDmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "CD")

psoriasis = df[grep("L40",df$X41270.0.0),]
df$psoriasis <- ifelse(grepl("L40", df$X41270.0.0), 1, 0)

Model = glm(psoriasis~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(psoriasis$X31.0.0)[1]
Males = table(psoriasis$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
psoriasismerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "psoriasis")

UC = df[grep("K51",df$X41270.0.0),]
df$UC <- ifelse(grepl("K51", df$X41270.0.0), 1, 0)

Model = glm(UC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(UC$X31.0.0)[1]
Males = table(UC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
UCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "UC")

T1D = df[grep("E10",df$X41270.0.0),]
df$T1D <- ifelse(grepl("E10", df$X41270.0.0), 1, 0)

Model = glm(T1D~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(T1D$X31.0.0)[1]
Males = table(T1D$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
T1Dmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "T1D")

ColorectalC = df[grep("C18",df$X41270.0.0),]
df$ColorectalC <- ifelse(grepl("C18", df$X41270.0.0), 1, 0)

Model = glm(ColorectalC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(ColorectalC$X31.0.0)[1]
Males = table(ColorectalC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
ColorectalCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "ColorectalC")

MelanomaC = df[grep("C43",df$X41270.0.0),]
df$MelanomaC <- ifelse(grepl("C43", df$X41270.0.0), 1, 0)

Model = glm(MelanomaC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(MelanomaC$X31.0.0)[1]
Males = table(MelanomaC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MelanomaCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "MelanomaC")

MyelomaC = df[grep("C90",df$X41270.0.0),]
df$MyelomaC <- ifelse(grepl("C90", df$X41270.0.0), 1, 0)

Model = glm(MyelomaC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(MyelomaC$X31.0.0)[1]
Males = table(MyelomaC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MyelomaCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "MyelomaC")

LungC = df[grep("C34",df$X41270.0.0),]
df$LungC <- ifelse(grepl("C34", df$X41270.0.0), 1, 0)

Model = glm(LungC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(LungC$X31.0.0)[1]
Males = table(LungC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
LungCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "LungC")

Ankspond = df[grep("M45",df$X41270.0.0),]
df$Ankspond <- ifelse(grepl("M45", df$X41270.0.0), 1, 0)

Model = glm(Ankspond~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(Ankspond$X31.0.0)[1]
Males = table(Ankspond$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Ankspondmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "Ankspond")

KidneyC = df[grep("C64",df$X41270.0.0),]
df$KidneyC <- ifelse(grepl("C64", df$X41270.0.0), 1, 0)

Model = glm(KidneyC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(KidneyC$X31.0.0)[1]
Males = table(KidneyC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
KidneyCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "KidneyC")

LiverC = df[grep("C22",df$X41270.0.0),]
df$LiverC <- ifelse(grepl("C22", df$X41270.0.0), 1, 0)

Model = glm(LiverC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(LiverC$X31.0.0)[1]
Males = table(LiverC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
LiverCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "LiverC")

BladderC = df[grep("C67",df$X41270.0.0),]
df$BladderC <- ifelse(grepl("C67", df$X41270.0.0), 1, 0)

Model = glm(BladderC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(BladderC$X31.0.0)[1]
Males = table(BladderC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
BladderCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "BladderC")

EsophogaelC = df[grep("C15",df$X41270.0.0),]
df$EsophogaelC <- ifelse(grepl("C15", df$X41270.0.0), 1, 0)

Model = glm(EsophogaelC~X31.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = df, family = binomial)
Females = table(EsophogaelC$X31.0.0)[1]
Males = table(EsophogaelC$X31.0.0)[2]
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
EsophogaelCmerge = data.frame(Females, Males, Betas, StErr, Pvals, row.names = "EsophogaelC")

Statistics = rbind(lupusmerge, hashimotomerge,gravesmerge,MSmerge, RAmerge, ThyroidCmerge,PBCmerge,
                   Alzheimersmerge,Celiacmerge,PAmerge,MGmerge,CDmerge,psoriasismerge,UCmerge,T1Dmerge,
                   ColorectalCmerge,MelanomaCmerge,MyelomaCmerge,LungCmerge,Ankspondmerge,KidneyCmerge,
                   LiverCmerge,BladderCmerge,EsophogaelCmerge)

#PCH hypothesis for women
Women = subset(df, X31.0.0 == 0)

Women$TotalBirths = Women$X2734.0.0 + Women$X3829.0.0
Women$TotalPreg = Women$X2734.0.0 + Women$X3829.0.0 + 0.5*Women$X3839.0.0 + 0.5*Women$X3849.0.0

Women$Autoimmune <- ifelse(grepl("M32|E06|E05|G35|M05|K743|G30|K900|L405|G700|K50|L40|K51|E10|M45", Women$X41270.0.0), 1, 0)
Women$SigAutoimmune <- ifelse(grepl("E05|G35|M32|E06|M05|K900|K743", Women$X41270.0.0), 1, 0)
  
#get rid of women who prefered not to answer (-3) or didn't know (-1)
Women1<-subset(Women, X2734.0.0!=-3)
Women1<-subset(Women1, X3829.0.0!=-3)
Women1<-subset(Women1, X3829.0.0!=-1)

Model = glm(Autoimmune ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
coef(summary(Model))

Model = glm(SigAutoimmune ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
coef(summary(Model))

Model = glm(graves ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
gravesmerge = data.frame(Betas, StErr, Pvals, row.names = "graves")

Model = glm(lupus ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
lupusmerge = data.frame(Betas, StErr, Pvals, row.names = "lupus")

Model = glm(hashimoto ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
hashimotomerge = data.frame(Betas, StErr, Pvals, row.names = "hashimoto")

Model = glm(graves ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MSmerge = data.frame(Betas, StErr, Pvals, row.names = "MS")

Model = glm(RA ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
RAmerge = data.frame(Betas, StErr, Pvals, row.names = "RA")

Model = glm(PBC ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PBCmerge = data.frame(Betas, StErr, Pvals, row.names = "PBC")

Model = glm(Alzheimers ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Alzheimersmerge = data.frame(Betas, StErr, Pvals, row.names = "Alzheimers")

Model = glm(Celiac ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Celiacmerge = data.frame(Betas, StErr, Pvals, row.names = "Celiac")

Model = glm(PA ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PAmerge = data.frame(Betas, StErr, Pvals, row.names = "PA")

Model = glm(MG ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MGmerge = data.frame(Betas, StErr, Pvals, row.names = "MG")

Model = glm(CD ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
CDmerge = data.frame(Betas, StErr, Pvals, row.names = "CD")

Model = glm(psoriasis ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
psoriasismerge = data.frame(Betas, StErr, Pvals, row.names = "psoriasis")

Model = glm(UC ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
UCmerge = data.frame(Betas, StErr, Pvals, row.names = "UC")

Model = glm(T1D ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
T1Dmerge = data.frame(Betas, StErr, Pvals, row.names = "T1D")

Model = glm(Ankspond ~ TotalBirths + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women1, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Ankspondmerge = data.frame(Betas, StErr, Pvals, row.names = "Ankspond")

StatisticsBirths = rbind(lupusmerge, hashimotomerge,gravesmerge,MSmerge, RAmerge, PBCmerge,
                     Alzheimersmerge,Celiacmerge,PAmerge,MGmerge,CDmerge,psoriasismerge,UCmerge,T1Dmerge,
                     Ankspondmerge)

#get rid of women who prefered not to answer (-3) or didn't know (-1) about spontaneous abortions
Women2<-subset(Women1, X3839.0.0!=-3)
Women2<-subset(Women2, X3839.0.0!=-1)
Women2<-subset(Women2, X3849.0.0!=-1)
Women2<-subset(Women2, X3849.0.0!=-3)

Model = glm(Autoimmune ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
coef(summary(Model))

Model = glm(SigAutoimmune ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
coef(summary(Model))

Model = glm(graves ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
gravesmerge = data.frame(Betas, StErr, Pvals, row.names = "graves")

Model = glm(lupus ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
lupusmerge = data.frame(Betas, StErr, Pvals, row.names = "lupus")

Model = glm(hashimoto ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
hashimotomerge = data.frame(Betas, StErr, Pvals, row.names = "hashimoto")

Model = glm(graves ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MSmerge = data.frame(Betas, StErr, Pvals, row.names = "MS")

Model = glm(RA ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
RAmerge = data.frame(Betas, StErr, Pvals, row.names = "RA")

Model = glm(PBC ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PBCmerge = data.frame(Betas, StErr, Pvals, row.names = "PBC")

Model = glm(Alzheimers ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Alzheimersmerge = data.frame(Betas, StErr, Pvals, row.names = "Alzheimers")

Model = glm(Celiac ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Celiacmerge = data.frame(Betas, StErr, Pvals, row.names = "Celiac")

Model = glm(PA ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PAmerge = data.frame(Betas, StErr, Pvals, row.names = "PA")

Model = glm(MG ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MGmerge = data.frame(Betas, StErr, Pvals, row.names = "MG")

Model = glm(CD ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
CDmerge = data.frame(Betas, StErr, Pvals, row.names = "CD")

Model = glm(psoriasis ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
psoriasismerge = data.frame(Betas, StErr, Pvals, row.names = "psoriasis")

Model = glm(UC ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
UCmerge = data.frame(Betas, StErr, Pvals, row.names = "UC")

Model = glm(T1D ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
T1Dmerge = data.frame(Betas, StErr, Pvals, row.names = "T1D")

Model = glm(Ankspond ~ TotalPreg + X34.0.0 + X21000.0.0 + X22189.0.0, data = Women2, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Ankspondmerge = data.frame(Betas, StErr, Pvals, row.names = "Ankspond")

StatisticsPregnancies = rbind(lupusmerge, hashimotomerge,gravesmerge,MSmerge, RAmerge, PBCmerge,
                         Alzheimersmerge,Celiacmerge,PAmerge,MGmerge,CDmerge,psoriasismerge,UCmerge,T1Dmerge,
                         Ankspondmerge)

Menarche<-subset(Women, X2714.0.0!=-1)
Menarche<-subset(Menarche, X2714.0.0!=-3)

Model = glm(Autoimmune ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
coef(summary(Model))

Model = glm(SigAutoimmune ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
coef(summary(Model))

Model = glm(graves ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
gravesmerge = data.frame(Betas, StErr, Pvals, row.names = "graves")

Model = glm(lupus ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
lupusmerge = data.frame(Betas, StErr, Pvals, row.names = "lupus")

Model = glm(hashimoto ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
hashimotomerge = data.frame(Betas, StErr, Pvals, row.names = "hashimoto")

Model = glm(graves ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MSmerge = data.frame(Betas, StErr, Pvals, row.names = "MS")

Model = glm(RA ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
RAmerge = data.frame(Betas, StErr, Pvals, row.names = "RA")

Model = glm(PBC ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PBCmerge = data.frame(Betas, StErr, Pvals, row.names = "PBC")

Model = glm(Alzheimers ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Alzheimersmerge = data.frame(Betas, StErr, Pvals, row.names = "Alzheimers")

Model = glm(Celiac ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Celiacmerge = data.frame(Betas, StErr, Pvals, row.names = "Celiac")

Model = glm(PA ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PAmerge = data.frame(Betas, StErr, Pvals, row.names = "PA")

Model = glm(MG ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MGmerge = data.frame(Betas, StErr, Pvals, row.names = "MG")

Model = glm(CD ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
CDmerge = data.frame(Betas, StErr, Pvals, row.names = "CD")

Model = glm(psoriasis ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
psoriasismerge = data.frame(Betas, StErr, Pvals, row.names = "psoriasis")

Model = glm(UC ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
UCmerge = data.frame(Betas, StErr, Pvals, row.names = "UC")

Model = glm(T1D ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
T1Dmerge = data.frame(Betas, StErr, Pvals, row.names = "T1D")

Model = glm(Ankspond ~ X2714.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menarche, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Ankspondmerge = data.frame(Betas, StErr, Pvals, row.names = "Ankspond")

StatisticsMenarche = rbind(lupusmerge, hashimotomerge,gravesmerge,MSmerge, RAmerge, PBCmerge,
                         Alzheimersmerge,Celiacmerge,PAmerge,MGmerge,CDmerge,psoriasismerge,UCmerge,T1Dmerge,
                         Ankspondmerge)

Menopause<-subset(Women, X3581.0.0!=-1)
Menopause<-subset(Menopause, X3581.0.0!=-3)

Model = glm(Autoimmune ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
coef(summary(Model))

Model = glm(SigAutoimmune ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
coef(summary(Model))

Model = glm(graves ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
gravesmerge = data.frame(Betas, StErr, Pvals, row.names = "graves")

Model = glm(lupus ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
lupusmerge = data.frame(Betas, StErr, Pvals, row.names = "lupus")

Model = glm(hashimoto ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
hashimotomerge = data.frame(Betas, StErr, Pvals, row.names = "hashimoto")

Model = glm(graves ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MSmerge = data.frame(Betas, StErr, Pvals, row.names = "MS")

Model = glm(RA ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
RAmerge = data.frame(Betas, StErr, Pvals, row.names = "RA")

Model = glm(PBC ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PBCmerge = data.frame(Betas, StErr, Pvals, row.names = "PBC")

Model = glm(Alzheimers ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Alzheimersmerge = data.frame(Betas, StErr, Pvals, row.names = "Alzheimers")

Model = glm(Celiac ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Celiacmerge = data.frame(Betas, StErr, Pvals, row.names = "Celiac")

Model = glm(PA ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PAmerge = data.frame(Betas, StErr, Pvals, row.names = "PA")

Model = glm(MG ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MGmerge = data.frame(Betas, StErr, Pvals, row.names = "MG")

Model = glm(CD ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
CDmerge = data.frame(Betas, StErr, Pvals, row.names = "CD")

Model = glm(psoriasis ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
psoriasismerge = data.frame(Betas, StErr, Pvals, row.names = "psoriasis")

Model = glm(UC ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
UCmerge = data.frame(Betas, StErr, Pvals, row.names = "UC")

Model = glm(T1D ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
T1Dmerge = data.frame(Betas, StErr, Pvals, row.names = "T1D")

Model = glm(Ankspond ~ X3581.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = Menopause, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Ankspondmerge = data.frame(Betas, StErr, Pvals, row.names = "Ankspond")

StatisticsMenopause = rbind(lupusmerge, hashimotomerge,gravesmerge,MSmerge, RAmerge, PBCmerge,
                   Alzheimersmerge,Celiacmerge,PAmerge,MGmerge,CDmerge,psoriasismerge,UCmerge,T1Dmerge,
                   Ankspondmerge)


BirthControlStart<-subset(Women, X2794.0.0!=-1)
BirthControlStart<-subset(BirthControlStart, X2794.0.0!=-3)

Model = glm(Autoimmune ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
coef(summary(Model))

Model = glm(SigAutoimmune ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
coef(summary(Model))

Model = glm(graves ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
gravesmerge = data.frame(Betas, StErr, Pvals, row.names = "graves")

Model = glm(lupus ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
lupusmerge = data.frame(Betas, StErr, Pvals, row.names = "lupus")

Model = glm(hashimoto ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
hashimotomerge = data.frame(Betas, StErr, Pvals, row.names = "hashimoto")

Model = glm(graves ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MSmerge = data.frame(Betas, StErr, Pvals, row.names = "MS")

Model = glm(RA ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
RAmerge = data.frame(Betas, StErr, Pvals, row.names = "RA")

Model = glm(PBC ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PBCmerge = data.frame(Betas, StErr, Pvals, row.names = "PBC")

Model = glm(Alzheimers ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Alzheimersmerge = data.frame(Betas, StErr, Pvals, row.names = "Alzheimers")

Model = glm(Celiac ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Celiacmerge = data.frame(Betas, StErr, Pvals, row.names = "Celiac")

Model = glm(PA ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PAmerge = data.frame(Betas, StErr, Pvals, row.names = "PA")

Model = glm(MG ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MGmerge = data.frame(Betas, StErr, Pvals, row.names = "MG")

Model = glm(CD ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
CDmerge = data.frame(Betas, StErr, Pvals, row.names = "CD")

Model = glm(psoriasis ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
psoriasismerge = data.frame(Betas, StErr, Pvals, row.names = "psoriasis")

Model = glm(UC ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
UCmerge = data.frame(Betas, StErr, Pvals, row.names = "UC")

Model = glm(T1D ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
T1Dmerge = data.frame(Betas, StErr, Pvals, row.names = "T1D")

Model = glm(Ankspond ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = BirthControlStart, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Ankspondmerge = data.frame(Betas, StErr, Pvals, row.names = "Ankspond")

StatisticsBC = rbind(lupusmerge, hashimotomerge,gravesmerge,MSmerge, RAmerge, PBCmerge,
                            Alzheimersmerge,Celiacmerge,PAmerge,MGmerge,CDmerge,psoriasismerge,UCmerge,T1Dmerge,
                            Ankspondmerge)

write.table(Statistics, "UKBBDiseaseSexRatio.txt", sep = "\t", quote = F)
write.table(StatisticsBC, "UKBBBirthControlStats.txt", sep = "\t", quote = F)
write.table(StatisticsBirths, "UKBBBirthsStats.txt", sep = "\t", quote = F)
write.table(StatisticsMenarche, "UKBBMenarcheStats.txt", sep = "\t", quote = F)
write.table(StatisticsMenopause, "UKBBMenopauseStats.txt", sep = "\t", quote = F)
write.table(StatisticsPregnancies, "UKBBMenopausePregnancies.txt", sep = "\t", quote = F)

Cycling = read.delim("PCH_more_cycling_participant.tsv")

Women = subset(Cycling, X31.0.0 == 0)
Women$Autoimmune <- ifelse(grepl("M32|E06|E05|G35|M05|K743|G30|K900|L405|G700|K50|L40|K51|E10|M45", Women$X41270.0.0), 1, 0)

EverOC<-subset(Women, X2784.0.0!=-1)
EverOC<-subset(EverOC, X2784.0.0!=-3)

Model = glm(Autoimmune ~ X2784.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = EverOC, family = binomial)
coef(summary(Model))

FirstBirth<-subset(Women, X2754.0.0!=-1)
FirstBirth<-subset(FirstBirth, X2754.0.0!=-3)

Model = glm(Autoimmune ~ X2754.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
coef(summary(Model))

Model = glm(graves ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
gravesmerge = data.frame(Betas, StErr, Pvals, row.names = "graves")

Model = glm(lupus ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
lupusmerge = data.frame(Betas, StErr, Pvals, row.names = "lupus")

Model = glm(hashimoto ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
hashimotomerge = data.frame(Betas, StErr, Pvals, row.names = "hashimoto")

Model = glm(graves ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MSmerge = data.frame(Betas, StErr, Pvals, row.names = "MS")

Model = glm(RA ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
RAmerge = data.frame(Betas, StErr, Pvals, row.names = "RA")

Model = glm(PBC ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PBCmerge = data.frame(Betas, StErr, Pvals, row.names = "PBC")

Model = glm(Alzheimers ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Alzheimersmerge = data.frame(Betas, StErr, Pvals, row.names = "Alzheimers")

Model = glm(Celiac ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Celiacmerge = data.frame(Betas, StErr, Pvals, row.names = "Celiac")

Model = glm(PA ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
PAmerge = data.frame(Betas, StErr, Pvals, row.names = "PA")

Model = glm(MG ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
MGmerge = data.frame(Betas, StErr, Pvals, row.names = "MG")

Model = glm(CD ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
CDmerge = data.frame(Betas, StErr, Pvals, row.names = "CD")

Model = glm(psoriasis ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
psoriasismerge = data.frame(Betas, StErr, Pvals, row.names = "psoriasis")

Model = glm(UC ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
UCmerge = data.frame(Betas, StErr, Pvals, row.names = "UC")

Model = glm(T1D ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
T1Dmerge = data.frame(Betas, StErr, Pvals, row.names = "T1D")

Model = glm(Ankspond ~ X2794.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = FirstBirth, family = binomial)
Betas = coef(summary(Model))[2,1]
StErr = coef(summary(Model))[2,2]
Pvals = coef(summary(Model))[2,4]
Ankspondmerge = data.frame(Betas, StErr, Pvals, row.names = "Ankspond")

StatisticsBC = rbind(lupusmerge, hashimotomerge,gravesmerge,MSmerge, RAmerge, PBCmerge,
                     Alzheimersmerge,Celiacmerge,PAmerge,MGmerge,CDmerge,psoriasismerge,UCmerge,T1Dmerge,
                     Ankspondmerge)



LastBirth<-subset(Women, X2764.0.0!=-1)
LastBirth<-subset(LastBirth, X2764.0.0!=-3)

Model = glm(Autoimmune ~ X2764.0.0 + X34.0.0 + X21000.0.0 + X22189.0.0, data = LastBirth, family = binomial)
coef(summary(Model))

LastBirth<-subset(LastBirth, X2754.0.0!=-3)
LastBirth<-subset(LastBirth, X2754.0.0!=-1)

LastBirth$BirthTime = LastBirth$X2764.0.0-LastBirth$X2754.0.0
Model = glm(Autoimmune ~ BirthTime + X34.0.0 + X21000.0.0 + X22189.0.0, data = LastBirth, family = binomial)
coef(summary(Model))



