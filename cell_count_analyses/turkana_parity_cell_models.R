setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")

library(dplyr)
library(ggplot2)
library(ggdark)

#datasets needed
Turkana <- read.delim("TurkanaAllInfo.txt", header = T, sep="\t")
BloodCount <- read.csv("BLOOD_SMEAR_COUNTS.csv", header =T)
Parity <- read.csv("FEMALE_REPRODUCTION.csv", header =T)

#select only the columns that we need
Parity <- Parity %>% dplyr::select(Unique.ID, Number.of.pregnancies)

#change NA individuals to 0 pregnancies
Parity[is.na(Parity)] <- 0

#select individuals that have blood count data
TurkanaBlood <- merge(Turkana, BloodCount)

#Only choose unique individuals
UniqueTurkana <- TurkanaBlood[!duplicated(TurkanaBlood$Unique.ID), ]

#make a column with total cells counted
UniqueTurkana=mutate(UniqueTurkana,TotalCells=UniqueTurkana$Monocytes + UniqueTurkana$Segmented.neutrophils + UniqueTurkana$Lymphocytes + UniqueTurkana$Eosinophils + UniqueTurkana$Basophils)

#make separate columns for the percent of each cell type 
UniqueTurkana=mutate(UniqueTurkana,MonocytesPercent=UniqueTurkana$Monocytes/UniqueTurkana$TotalCells)
UniqueTurkana=mutate(UniqueTurkana,Segmented.neutrophilsPercent=UniqueTurkana$Segmented.neutrophils/UniqueTurkana$TotalCells)
UniqueTurkana=mutate(UniqueTurkana,LymphocytesPercent=UniqueTurkana$Lymphocytes/UniqueTurkana$TotalCells)
UniqueTurkana=mutate(UniqueTurkana,EosinophilsPercent=UniqueTurkana$Eosinophils/UniqueTurkana$TotalCells)
UniqueTurkana=mutate(UniqueTurkana,BasophilsPercent=UniqueTurkana$Basophils/UniqueTurkana$TotalCells)

#check for outliers
SDOver <- mean(UniqueTurkana$MonocytesPercent) + sd(UniqueTurkana$MonocytesPercent)*5 
SDUnder <- mean(UniqueTurkana$MonocytesPercent) - sd(UniqueTurkana$MonocytesPercent)*5 

UniqueTurkana1 <- UniqueTurkana %>% filter((MonocytesPercent<SDOver)&(MonocytesPercent>SDUnder))

SDOver <- mean(UniqueTurkana$Segmented.neutrophilsPercent) + sd(UniqueTurkana$Segmented.neutrophilsPercent)*5 
SDUnder <- mean(UniqueTurkana$Segmented.neutrophilsPercent) - sd(UniqueTurkana$Segmented.neutrophilsPercent)*5 

UniqueTurkana2 <- UniqueTurkana %>% filter((Segmented.neutrophilsPercent<SDOver)&(Segmented.neutrophilsPercent>SDUnder))

SDOver <- mean(UniqueTurkana$LymphocytesPercent) + sd(UniqueTurkana$LymphocytesPercent)*5 
SDUnder <- mean(UniqueTurkana$LymphocytesPercent) - sd(UniqueTurkana$LymphocytesPercent)*5 

UniqueTurkana3 <- UniqueTurkana %>% filter((LymphocytesPercent<SDOver)&(LymphocytesPercent>SDUnder))

SDOver <- mean(UniqueTurkana$BasophilsPercent) + sd(UniqueTurkana$BasophilsPercent)*5 
SDUnder <- mean(UniqueTurkana$BasophilsPercent) - sd(UniqueTurkana$BasophilsPercent)*5 

UniqueTurkana4 <- UniqueTurkana %>% filter((BasophilsPercent<SDOver)&(BasophilsPercent>SDUnder))

SDOver <- mean(UniqueTurkana$EosinophilsPercent) + sd(UniqueTurkana$EosinophilsPercent)*5 
SDUnder <- mean(UniqueTurkana$EosinophilsPercent) - sd(UniqueTurkana$EosinophilsPercent)*5 

UniqueTurkana5 <- UniqueTurkana %>% filter((EosinophilsPercent<SDOver)&(EosinophilsPercent>SDUnder))

UniqueTurkana <- merge(UniqueTurkana,UniqueTurkana4)

#make data frame for population combinations interested in
NoMedium <- UniqueTurkana %>% filter(Lifestyle != "RuralNotPastoral")

#analyses with both urban and Urban individuals, merge with parity information
NoMedium2 <- merge(NoMedium, Parity)

#select females
Women <- NoMedium2 %>% filter(Gender =="F")

#scale cell count proportions
ScaleWomen <- Women
ScaleWomen[c(24,25,26,27,28,29)] <- scale(Women[c(24,25,26,27,28,29)])

#linear model to identify if parity can be predicted by lifestyle
ParityGLM = glm(Number.of.pregnancies ~ Lifestyle + Age, data = Women, family = gaussian)

coef(summary(ParityGLM))

p <- ggplot(Women, aes(x=Number.of.pregnancies,y=Lifestyle,  fill = Lifestyle, color = Lifestyle)) + 
  geom_violin(linewidth = 2) + 
  dark_theme_classic() +
  scale_color_manual(values=c("#ff8fa0", "#388e86")) +
  scale_fill_manual(values=c("#ffb6c1", "#43aaa0")) +
  geom_boxplot(width=0.2, linewidth = 2) + 
  theme(legend.position = "none", text = element_text(size = 40)) +
  theme(axis.text.y = element_blank()) 

ggsave("PregnancyViolins.pdf", plot = p, width = 8, height = 4)

#linear model for cell counts ~ parity in all Turkana together
Eosinophils = glm(EosinophilsPercent ~ Number.of.pregnancies + Age + Lifestyle, data = Women, family = gaussian)

Betas = coef(summary(Eosinophils))[2,1]
StErr = coef(summary(Eosinophils))[2,2]
Pvals = coef(summary(Eosinophils))[2,4]

TurkanaEosinophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaEosinophils")

Monocytes = glm(MonocytesPercent ~ Number.of.pregnancies + Lifestyle + Age, data = Women, family = gaussian)

Betas = coef(summary(Monocytes))[2,1]
StErr = coef(summary(Monocytes))[2,2]
Pvals = coef(summary(Monocytes))[2,4]

TurkanaMonocytes = data.frame(Betas, StErr, Pvals, row.names = "TurkanaMonocytes")

Basophils = glm(BasophilsPercent ~ Number.of.pregnancies + Age + Lifestyle, data = Women, family = gaussian)

Betas = coef(summary(Basophils))[2,1]
StErr = coef(summary(Basophils))[2,2]
Pvals = coef(summary(Basophils))[2,4]

TurkanaBasophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaBasophils")

Segmented.neutrophils = glm(Segmented.neutrophilsPercent ~ Number.of.pregnancies + Age + Lifestyle, data = Women, family = gaussian)

Betas = coef(summary(Segmented.neutrophils))[2,1]
StErr = coef(summary(Segmented.neutrophils))[2,2]
Pvals = coef(summary(Segmented.neutrophils))[2,4]

TurkanaSegmented.neutrophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaSegmented.neutrophils")

Lymphocytes = glm(LymphocytesPercent ~ Number.of.pregnancies + Age + Lifestyle, data = Women, family = gaussian)

Betas = coef(summary(Lymphocytes))[2,1]
StErr = coef(summary(Lymphocytes))[2,2]
Pvals = coef(summary(Lymphocytes))[2,4]

TurkanaLymphocytes = data.frame(Betas, StErr, Pvals, row.names = "TurkanaLymphocytes")

Statisticsdf = rbind(TurkanaBasophils, TurkanaEosinophils, TurkanaSegmented.neutrophils, TurkanaMonocytes, TurkanaLymphocytes)

write.table(Statisticsdf, "TurkanaParityModelStatistics.txt", quote = F)

#analyses only in Pastoralist Turkana
#scale cell count proportions
PastoralistWomen <- subset(Women, Lifestyle == "Pastoralist")

ScalePastoralistWomen <- PastoralistWomen
ScalePastoralistWomen[c(24,25,26,27,28,29)] <- scale(PastoralistWomen[c(24,25,26,27,28,29)])

#linear model for cell counts ~ parity in all Turkana together
Eosinophils = glm(EosinophilsPercent ~ Number.of.pregnancies + Age, data = ScalePastoralistWomen, family = gaussian)

Betas = coef(summary(Eosinophils))[2,1]
StErr = coef(summary(Eosinophils))[2,2]
Pvals = coef(summary(Eosinophils))[2,4]

TurkanaEosinophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaEosinophils")

Monocytes = glm(MonocytesPercent ~ Number.of.pregnancies + Age, data = ScalePastoralistWomen, family = gaussian)

Betas = coef(summary(Monocytes))[2,1]
StErr = coef(summary(Monocytes))[2,2]
Pvals = coef(summary(Monocytes))[2,4]

TurkanaMonocytes = data.frame(Betas, StErr, Pvals, row.names = "TurkanaMonocytes")

Basophils = glm(BasophilsPercent ~ Number.of.pregnancies + Age, data = ScalePastoralistWomen, family = gaussian)

Betas = coef(summary(Basophils))[2,1]
StErr = coef(summary(Basophils))[2,2]
Pvals = coef(summary(Basophils))[2,4]

TurkanaBasophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaBasophils")

Segmented.neutrophils = glm(Segmented.neutrophilsPercent ~ Number.of.pregnancies + Age, data = ScalePastoralistWomen, family = gaussian)

Betas = coef(summary(Segmented.neutrophils))[2,1]
StErr = coef(summary(Segmented.neutrophils))[2,2]
Pvals = coef(summary(Segmented.neutrophils))[2,4]

TurkanaSegmented.neutrophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaSegmented.neutrophils")

Lymphocytes = glm(LymphocytesPercent ~ Number.of.pregnancies + Age, data = ScalePastoralistWomen, family = gaussian)

Betas = coef(summary(Lymphocytes))[2,1]
StErr = coef(summary(Lymphocytes))[2,2]
Pvals = coef(summary(Lymphocytes))[2,4]

TurkanaLymphocytes = data.frame(Betas, StErr, Pvals, row.names = "TurkanaLymphocytes")

Statisticsdf = rbind(TurkanaBasophils, TurkanaEosinophils, TurkanaSegmented.neutrophils, TurkanaMonocytes, TurkanaLymphocytes)

write.table(Statisticsdf, "PastoralistTurkanaParityModelStatistics.txt", quote = F)

#analyses only in Urban Turkana
#scale cell count proportions
UrbanWomen <- subset(Women, Lifestyle == "Urban")

ScaleUrbanWomen <- UrbanWomen
ScaleUrbanWomen[c(24,25,26,27,28,29)] <- scale(UrbanWomen[c(24,25,26,27,28,29)])

#linear model for cell counts ~ parity in all Turkana together
Eosinophils = glm(EosinophilsPercent ~ Number.of.pregnancies + Age, data = ScaleUrbanWomen, family = gaussian)

Betas = coef(summary(Eosinophils))[2,1]
StErr = coef(summary(Eosinophils))[2,2]
Pvals = coef(summary(Eosinophils))[2,4]

TurkanaEosinophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaEosinophils")

Monocytes = glm(MonocytesPercent ~ Number.of.pregnancies + Age, data = ScaleUrbanWomen, family = gaussian)

Betas = coef(summary(Monocytes))[2,1]
StErr = coef(summary(Monocytes))[2,2]
Pvals = coef(summary(Monocytes))[2,4]

TurkanaMonocytes = data.frame(Betas, StErr, Pvals, row.names = "TurkanaMonocytes")

Basophils = glm(BasophilsPercent ~ Number.of.pregnancies + Age, data = ScaleUrbanWomen, family = gaussian)

Betas = coef(summary(Basophils))[2,1]
StErr = coef(summary(Basophils))[2,2]
Pvals = coef(summary(Basophils))[2,4]

TurkanaBasophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaBasophils")

Segmented.neutrophils = glm(Segmented.neutrophilsPercent ~ Number.of.pregnancies + Age, data = ScaleUrbanWomen, family = gaussian)

Betas = coef(summary(Segmented.neutrophils))[2,1]
StErr = coef(summary(Segmented.neutrophils))[2,2]
Pvals = coef(summary(Segmented.neutrophils))[2,4]

TurkanaSegmented.neutrophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaSegmented.neutrophils")

Lymphocytes = glm(LymphocytesPercent ~ Number.of.pregnancies + Age, data = ScaleUrbanWomen, family = gaussian)

Betas = coef(summary(Lymphocytes))[2,1]
StErr = coef(summary(Lymphocytes))[2,2]
Pvals = coef(summary(Lymphocytes))[2,4]

TurkanaLymphocytes = data.frame(Betas, StErr, Pvals, row.names = "TurkanaLymphocytes")

Statisticsdf = rbind(TurkanaBasophils, TurkanaEosinophils, TurkanaSegmented.neutrophils, TurkanaMonocytes, TurkanaLymphocytes)

write.table(Statisticsdf, "UrbanTurkanaParityModelStatistics.txt", quote = F)
