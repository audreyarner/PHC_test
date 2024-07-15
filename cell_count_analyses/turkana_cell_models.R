setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")
library(dplyr)

turkana = read.delim("Turkana_lifestyle_metrics.txt", sep = "\t")

turkana_cellcount = turkana[!is.na(turkana$cell_count_total),]
turkana_cellcount = subset(turkana_cellcount, cell_count_total != 0)
turkana_cellcount = turkana_cellcount[!is.na(turkana_cellcount$lifestyle),]


#select only one measurement if individual repeated
UniqueTurkanaBlood <- turkana_cellcount[!duplicated(turkana_cellcount$unique_id), ]

#check for outliers
SDOver <- mean(UniqueTurkanaBlood$monocyte_perc) + sd(UniqueTurkanaBlood$monocyte_perc)*5 
SDUnder <- mean(UniqueTurkanaBlood$monocyte_perc) - sd(UniqueTurkanaBlood$monocyte_perc)*5 

UniqueTurkanaBlood1 <- UniqueTurkanaBlood %>% filter((monocyte_perc<SDOver)&(monocyte_perc>SDUnder))

SDOver <- mean(UniqueTurkanaBlood$neutrophil_perc) + sd(UniqueTurkanaBlood$neutrophil_perc)*5 
SDUnder <- mean(UniqueTurkanaBlood$neutrophil_perc) - sd(UniqueTurkanaBlood$neutrophil_perc)*5 

UniqueTurkanaBlood2 <- UniqueTurkanaBlood %>% filter((neutrophil_perc<SDOver)&(neutrophil_perc>SDUnder))

SDOver <- mean(UniqueTurkanaBlood$lymphocyte_perc) + sd(UniqueTurkanaBlood$neutrophil_perc)*5 
SDUnder <- mean(UniqueTurkanaBlood$neutrophil_perc) - sd(UniqueTurkanaBlood$neutrophil_perc)*5 

UniqueTurkanaBlood3 <- UniqueTurkanaBlood %>% filter((neutrophil_perc<SDOver)&(neutrophil_perc>SDUnder))

SDOver <- mean(UniqueTurkanaBlood$basophil_perc) + sd(UniqueTurkanaBlood$basophil_perc)*5 
SDUnder <- mean(UniqueTurkanaBlood$basophil_perc) - sd(UniqueTurkanaBlood$basophil_perc)*5 

UniqueTurkanaBlood4 <- UniqueTurkanaBlood %>% filter((basophil_perc<SDOver)&(basophil_perc>SDUnder))

SDOver <- mean(UniqueTurkanaBlood$eosinophil_perc) + sd(UniqueTurkanaBlood$eosinophil_perc)*5 
SDUnder <- mean(UniqueTurkanaBlood$eosinophil_perc) - sd(UniqueTurkanaBlood$eosinophil_perc)*5 

UniqueTurkanaBlood5 <- UniqueTurkanaBlood %>% filter((eosinophil_perc<SDOver)&(eosinophil_perc>SDUnder))

UniqueTurkanaBlood <- merge(UniqueTurkanaBlood,UniqueTurkanaBlood4)

#scale cell proportions
scale_turkana_blood <- UniqueTurkanaBlood
scale_turkana_blood[c(176,177,178,179,180)] <- scale(scale_turkana_blood[c(176,177,178,179,180)])

scale_turkana_blood = subset(scale_turkana_blood, lifestyle != "PeriUrban")

#models to see if lifestyle predicts cell count proportions
Eosinophils = glm(eosinophil_perc ~ lifestyle + age + gender, data = scale_turkana_blood, family = gaussian)
coef(summary(Eosinophils))[2,4]
Monocytes = glm(monocyte_perc ~ lifestyle + age + gender, data = scale_turkana_blood, family = gaussian)
coef(summary(Monocytes))[2,4]
Basophils = glm(basophil_perc ~ lifestyle + age + gender, data = scale_turkana_blood, family = gaussian)
coef(summary(Basophils))[2,4]
Segmented.neutrophils = glm(neutrophil_perc ~ lifestyle + age + gender, data = scale_turkana_blood, family = gaussian)
coef(summary(Segmented.neutrophils))[2,4]
Lymphocytes = glm(lymphocyte_perc ~ lifestyle + age + gender, data = scale_turkana_blood, family = gaussian)
coef(summary(Lymphocytes))[2,4]

#GLM on interaction of lifestyle and gender
Eosinophils = glm(eosinophil_perc ~ gender*lifestyle + age, data = scale_turkana_blood, family = gaussian)

Betas = coef(summary(Eosinophils))[5,1]
StErr = coef(summary(Eosinophils))[5,2]
Pvals = coef(summary(Eosinophils))[5,4]

TurkanaEosinophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaEosinophils")

Monocytes = glm(monocyte_perc ~ gender*lifestyle + age, data = scale_turkana_blood, family = gaussian)

Betas = coef(summary(Monocytes))[5,1]
StErr = coef(summary(Monocytes))[5,2]
Pvals = coef(summary(Monocytes))[5,4]

TurkanaMonocytes = data.frame(Betas, StErr, Pvals, row.names = "TurkanaMonocytes")

Basophils = glm(basophil_perc ~ gender*lifestyle + age, data = scale_turkana_blood, family = gaussian)

Betas = coef(summary(Basophils))[5,1]
StErr = coef(summary(Basophils))[5,2]
Pvals = coef(summary(Basophils))[5,4]

TurkanaBasophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaBasophils")

Segmented.neutrophils = glm(neutrophil_perc ~ gender*lifestyle + age, data = scale_turkana_blood, family = gaussian)

Betas = coef(summary(Segmented.neutrophils))[5,1]
StErr = coef(summary(Segmented.neutrophils))[5,2]
Pvals = coef(summary(Segmented.neutrophils))[5,4]

TurkanaSegmented.neutrophils = data.frame(Betas, StErr, Pvals, row.names = "TurkanaSegmented.neutrophils")

Lymphocytes = glm(neutrophil_perc ~ gender*lifestyle + age, data = scale_turkana_blood, family = gaussian)

Betas = coef(summary(Lymphocytes))[5,1]
StErr = coef(summary(Lymphocytes))[5,2]
Pvals = coef(summary(Lymphocytes))[5,4]

TurkanaLymphocytes = data.frame(Betas, StErr, Pvals, row.names = "TurkanaLymphocytes")

Statisticsdf = rbind(TurkanaBasophils, TurkanaEosinophils, TurkanaSegmented.neutrophils, TurkanaMonocytes, TurkanaLymphocytes)

Statisticsdf$FDR = p.adjust(Statisticsdf$Pvals, method = "fdr")
View(Statisticsdf)

#write.table(Statisticsdf, "TurkanaInteractionModelStatistics.txt", quote = F)
