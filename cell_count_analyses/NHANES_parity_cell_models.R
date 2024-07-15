setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")

library(nhanesA)
library(dplyr)
library(ggplot2)
library(tidyr)

#import CBC data and demographic information
cbc_d <- nhanes("CBC_J")

demo_d <- nhanes("DEMO_J")

preg_d <- nhanes("RHQ_J")

#change so that coded info is set in words
demo_d <- nhanesTranslate('DEMO_J', 'RIAGENDR', data=demo_d)
demo_d <- nhanesTranslate('DEMO_J', 'RIDEXPRG', data=demo_d)

#omit lines with NA
cbc <- na.omit(cbc_d)

#select only columns needed
demog <- demo_d %>% select("SEQN", "RIAGENDR", "RIDEXPRG", 'RIDAGEYR', "RIDAGEMN")

#merge biography and cell count data
CellDemogData <- merge(demog, cbc)

#filter out all individuals younger than 18
CellDemogAge <- CellDemogData %>% filter(RIDAGEYR>=18)

#Only include females
Females <- CellDemogAge %>% filter(RIAGENDR=="Female")

#collect pregnancy info and take out individuals without data
ParityInfo <- preg_d %>% select(SEQN, RHQ160)
ParityInfo <- na.omit(ParityInfo)

FemaleParity <- merge(Females, ParityInfo)

#filter out +/- 5 standard deviations for each cell type
SDOver <- mean(FemaleParity$LBXLYPCT) + sd(FemaleParity$LBXLYPCT)*5 
SDUnder <- mean(FemaleParity$LBXLYPCT) - sd(FemaleParity$LBXLYPCT)*5 

FemaleParity1 <- FemaleParity %>% filter((LBXLYPCT<SDOver)&LBXLYPCT>SDUnder)

SDOver <- mean(FemaleParity$LBXMOPCT) + sd(FemaleParity$LBXMOPCT)*5 
SDUnder <- mean(FemaleParity$LBXMOPCT) - sd(FemaleParity$LBXMOPCT)*5 

FemaleParity2 <- FemaleParity %>% filter((LBXMOPCT<SDOver)&LBXMOPCT>SDUnder)

SDOver <- mean(FemaleParity$LBXNEPCT) + sd(FemaleParity$LBXNEPCT)*5 
SDUnder <- mean(FemaleParity$LBXNEPCT) - sd(FemaleParity$LBXNEPCT)*5 

FemaleParity3 <- FemaleParity %>% filter((LBXNEPCT<SDOver)&LBXNEPCT>SDUnder)

SDOver <- mean(FemaleParity$LBXBAPCT) + sd(FemaleParity$LBXBAPCT)*5 
SDUnder <- mean(FemaleParity$LBXBAPCT) - sd(FemaleParity$LBXBAPCT)*5 

FemaleParity4 <- FemaleParity %>% filter((LBXBAPCT<SDOver)&LBXBAPCT>SDUnder)

SDOver <- mean(FemaleParity$LBXEOPCT) + sd(FemaleParity$LBXEOPCT)*5 
SDUnder <- mean(FemaleParity$LBXEOPCT) - sd(FemaleParity$LBXEOPCT)*5 

FemaleParity5 <- FemaleParity %>% filter((LBXEOPCT<SDOver)&LBXEOPCT>SDUnder)

FemaleParityTogether <- merge(FemaleParity1, FemaleParity2)
FemaleParityTogether <- merge(FemaleParityTogether, FemaleParity3)
FemaleParityTogether <- merge(FemaleParityTogether, FemaleParity4)
FemaleParityTogether <- merge(FemaleParityTogether, FemaleParity5)

#scale cell count and parity data
ScaleParity <- FemaleParityTogether

#scale cell data
ScaleParity[c(7,8,9,10,11,12,27)] <- scale(FemaleParityTogether[c(7,8,9,10,11,12,27)])

#parity models
#linear models on cell percentages
Lymphocytes = glm( LBXLYPCT ~ RHQ160 +  RIDAGEYR, data = ScaleParity, family = gaussian)

Betas = coef(summary(Lymphocytes))[2,1]
StErr = coef(summary(Lymphocytes))[2,2]
Pvals = coef(summary(Lymphocytes))[2,4]

NHANESLymphocytes = data.frame(Betas, StErr, Pvals, row.names = "NHANESLymphocytes")

Segmented.neutrophils = glm( LBXNEPCT ~ RHQ160 +  RIDAGEYR, data = ScaleParity, family = gaussian)

Betas = coef(summary(Segmented.neutrophils))[2,1]
StErr = coef(summary(Segmented.neutrophils))[2,2]
Pvals = coef(summary(Segmented.neutrophils))[2,4]

NHANESSegmented.neutrophils = data.frame(Betas, StErr, Pvals, row.names = "NHANESSegmented.neutrophils")

Monocytes = glm( LBXMOPCT ~ RHQ160 +  RIDAGEYR, data = ScaleParity, family = gaussian)

Betas = coef(summary(Monocytes))[2,1]
StErr = coef(summary(Monocytes))[2,2]
Pvals = coef(summary(Monocytes))[2,4]

NHANESMonocytes = data.frame(Betas, StErr, Pvals, row.names = "NHANESMonocytes")

Basophils = glm( LBXBAPCT ~ RHQ160 +  RIDAGEYR, data = ScaleParity, family = gaussian)

Betas = coef(summary(Basophils))[2,1]
StErr = coef(summary(Basophils))[2,2]
Pvals = coef(summary(Basophils))[2,4]

NHANESBasophils = data.frame(Betas, StErr, Pvals, row.names = "NHANESBasophils")

Eosinophils = glm( LBXEOPCT ~ RHQ160 +  RIDAGEYR, data = ScaleParity, family = gaussian)

Betas = coef(summary(Eosinophils))[2,1]
StErr = coef(summary(Eosinophils))[2,2]
Pvals = coef(summary(Eosinophils))[2,4]

NHANESEosinophils = data.frame(Betas, StErr, Pvals, row.names = "NHANESEosinophils")

Statisticsdf = rbind(NHANESBasophils, NHANESEosinophils, NHANESSegmented.neutrophils, NHANESMonocytes, NHANESLymphocytes)

write.table(Statisticsdf, "NHANESParityCellModelStatistics.txt", quote = F)
