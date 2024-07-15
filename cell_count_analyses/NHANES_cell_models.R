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

#filter out +/- 5 standard deviations for each cell type
SDOver <- mean(CellDemogAge$LBXLYPCT) + sd(CellDemogAge$LBXLYPCT)*5 
SDUnder <- mean(CellDemogAge$LBXLYPCT) - sd(CellDemogAge$LBXLYPCT)*5 

CellDemogAge1 <- CellDemogAge %>% filter((LBXLYPCT<SDOver)&LBXLYPCT>SDUnder)

SDOver <- mean(CellDemogAge$LBXMOPCT) + sd(CellDemogAge$LBXMOPCT)*5 
SDUnder <- mean(CellDemogAge$LBXMOPCT) - sd(CellDemogAge$LBXMOPCT)*5 

CellDemogAge2 <- CellDemogAge %>% filter((LBXMOPCT<SDOver)&LBXMOPCT>SDUnder)

SDOver <- mean(CellDemogAge$LBXNEPCT) + sd(CellDemogAge$LBXNEPCT)*5 
SDUnder <- mean(CellDemogAge$LBXNEPCT) - sd(CellDemogAge$LBXNEPCT)*5 

CellDemogAge3 <- CellDemogAge %>% filter((LBXNEPCT<SDOver)&LBXNEPCT>SDUnder)

SDOver <- mean(CellDemogAge$LBXBAPCT) + sd(CellDemogAge$LBXBAPCT)*5 
SDUnder <- mean(CellDemogAge$LBXBAPCT) - sd(CellDemogAge$LBXBAPCT)*5 

CellDemogAge4 <- CellDemogAge %>% filter((LBXBAPCT<SDOver)&LBXBAPCT>SDUnder)

SDOver <- mean(CellDemogAge$LBXEOPCT) + sd(CellDemogAge$LBXEOPCT)*5 
SDUnder <- mean(CellDemogAge$LBXEOPCT) - sd(CellDemogAge$LBXEOPCT)*5 

CellDemogAge5 <- CellDemogAge %>% filter((LBXEOPCT<SDOver)&LBXEOPCT>SDUnder)

CellDemogAgeTogether <- merge(CellDemogAge1, CellDemogAge2)
CellDemogAgeTogether <- merge(CellDemogAgeTogether, CellDemogAge3)
CellDemogAgeTogether <- merge(CellDemogAgeTogether, CellDemogAge4)
CellDemogAgeTogether <- merge(CellDemogAgeTogether, CellDemogAge5)

#change name of gender 
CellDemogAgeTogether <- CellDemogAgeTogether %>%
  mutate(RIAGENDR = relevel(RIAGENDR, ref = "Female"))

ScaleCell <- CellDemogAgeTogether

#scale cell data
ScaleCell[c(7,8,9,10,11,12)] <- scale(CellDemogAgeTogether[c(7,8,9,10,11,12)])

#linear models on cell percentages
Lymphocytes = glm( LBXLYPCT ~ RIAGENDR +  RIDAGEYR, data = ScaleCell, family = gaussian)

Betas = coef(summary(Lymphocytes))[2,1]
StErr = coef(summary(Lymphocytes))[2,2]
Pvals = coef(summary(Lymphocytes))[2,4]

NHANESLymphocytes = data.frame(Betas, StErr, Pvals, row.names = "NHANESLymphocytes")

Segmented.neutrophils = glm( LBXNEPCT ~ RIAGENDR +  RIDAGEYR, data = ScaleCell, family = gaussian)

Betas = coef(summary(Segmented.neutrophils))[2,1]
StErr = coef(summary(Segmented.neutrophils))[2,2]
Pvals = coef(summary(Segmented.neutrophils))[2,4]

NHANESSegmented.neutrophils = data.frame(Betas, StErr, Pvals, row.names = "NHANESSegmented.neutrophils")

Monocytes = glm( LBXMOPCT ~ RIAGENDR +  RIDAGEYR, data = ScaleCell, family = gaussian)

Betas = coef(summary(Monocytes))[2,1]
StErr = coef(summary(Monocytes))[2,2]
Pvals = coef(summary(Monocytes))[2,4]

NHANESMonocytes = data.frame(Betas, StErr, Pvals, row.names = "NHANESMonocytes")

Basophils = glm( LBXBAPCT ~ RIAGENDR +  RIDAGEYR, data = ScaleCell, family = gaussian)

Betas = coef(summary(Basophils))[2,1]
StErr = coef(summary(Basophils))[2,2]
Pvals = coef(summary(Basophils))[2,4]

NHANESBasophils = data.frame(Betas, StErr, Pvals, row.names = "NHANESBasophils")

Eosinophils = glm( LBXEOPCT ~ RIAGENDR +  RIDAGEYR, data = ScaleCell, family = gaussian)

Betas = coef(summary(Eosinophils))[2,1]
StErr = coef(summary(Eosinophils))[2,2]
Pvals = coef(summary(Eosinophils))[2,4]

NHANESEosinophils = data.frame(Betas, StErr, Pvals, row.names = "NHANESEosinophils")

Statisticsdf = rbind(NHANESBasophils, NHANESEosinophils, NHANESSegmented.neutrophils, NHANESMonocytes, NHANESLymphocytes)
write.table(Statisticsdf, "NHANESCellModelStatistics.txt", quote = F)


