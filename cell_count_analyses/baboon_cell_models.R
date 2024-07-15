library(dplyr)
library(ggplot2)
library(tidyr)

#data from Lea et al 2018
Baboon <- read.delim("pnas.1811967115.sd02.txt", header =T)

#only include cell counts from unique individuals
UniqueBaboon <- Baboon[!duplicated(Baboon$Age),]

UniqueBaboon <- na.omit(UniqueBaboon)

#scale baboon data
ScaleBaboon <- UniqueBaboon

ScaleBaboon[c(2,3,4,5,6)] <- scale(UniqueBaboon[c(2,3,4,5,6)])

#take out any indiviudals outside of 5 standard deviations from the mean
SDOver <- mean(ScaleBaboon$CD4proportion) + sd(ScaleBaboon$CD4proportion)*5 
SDUnder <- mean(ScaleBaboon$CD4proportion) - sd(ScaleBaboon$CD4proportion)*5 

ScaleBaboon1 <- ScaleBaboon %>% filter((CD4proportion<SDOver)&(CD4proportion>SDUnder))

SDOver <- mean(ScaleBaboon$CD8proportion) + sd(ScaleBaboon$CD8proportion)*5 
SDUnder <- mean(ScaleBaboon$CD8proportion) - sd(ScaleBaboon$CD8proportion)*5 

ScaleBaboon2 <- ScaleBaboon %>% filter((CD8proportion<SDOver)&(CD8proportion>SDUnder))

SDOver <- mean(ScaleBaboon$CD20proportion) + sd(ScaleBaboon$CD20proportion)*5 
SDUnder <- mean(ScaleBaboon$CD20proportion) - sd(ScaleBaboon$CD20proportion)*5 

ScaleBaboon3 <- ScaleBaboon %>% filter((CD20proportion<SDOver)&(CD20proportion>SDUnder))

SDOver <- mean(ScaleBaboon$CD14proportion) + sd(ScaleBaboon$CD14proportion)*5 
SDUnder <- mean(ScaleBaboon$CD14proportion) - sd(ScaleBaboon$CD14proportion)*5 

ScaleBaboon4 <- ScaleBaboon %>% filter((CD14proportion<SDOver)&(CD14proportion>SDUnder))

SDOver <- mean(ScaleBaboon$CD16proportion) + sd(ScaleBaboon$CD16proportion)*5 
SDUnder <- mean(ScaleBaboon$CD16proportion) - sd(ScaleBaboon$CD16proportion)*5 

ScaleBaboon5 <- ScaleBaboon %>% filter((CD16proportion<SDOver)&(CD16proportion>SDUnder))

ScaleBaboonTogether <- merge(ScaleBaboon1, ScaleBaboon2)
ScaleBaboonTogether <- merge(ScaleBaboonTogether, ScaleBaboon3)
ScaleBaboonTogether <- merge(ScaleBaboonTogether, ScaleBaboon4)
ScaleBaboonTogether <- merge(ScaleBaboonTogether, ScaleBaboon5)

#linear models on cell count proportions
CD4 = glm( CD4proportion ~ Sex +  Age, data = ScaleBaboonTogether, family = gaussian)

Betas = coef(summary(CD4))[2,1]
StErr = coef(summary(CD4))[2,2]
Pvals = coef(summary(CD4))[2,4]

BaboonCD4 = data.frame(Betas, StErr, Pvals, row.names = "BaboonCD4")

CD8 = glm( CD8proportion ~ Sex +  Age, data = ScaleBaboonTogether, family = gaussian)

Betas = coef(summary(CD8))[2,1]
StErr = coef(summary(CD8))[2,2]
Pvals = coef(summary(CD8))[2,4]

BaboonCD8 = data.frame(Betas, StErr, Pvals, row.names = "BaboonCD8")

CD20 = glm( CD20proportion ~ Sex +  Age, data = ScaleBaboonTogether, family = gaussian)

Betas = coef(summary(CD20))[2,1]
StErr = coef(summary(CD20))[2,2]
Pvals = coef(summary(CD20))[2,4]

BaboonCD20 = data.frame(Betas, StErr, Pvals, row.names = "BaboonCD20")

CD14 = glm( CD14proportion ~ Sex +  Age, data = ScaleBaboonTogether, family = gaussian)

Betas = coef(summary(CD14))[2,1]
StErr = coef(summary(CD14))[2,2]
Pvals = coef(summary(CD14))[2,4]

BaboonCD14 = data.frame(Betas, StErr, Pvals, row.names = "BaboonCD14")

CD16 = glm( CD16proportion ~ Sex +  Age, data = ScaleBaboonTogether, family = gaussian)

Betas = coef(summary(CD16))[2,1]
StErr = coef(summary(CD16))[2,2]
Pvals = coef(summary(CD16))[2,4]

BaboonCD16 = data.frame(Betas, StErr, Pvals, row.names = "BaboonCD16")

Statisticsdf = rbind(BaboonCD4, BaboonCD8, BaboonCD20, BaboonCD14, BaboonCD16)

write.table(Statisticsdf, "BaboonCellModelStatistics.txt", quote = F)


