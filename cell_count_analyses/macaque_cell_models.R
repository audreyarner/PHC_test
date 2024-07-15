setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")
library(tidyverse)

data <- read.delim("flow_dataset_Sept2021.txt")

#eliminate a few columns and make batch column the last one in the dataset
data=data[,c(-9, -15:-17)]
data=data[,c(1:12, 14:16, 13)]
data$batch=as.factor(data$batch)

#make ratio of total cell count
data$TotalCells <- data$cd3tcells + data$bcells + data$nkcells + data$cd8tregs + data$cd4tcells + data$cd8tcells + data$tregs + data$cd16hladrpos + data$cd14cd16hladrpos + data$cd14hladrpos
data=mutate(data,cd3tcellsPercent=data$cd3tcells/data$TotalCells)
data=mutate(data,bcellsPercent=data$bcells/data$TotalCells)
data=mutate(data,nkcellsPercent=data$nkcells/data$TotalCells)
data=mutate(data,cd8tregsPercent=data$cd8tregs/data$TotalCells)
data=mutate(data,cd4tcellsPercent=data$cd4tcells/data$TotalCells)
data=mutate(data,cd8tcellsPercent=data$cd8tcells/data$TotalCells)
data=mutate(data,tregsPercent=data$tregs/data$TotalCells)
data=mutate(data,cd16hladrposPercent=data$cd16hladrpos/data$TotalCells)
data=mutate(data,cd14hladrposPercent=data$cd14hladrpos/data$TotalCells)
data=mutate(data,cd14cd16hladrposPercent=data$cd14cd16hladrpos/data$TotalCells)

#scale data
ScaleMacaque <- data

ScaleMacaque[c(18:27)] <- scale(data[c(18:27)])

#calculate over/under 5 standard deviations
SDOver <- mean(ScaleMacaque$cd3tcellsPercent) + sd(ScaleMacaque$cd3tcellsPercent)*5 
SDUnder <- mean(ScaleMacaque$cd3tcellsPercent) - sd(ScaleMacaque$cd3tcellsPercent)*5 

ScaleMacaque1 <- ScaleMacaque %>% filter((cd3tcellsPercent<SDOver)&(cd3tcellsPercent>SDUnder))

SDOver <- mean(ScaleMacaque$bcells) + sd(ScaleMacaque$bcells)*5 
SDUnder <- mean(ScaleMacaque$bcells) - sd(ScaleMacaque$bcells)*5 

ScaleMacaque2 <- ScaleMacaque %>% filter((bcells<SDOver)&(bcells>SDUnder))

SDOver <- mean(ScaleMacaque$nkcells) + sd(ScaleMacaque$nkcells)*5 
SDUnder <- mean(ScaleMacaque$nkcells) - sd(ScaleMacaque$nkcells)*5 

ScaleMacaque3 <- ScaleMacaque %>% filter((nkcells<SDOver)&(nkcells>SDUnder))

SDOver <- mean(ScaleMacaque$cd8tregs) + sd(ScaleMacaque$cd8tregs)*5 
SDUnder <- mean(ScaleMacaque$cd8tregs) - sd(ScaleMacaque$cd8tregs)*5 

ScaleMacaque4 <- ScaleMacaque %>% filter((cd8tregs<SDOver)&(cd8tregs>SDUnder))

SDOver <- mean(ScaleMacaque$cd4tcells) + sd(ScaleMacaque$cd4tcells)*5 
SDUnder <- mean(ScaleMacaque$cd4tcells) - sd(ScaleMacaque$cd4tcells)*5 

ScaleMacaque5 <- ScaleMacaque %>% filter((cd4tcells<SDOver)&(cd4tcells>SDUnder))

SDOver <- mean(ScaleMacaque$cd8tcells) + sd(ScaleMacaque$cd8tcells)*5 
SDUnder <- mean(ScaleMacaque$cd8tcells) - sd(ScaleMacaque$cd8tcells)*5 

ScaleMacaque6 <- ScaleMacaque %>% filter((cd8tcells<SDOver)&(cd8tcells>SDUnder))

SDOver <- mean(ScaleMacaque$tregs) + sd(ScaleMacaque$tregs)*5 
SDUnder <- mean(ScaleMacaque$tregs) - sd(ScaleMacaque$tregs)*5 

ScaleMacaque7 <- ScaleMacaque %>% filter((tregs<SDOver)&(tregs>SDUnder))

SDOver <- mean(ScaleMacaque$cd16hladrpos) + sd(ScaleMacaque$cd16hladrpos)*5 
SDUnder <- mean(ScaleMacaque$cd16hladrpos) - sd(ScaleMacaque$cd16hladrpos)*5 

ScaleMacaque8 <- ScaleMacaque %>% filter((cd16hladrpos<SDOver)&(cd16hladrpos>SDUnder))

SDOver <- mean(ScaleMacaque$cd14hladrpos) + sd(ScaleMacaque$cd14hladrpos)*5 
SDUnder <- mean(ScaleMacaque$cd14hladrpos) - sd(ScaleMacaque$cd14hladrpos)*5 

ScaleMacaque9 <- ScaleMacaque %>% filter((cd14hladrpos<SDOver)&(cd14hladrpos>SDUnder))

SDOver <- mean(ScaleMacaque$cd14cd16hladrpos) + sd(ScaleMacaque$cd14cd16hladrpos)*5 
SDUnder <- mean(ScaleMacaque$cd14cd16hladrpos) - sd(ScaleMacaque$cd14cd16hladrpos)*5 

ScaleMacaque10 <- ScaleMacaque %>% filter((cd14cd16hladrpos<SDOver)&(cd14cd16hladrpos>SDUnder))

ScaleMacaqueTogether <- merge(ScaleMacaque1, ScaleMacaque2)
ScaleMacaqueTogether <- merge(ScaleMacaqueTogether, ScaleMacaque3)
ScaleMacaqueTogether <- merge(ScaleMacaqueTogether, ScaleMacaque4)
ScaleMacaqueTogether <- merge(ScaleMacaqueTogether, ScaleMacaque5)
ScaleMacaqueTogether <- merge(ScaleMacaqueTogether, ScaleMacaque6)
ScaleMacaqueTogether <- merge(ScaleMacaqueTogether, ScaleMacaque7)
ScaleMacaqueTogether <- merge(ScaleMacaqueTogether, ScaleMacaque8)
ScaleMacaqueTogether <- merge(ScaleMacaqueTogether, ScaleMacaque9)
ScaleMacaqueTogether <- merge(ScaleMacaqueTogether, ScaleMacaque10)

#linear models for each cell type
cd3tcells = glm(cd3tcellsPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(cd3tcells))[2,1]
StErr = coef(summary(cd3tcells))[2,2]
Pvals = coef(summary(cd3tcells))[2,4]

Macaquecd3tcells = data.frame(Betas, StErr, Pvals, row.names = "Macaquecd3tcells")

bcells = glm(bcellsPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(bcells))[2,1]
StErr = coef(summary(bcells))[2,2]
Pvals = coef(summary(bcells))[2,4]

Macaquebcells = data.frame(Betas, StErr, Pvals, row.names = "Macaquebcells")

nkcells = glm(nkcellsPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(nkcells))[2,1]
StErr = coef(summary(nkcells))[2,2]
Pvals = coef(summary(nkcells))[2,4]

Macaquenkcells = data.frame(Betas, StErr, Pvals, row.names = "Macaquenkcells")

cd8tregs = glm(cd8tregsPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(cd8tregs))[2,1]
StErr = coef(summary(cd8tregs))[2,2]
Pvals = coef(summary(cd8tregs))[2,4]

Macaquecd8tregs = data.frame(Betas, StErr, Pvals, row.names = "Macaquecd8tregs")

cd4tcells = glm(cd4tcellsPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(cd4tcells))[2,1]
StErr = coef(summary(cd4tcells))[2,2]
Pvals = coef(summary(cd4tcells))[2,4]

Macaquecd4tcells = data.frame(Betas, StErr, Pvals, row.names = "Macaquecd4tcells")

cd8tcells = glm(cd8tcellsPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(cd8tcells))[2,1]
StErr = coef(summary(cd8tcells))[2,2]
Pvals = coef(summary(cd8tcells))[2,4]

Macaquecd8tcells = data.frame(Betas, StErr, Pvals, row.names = "Macaquecd8tcells")

tregs = glm(tregsPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(tregs))[2,1]
StErr = coef(summary(tregs))[2,2]
Pvals = coef(summary(tregs))[2,4]

Macaquetregs = data.frame(Betas, StErr, Pvals, row.names = "Macaquetregs")

cd16hladrpos = glm(cd16hladrposPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(cd16hladrpos))[2,1]
StErr = coef(summary(cd16hladrpos))[2,2]
Pvals = coef(summary(cd16hladrpos))[2,4]

Macaquecd16hladrpos = data.frame(Betas, StErr, Pvals, row.names = "Macaquecd16hladrpos")

cd14cd16hladrpos = glm(cd14cd16hladrposPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(cd14cd16hladrpos))[2,1]
StErr = coef(summary(cd14cd16hladrpos))[2,2]
Pvals = coef(summary(cd14cd16hladrpos))[2,4]

Macaquecd14cd16hladrpos = data.frame(Betas, StErr, Pvals, row.names = "Macaquecd14cd16hladrpos")

cd14hladrpos = glm(cd14hladrposPercent ~ sex + age + batch, data = ScaleMacaqueTogether, family = gaussian)

Betas = coef(summary(cd14hladrpos))[2,1]
StErr = coef(summary(cd14hladrpos))[2,2]
Pvals = coef(summary(cd14hladrpos))[2,4]

Macaquecd14hladrpos = data.frame(Betas, StErr, Pvals, row.names = "Macaquecd14hladrpos")

data$TotalCells <- data$cd3tcells + data$bcells + data$nkcells + data$cd8tregs + data$cd4tcells + data$cd8tcells + data$tregs + data$cd16hladrpos + data$cd14cd16hladrpos + data$cd14hladrpos

Statisticsdf = rbind(Macaquecd3tcells, Macaquebcells, Macaquenkcells, Macaquecd8tregs, Macaquecd4tcells, Macaquecd8tcells, Macaquetregs,Macaquecd16hladrpos,Macaquecd14cd16hladrpos,Macaquecd14hladrpos)

write.table(Statisticsdf, "MacaqueCellModelStatistics.txt", quote = F)





