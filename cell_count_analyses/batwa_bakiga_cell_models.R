library("readxl")

#read data from Harrison et al
my_data <- read_excel("HarrisonS1_SampleInfo.xlsx")

#only include unique samples
UniqueTwaKiga <- my_data[!duplicated(my_data$Genotyping_ID), ]

#take out individuals over 5 standard deviations away from the mean
SDOver <- mean(UniqueTwaKiga$CD3) + sd(UniqueTwaKiga$CD3)*5 
SDUnder <- mean(UniqueTwaKiga$CD3) - sd(UniqueTwaKiga$CD3)*5 

UniqueTwaKiga1 <- UniqueTwaKiga %>% filter((CD3<SDOver)&(CD3>SDUnder))

SDOver <- mean(UniqueTwaKiga$CD4) + sd(UniqueTwaKiga$CD4)*5 
SDUnder <- mean(UniqueTwaKiga$CD4) - sd(UniqueTwaKiga$CD4)*5 

UniqueTwaKiga2 <- UniqueTwaKiga %>% filter((CD4<SDOver)&(CD4>SDUnder))

SDOver <- mean(UniqueTwaKiga$CD8) + sd(UniqueTwaKiga$CD8)*5 
SDUnder <- mean(UniqueTwaKiga$CD8) - sd(UniqueTwaKiga$CD8)*5 

UniqueTwaKiga3 <- UniqueTwaKiga %>% filter((CD8<SDOver)&(CD8>SDUnder))

SDOver <- mean(UniqueTwaKiga$CD14) + sd(UniqueTwaKiga$CD14)*5 
SDUnder <- mean(UniqueTwaKiga$CD14) - sd(UniqueTwaKiga$CD14)*5 

UniqueTwaKiga4 <- UniqueTwaKiga %>% filter((CD14<SDOver)&(CD14>SDUnder))

SDOver <- mean(UniqueTwaKiga$CD20) + sd(UniqueTwaKiga$CD20)*5 
SDUnder <- mean(UniqueTwaKiga$CD20) - sd(UniqueTwaKiga$CD20)*5 

UniqueTwaKiga5 <- UniqueTwaKiga %>% filter((CD20<SDOver)&(CD20>SDUnder))

SDOver <- mean(UniqueTwaKiga$CD56) + sd(UniqueTwaKiga$CD56)*5 
SDUnder <- mean(UniqueTwaKiga$CD56) - sd(UniqueTwaKiga$CD56)*5 

UniqueTwaKiga6 <- UniqueTwaKiga %>% filter((CD56<SDOver)&(CD56>SDUnder))

UniqueTwaKigaTogether <- merge(UniqueTwaKiga1, UniqueTwaKiga2)
UniqueTwaKigaTogether <- merge(UniqueTwaKigaTogether, UniqueTwaKiga3)
UniqueTwaKigaTogether <- merge(UniqueTwaKigaTogether, UniqueTwaKiga4)
UniqueTwaKigaTogether <- merge(UniqueTwaKigaTogether, UniqueTwaKiga5)
UniqueTwaKigaTogether <- merge(UniqueTwaKigaTogether, UniqueTwaKiga6)

ScaleTwaKiga <- UniqueTwaKigaTogether

#scale cell count data
ScaleTwaKiga[c(9,10,11,12,13,14)] <- scale(UniqueTwaKiga[c(9,10,11,12,13,14)])

#models using guassian distribution
CD3 = glm( CD3 ~ Sex +  Population, data = ScaleTwaKiga, family = gaussian)

Betas = coef(summary(CD3))[2,1]
StErr = coef(summary(CD3))[2,2]
Pvals = coef(summary(CD3))[2,4]

TwaKigaCD3 = data.frame(Betas, StErr, Pvals, row.names = "TwaKigaCD3")

CD4 = glm( CD4 ~ Sex +  Population, data = ScaleTwaKiga, family = gaussian)

Betas = coef(summary(CD4))[2,1]
StErr = coef(summary(CD4))[2,2]
Pvals = coef(summary(CD4))[2,4]

TwaKigaCD4 = data.frame(Betas, StErr, Pvals, row.names = "TwaKigaCD4")

CD8 = glm( CD8 ~ Sex +  Population, data = ScaleTwaKiga, family = gaussian)

Betas = coef(summary(CD8))[2,1]
StErr = coef(summary(CD8))[2,2]
Pvals = coef(summary(CD8))[2,4]

TwaKigaCD8 = data.frame(Betas, StErr, Pvals, row.names = "TwaKigaCD8")

CD14 = glm( CD14 ~ Sex +  Population, data = ScaleTwaKiga, family = gaussian)

Betas = coef(summary(CD14))[2,1]
StErr = coef(summary(CD14))[2,2]
Pvals = coef(summary(CD14))[2,4]

TwaKigaCD14 = data.frame(Betas, StErr, Pvals, row.names = "TwaKigaCD14")

CD20 = glm( CD20 ~ Sex +  Population, data = ScaleTwaKiga, family = gaussian)

Betas = coef(summary(CD20))[2,1]
StErr = coef(summary(CD20))[2,2]
Pvals = coef(summary(CD20))[2,4]

TwaKigaCD20 = data.frame(Betas, StErr, Pvals, row.names = "TwaKigaCD20")

CD56 = glm( CD56 ~ Sex + Population, data = ScaleTwaKiga, family = gaussian)

Betas = coef(summary(CD56))[2,1]
StErr = coef(summary(CD56))[2,2]
Pvals = coef(summary(CD56))[2,4]

TwaKigaCD56 = data.frame(Betas, StErr, Pvals, row.names = "TwaKigaCD56")

Statisticsdf = rbind(TwaKigaCD3, TwaKigaCD4, TwaKigaCD8, TwaKigaCD20, TwaKigaCD14, TwaKigaCD56)

write.table(Statisticsdf, "BatwaBakigaCellModelStatistics.txt", quote = F)

