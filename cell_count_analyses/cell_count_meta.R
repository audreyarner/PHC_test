setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")
library(metafor)
library(reshape2)
library(ggdark)

baboon <- read.delim("BaboonCellModelStatistics.txt", sep = " ")
baboon$Lifestyle = "Traditional"
baboon$Immunity = c("Adaptive", "Adaptive", "Innate", "Innate", "Adaptive")

macaque <- read.delim("MacaqueCellModelStatistics.txt", sep=" ")
macaque$Lifestyle = "Traditional"
macaque$Immunity = c("Adaptive", "Adaptive", "Innate", "Adaptive", "Adaptive", "Adaptive", "Adaptive", "Innate", "Innate", "Innate")

NHANES <- read.delim("NHANESCellModelStatistics.txt", sep = " ")
NHANES$Lifestyle = "Novel"
NHANES$Immunity = c("Innate", "Innate", "Innate", "Innate", "Adaptive")

twakiga <- read.delim("BatwaBakigaCellModelStatistics.txt", sep = " ")
twakiga$Lifestyle = "Traditional"
twakiga$Immunity = c("Adaptive", "Adaptive", "Adaptive", "Adaptive", "Innate", "Innate")

UrbanT <- read.delim("UrbanTurkanaCellModelStatistics.txt", sep = " ")
UrbanT$Lifestyle = "Novel"
UrbanT$Immunity = c("Innate", "Innate", "Innate", "Innate", "Adaptive")

PastoralistT <- read.delim("PastoralistTurkanaCellModelStatistics.txt", sep = " ")
PastoralistT$Lifestyle = "Traditional"
PastoralistT$Immunity = c("Innate", "Innate", "Innate", "Innate", "Adaptive")

AllData = rbind(baboon, macaque, NHANES, twakiga, UrbanT, PastoralistT)
AllData$CellType = row.names(AllData)

m.qual <- rma(yi = Betas,
               sei = StErr,
               data = AllData,
               method = "ML",
               mods = ~ Lifestyle*Immunity)
m.qual

m.qual2 <- rma(yi = Betas,
              sei = StErr,
              data = AllData,
              method = "ML",
              mods = ~ Lifestyle+Immunity)
m.qual2

Adaptive = subset(AllData, Immunity == "Adaptive")
m.qual3 <- rma(yi = Betas,
               sei = StErr,
               data = Adaptive,
               method = "ML",
               mods = ~ Lifestyle)
m.qual3

Innate = subset(AllData, Immunity == "Innate")
m.qual4 <- rma(yi = Betas,
               sei = StErr,
               data = Innate,
               method = "ML",
               mods = ~ Lifestyle)
m.qual4

Innate$AbsBetas = abs(Innate$Betas)
Innate$SE1 = Innate$AbsBeta - Innate$StErr
Innate$SE2 = Innate$AbsBeta + Innate$StErr

Innate$CellType <- factor(Innate$CellType,levels = c("BaboonCD20",
                                                     "BaboonCD14",
                                                    "Macaquenkcells",
                                                    "Macaquecd8tregs",              
                                                     "Macaquetregs",
                                                    "Macaquecd16hladrpos",          
                                                     "Macaquecd14cd16hladrpos",
                                                    "Macaquecd14hladrpos",          
                                                    "TwaKigaCD14",
                                                    "TwaKigaCD56",
                                                    "TurkanaBasophils1",
                                                    "TurkanaEosinophils1",          
                                                    "TurkanaSegmented.neutrophils1",
                                                    "TurkanaMonocytes1",
                                                    "TurkanaBasophils",
                                                    "TurkanaEosinophils",          
                                                    "TurkanaSegmented.neutrophils",
                                                    "TurkanaMonocytes",
                                                    "NHANESBasophils",
                                                    "NHANESEosinophils",
                                                    "NHANESSegmented.neutrophils",
                                                    "NHANESMonocytes")) 


Plot <- ggplot(Innate, aes(x=CellType,y=AbsBetas, color = Lifestyle)) + 
  geom_pointrange(aes(ymin = SE1, ymax = SE2), shape = 17, size = 1) +
  dark_theme_classic() +
  scale_color_manual(values=c("#ffb6c1","#43aaa0")) +
  theme(legend.position = "none", 
        text = element_text(size = 30),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

ggsave(Plot, file="plot_AABAInnateMetaSexBetas.pdf", height=5, width =6, useDingbats=FALSE)

Adaptive$AbsBetas = abs(Adaptive$Betas)
Adaptive$SE1 = Adaptive$AbsBeta - Adaptive$StErr
Adaptive$SE2 = Adaptive$AbsBeta + Adaptive$StErr

Adaptive$CellType <- factor(Adaptive$CellType,levels = c("BaboonCD4",
                                                         "BaboonCD8" ,
                                                         "BaboonCD16",         
                                                         "Macaquecd3tcells",
                                                         "Macaquebcells",
                                                         "Macaquecd4tcells",
                                                         "Macaquecd8tcells",
                                                         "TwaKigaCD3",         
                                                          "TwaKigaCD4",
                                                         "TwaKigaCD8",
                                                         "TwaKigaCD20",
                                                         "TurkanaLymphocytes1",
                                                         "TurkanaLymphocytes",
                                                         "NHANESLymphocytes")) 



Plot <- ggplot(Adaptive, aes(x=CellType,y=AbsBetas, color = Lifestyle)) + 
  geom_pointrange(aes(ymin = SE1, ymax = SE2), shape = 17, size = 1) +
  dark_theme_classic() +
  scale_color_manual(values=c("#ffb6c1","#43aaa0")) +
  theme(legend.position = "none", text = element_text(size=40),
                                            axis.ticks.x = element_blank(),
                                            axis.text.x = element_blank()) 

ggsave(Plot, file="plot_AABAAdaptiveMetaSexBetas.pdf", height=6, width =5, useDingbats=FALSE)

