setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")
library(metafor)

#read in each of the populations statistics
macaque <- read.delim("MacaqueParityCellModelStatistics.txt", sep=" ")
macaque$Immunity = c("Adaptive", "Adaptive", "Innate", "Adaptive", "Adaptive", "Adaptive", "Adaptive", "Innate", "Innate", "Innate")

NHANES <- read.delim("NHANESParityCellModelStatistics.txt", sep = " ")
NHANES$Immunity = c("Innate", "Innate", "Innate", "Innate", "Adaptive")

Turkana = read.delim("TurkanaParityModelStatistics.txt", sep = " ")
Turkana$Immunity = c("Innate", "Innate", "Innate", "Innate", "Adaptive")

AllData = rbind(macaque, NHANES, Turkana)

#run meta-analysis models
m.qual <- rma(yi = Betas,
              sei = StErr,
              data = AllData,
              method = "ML",
              mods = ~ 1)
m.qual

Adaptive = subset(AllData, Immunity == "Adaptive")
m.qual3 <- rma(yi = Betas,
               sei = StErr,
               data = Adaptive,
               method = "ML",
               mods = ~ 1)
m.qual3

Innate = subset(AllData, Immunity == "Innate")
m.qual4 <- rma(yi = Betas,
               sei = StErr,
               data = Innate,
               method = "ML",
               mods = ~ 1)
m.qual4

#graph results
AllData$CellType = row.names(AllData)
AllData$SE1 = AllData$Betas - AllData$StErr
AllData$SE2 = AllData$Betas + AllData$StErr

AllData$CellType <- factor(AllData$CellType,levels = c("Macaquenkcells",
                                                     "Macaquecd8tregs",              
                                                     "Macaquetregs",
                                                     "Macaquecd16hladrpos",          
                                                     "Macaquecd14cd16hladrpos",
                                                     "Macaquecd14hladrpos",          
                                                     "Macaquecd3tcells",
                                                     "Macaquebcells",
                                                     "Macaquecd4tcells",
                                                     "Macaquecd8tcells",
                                                     "TurkanaBasophils1",
                                                     "TurkanaEosinophils1",          
                                                     "TurkanaSegmented.neutrophils1",
                                                     "TurkanaMonocytes1",
                                                     "TurkanaLymphocytes1",
                                                     "TurkanaLymphocytes",
                                                     "TurkanaBasophils",
                                                     "TurkanaEosinophils",          
                                                     "TurkanaSegmented.neutrophils",
                                                     "TurkanaMonocytes",
                                                     "NHANESBasophils",
                                                     "NHANESEosinophils",
                                                     "NHANESSegmented.neutrophils",
                                                     "NHANESMonocytes",
                                                     "NHANESLymphocytes")) 


Plot <- ggplot(AllData, aes(x=CellType,y=Betas)) +
  geom_hline(yintercept=0,linetype=2, linewidth = 2, color = "#808080") +
  geom_pointrange(aes(ymin = SE1, ymax = SE2), size = 2, colour = "#ec9166") +
  dark_theme_classic() +
  theme(legend.position = "none", 
        text = element_text(size = 30),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

ggsave(Plot, file="plot_AABAParityMetaBetas.pdf", height=6, width =10, useDingbats=FALSE)

