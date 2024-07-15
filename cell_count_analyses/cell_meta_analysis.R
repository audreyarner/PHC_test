setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")

library(dplyr)
library(biomaRt)
library(metafor)

baboondata <- read.table("BaboonLMEOutput2-1-2023.txt", header =T)

#remove genes whose variance is less than 0, get SE
baboondata <- baboondata %>% filter(Together2.sexM >= 0)
baboondata <- mutate(baboondata, SE = sqrt(baboondata$Together2.sexM)/sqrt(98))

#get baboon homologs and clean data, make ENSEMBL ID the row names
BaboonHomologs <- read.table("mart_BaboonOrthologs.txt", fill=T, skip = 1, na.strings=c("","NA"))
BaboonHomologs2 <- na.omit(BaboonHomologs)

BaboonHomologs3 <- BaboonHomologs2[!duplicated(BaboonHomologs2$V3), ]
colnames(BaboonHomologs3)[2] <- "genename"
colnames(BaboonHomologs3)[3] <- "ensembl_gene_id"

BaboonHomologsSmall <- BaboonHomologs3 %>% dplyr::select(ensembl_gene_id, genename)

BaboonLargeData <- baboondata
BaboonLargeData$genename <- row.names(BaboonLargeData)
BaboonLargeData <- merge(BaboonLargeData, BaboonHomologsSmall)
row.names(BaboonLargeData) <- BaboonLargeData$ensembl_gene_id
colnames(BaboonLargeData)[3] <- "beta"
colnames(BaboonLargeData)[15] <- "pval"

#list of baboon ensembl IDs
BaboonEnsg <- data.frame(row.names(BaboonLargeData))

BaboonLargeData <- mutate(BaboonLargeData, Lifestyle = "Traditional")
BaboonLargeData <- mutate(BaboonLargeData, Population = "Baboon")

BaboonLargeData <- BaboonLargeData %>% dplyr::select(beta, SE, Lifestyle, Population)

colnames(BaboonEnsg)[1] <- "ensembl_gene_id"

#macaque data read in
macaquedata <- read.table("MacaqueRawDEAnalysis2-2.txt", header = T)

#Macaque homologs
Macaca1 <- read.table("mart_macacaEnsg.txt", fill = T, skip = 1, na.strings = c("", "NA"))
Macaca2 <- read.table("mart_macacaName.txt", fill = T, skip = 1, na.strings = c("", "NA"))

Macacahomologs <- rbind(Macaca1, Macaca2)
MacacaHomologs2 <- na.omit(Macacahomologs)
MacacaHomologs3 <- MacacaHomologs2[!duplicated(MacacaHomologs2$V3), ]
MacacaMerge <- MacacaHomologs3 %>% dplyr::select(V3, V2,V1)
MacacaData1 <- merge(MacacaMerge, macaquedata, by.x="V1", by.y="Row.names")
MacacaData2 <- merge(MacacaMerge, macaquedata, by.x="V2", by.y="Row.names")
MacacaData3 <- rbind(MacacaData1, MacacaData2)

colnames(MacacaData3)[12] <- "beta"

MacacaData3 <- mutate(MacacaData3, Lifestyle = "Traditional")
MacacaData3 <- mutate(MacacaData3, Population = "Macaque")

#remove duplicates & make ensembl ids the row names
MacacaData3 <- MacacaData3[!duplicated(MacacaData3$V3), ]
row.names(MacacaData3) <- MacacaData3$V3

MacacaData4 <- MacacaData3 %>% dplyr::select(beta, SE, Lifestyle, Population)
MacacaEnsg <- MacacaHomologs3 %>% dplyr::select(V3)

colnames(MacacaEnsg)[1] <- "ensembl_gene_id"

#Batwa/BakigaData convert to ENSG
TwaKiga <- read.table("TwaKigaRawDEAnalysis1-31-22.txt", header = T)

#get ensembl ids
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
TwaKigaEnsg<-getBM(attributes=c('ensembl_gene_id',"hgnc_symbol"), mart = ensembl, filters="hgnc_symbol", values =TwaKiga$Row.names)                 

TwaKiga2 <- TwaKiga[TwaKiga$Row.names %in% TwaKigaEnsg$hgnc_symbol, ]
TwaKigaEnsg <- TwaKigaEnsg[!duplicated(TwaKigaEnsg$hgnc_symbol), ]

row.names(TwaKigaEnsg) <- TwaKigaEnsg$ensembl_gene_id

colnames(TwaKiga2)[10] <- "beta"

TwaKiga2 <- mutate(TwaKiga2, Lifestyle = "Traditional")
TwaKiga2 <- mutate(TwaKiga2, Population = "TwaKiga")

row.names(TwaKiga2) <- TwaKigaEnsg$ensembl_gene_id

TwaKiga22 <- TwaKiga2 %>% dplyr::select(beta, SE, Lifestyle, Population)


#GTEx data pre-calculated effect sizes and betas
GTExBeta <- read.table("effect_size.tsv", header = T)
GTExBeta <- na.omit(GTExBeta)
GTExWB <- GTExBeta %>% dplyr::select(Whole_Blood)
colnames(GTExWB)[1] <- "beta"

GTExSE <- read.table("effect_size_se.tsv", header =T)
GTExSE <- na.omit(GTExSE)
GTExSE <- GTExSE %>% dplyr::select(Whole_Blood)
colnames(GTExSE)[1] <- "SE"
GTEx <- cbind(GTExWB, GTExSE)

GTEx <- mutate(GTEx, Lifestyle = "Novel")
GTEx <- mutate(GTEx, Population = "GTEX")

GeneName = function(information){
  return(gsub("\\..*", "", information))
}

#get smaller ensembl id so can merge
rownames(GTEx) <- GeneName((rownames(GTEx)))
GTEx <- mutate(GTEx, ensembl_gene_id = row.names(GTEx))
GTEx2 <- GTEx %>% dplyr::select(beta, SE, Lifestyle, Population)
GTExEnsg <- data.frame(GTEx$ensembl_gene_id)

colnames(GTExEnsg)[1] <- "ensembl_gene_id"

#Turkana Pastoralist Data
TurkanaPastoral=read.delim('TurkanaPastoralistLMEOutput.txt')

#get SE and clean data
TurkanaPastoral <- mutate(TurkanaPastoral, SE = sqrt(TurkanaPastoral$MergeAll2.GenderM)/sqrt(43))
colnames(TurkanaPastoral)[3] <- "beta"

TurkanaPastoral <- mutate(TurkanaPastoral, Lifestyle = "Traditional")
TurkanaPastoral <- mutate(TurkanaPastoral, Population = "PastoralTurkana")


TurkanaPastoral2 <- TurkanaPastoral %>% dplyr::select(beta, SE, Lifestyle, Population)

TurkanaPastoralEnsg <- data.frame(row.names(TurkanaPastoral))

colnames(TurkanaPastoralEnsg)[1] <- "ensembl_gene_id"

#urban Turkana data cleaning
TurkanaUrban=read.delim('TurkanaUrbanLMEOutput.txt')

TurkanaUrban <- mutate(TurkanaUrban, SE = sqrt(TurkanaUrban$MergeAll2.GenderM)/sqrt(201))
colnames(TurkanaUrban)[3] <- "beta"

TurkanaUrban <- mutate(TurkanaUrban, Lifestyle = "Novel")
TurkanaUrban <- mutate(TurkanaUrban,Population = "UrbanTurkana")

TurkanaUrban2 <- TurkanaUrban %>% dplyr::select(beta, SE, Lifestyle, Population)
TurkanaUrbanEnsg <- data.frame(row.names(TurkanaUrban))

colnames(TurkanaUrbanEnsg)[1] <- "ensembl_gene_id"

#merge by ensembl gene id
TurkanaEnsg <- merge(TurkanaPastoralEnsg, TurkanaUrbanEnsg)

MergeBaboon <- merge(TurkanaEnsg, BaboonEnsg)

MergeMacaca <- merge(MacacaEnsg, MergeBaboon)

MergeGTEx <- merge(MergeMacaca, GTExEnsg)

MergeAll <- merge(MergeGTEx, TwaKigaEnsg)

AllGenes <- as.character(MergeAll[["ensembl_gene_id"]])

populations <- list(BaboonLargeData)


GeneNameDF <- data.frame(Gene = character(),
                        Pvalue = numeric(),
                        Beta = numeric(),
                        stringsAsFactors=FALSE) 

#meta-analysis for each gene as a permutation
for(i in 1:length(AllGenes)){
  GenePopDF <- data.frame(beta = numeric(),
                          SE = numeric(), 
                          Lifestyle = character(),
                          stringsAsFactors=FALSE) 
  GeneName<-AllGenes[i]
  TurkanaTrad <- TurkanaPastoral2[GeneName,]
  TurkanaUrb <- TurkanaUrban2[GeneName,]
  BaboonTrad <- BaboonLargeData[GeneName,]
  MacacaTrad <- MacacaData4[GeneName,]
  GTeXNovel <- GTEx2[GeneName,]
  TwaKigaTrad <- TwaKiga22[GeneName,]
  New <- rbind(GenePopDF, TurkanaTrad, BaboonTrad, TurkanaUrb, MacacaTrad, GTeXNovel, TwaKigaTrad)
  m.qual1 <- rma(yi = beta,
                   sei = SE,
                   data = New,
                   method = "ML",
                   mods = ~ Lifestyle,
                  control=list(maxiter=1000))
  pval <- m.qual1[["pval"]][2]
  beta <- m.qual1[["beta"]][2]
  
  
  Gene <- c(GeneName)
  Pvalue <- c(pval)
  Beta <- c(beta)
  
  df <- data.frame(Gene, Pvalue, Beta)

  GeneNameDF <- rbind(GeneNameDF, df)
}

LifestyleSexFDR <- unname(unlist(GeneNameDF[,"Pvalue"]))
fdr <- p.adjust(LifestyleSexFDR, method = "fdr")
GeneNameDF$pdiff.FDR <- fdr

#how many genes fall into each category?
length(which(GeneNameDF$pdiff.FDR<0.1))

length(which(GeneNameDF$pdiff.FDR<0.1 & GeneNameDF$Beta <= 0))
length(which(GeneNameDF$pdiff.FDR<0.1 & GeneNameDF$Beta > 0))

#make dataframes with different categories of genes
#higher sex difference in novel environment
NovelSignificant = GeneNameDF %>% filter( pdiff.FDR<0.1 & GeneNameDF$Beta <= 0)
#higher sex difference in traditional environment
TraditionalSignificant = GeneNameDF %>% filter( pdiff.FDR<0.1 & GeneNameDF$Beta > 0)

NovelSignificant2 <- NovelSignificant %>% dplyr::select(Gene)
TraditionalSignificant2 <- TraditionalSignificant %>% dplyr::select(Gene)

write.table(NovelSignificant2, "MetaAnalysisNovelGenes02062023", quote = F, row.names = F)
write.table(TraditionalSignificant2, "MetaAnalysisATraditionalGenes02062023", quote = F, row.names = F)

write.table(GeneNameDF, "MetaAnalysisAllGenes02062023", quote = F, row.names = F)

#plot betas for significant genes
SignificantGenes <- GeneNameDF[order(GeneNameDF$pdiff.FDR),]

SignificantGenes <- head(SignificantGenes, 10)

GenePopDF <- data.frame(beta = numeric(),
                        SE = numeric(), 
                        Lifestyle = character(),
                        Population = character(),
                        Gene = character(),
                        stringsAsFactors=FALSE) 

for(i in 1:nrow(SignificantGenes)){
  GeneName<-SignificantGenes$Gene[i]
  TurkanaTrad <- TurkanaPastoral2[GeneName,]
  TurkanaUrb <- TurkanaUrban2[GeneName,]
  BaboonTrad <- BaboonLargeData[GeneName,]
  MacacaTrad <- MacacaData4[GeneName,]
  GTeXNovel <- GTEx2[GeneName,]
  TwaKigaTrad <- TwaKiga22[GeneName,]
  New <- rbind(TurkanaTrad, BaboonTrad, TurkanaUrb, MacacaTrad, GTeXNovel, TwaKigaTrad)
  New <- mutate(New, Gene = GeneName)
  
  GenePopDF <- rbind(GenePopDF, New)
}

write.table(GenePopDF, "SignificantMetaAllPopBetas.txt", quote=F, row.names=F)

library(ggplot2)

P = ggplot(GenePopDF, aes(x=Gene, y=beta, color=Lifestyle)) +
  geom_hline(yintercept=0,linetype=2, size = 2) +
  geom_pointrange((aes(ymin=beta-SE, ymax=beta+SE)), shape = 17, size = 1.5) +
  dark_theme_classic() +
  scale_color_manual(values=c("#ffb6c1","#43aaa0")) +
  theme(legend.position = "none", axis.text.x=element_text(angle = -90, hjust = 0))

ggsave(P, file="plot_MetaPCHTop10.pdf", height=6, width =8, useDingbats=FALSE)


#background genes for use in further analyses
BackgroundGenes <- GeneNameDF %>% filter(pdiff.FDR > 0.1)

NovelBackground = GeneNameDF %>% filter( pdiff.FDR>0.1 & GeneNameDF$Beta >= 0)
TraditionalBackground = GeneNameDF %>% filter( pdiff.FDR>0.1 & GeneNameDF$Beta < 0)

NovelBackground2 <- NovelBackground %>% dplyr::select(Gene)
write.table(NovelBackground2, "NovelBackgroundMetaAnalysis.txt", quote=F, row.names = F)

TraditionalBackground2 <- TraditionalBackground %>% dplyr::select(Gene)
write.table(TraditionalBackground2, "TraditionalBackgroundMetaAnalysis.txt", quote=F, row.names = F)

BackgroundGenes2 <- BackgroundGenes %>% dplyr::select(Gene)
write.table(BackgroundGenes2, "BackgroundGenesMetaAnalysis.txt", quote=F, row.names = F)

FocalGenes <- GeneNameDF %>% filter(pdiff.FDR <= 0.1)

FocalGenes2 <- FocalGenes %>% dplyr::select(Gene)
write.table(FocalGenes2, "FocalGenesMetaAnalysis.txt", quote=F, row.names = F)




