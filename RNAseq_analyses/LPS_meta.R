setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")
library(dplyr)

#lps genes
TwaKiga <- read.csv("TwaKigaLPSGenes.csv", header = T)
Baboon <- read.csv("BaboonLPSGenes.csv", header =T)
Macaca <- read.csv("MacaqueLPSGenes.csv", header = T)

TwaBaboon <- inner_join(TwaKiga, Baboon, by=c("X"="Gene"))

All <- inner_join(TwaBaboon, Macaca, by=c("X"="gene"))

#fdr correction
LPSFDR <- unname(unlist(All[,"Condition.effect.in.both.sexes..p.value"]))
fdr <- p.adjust(LPSFDR, method = "fdr")
All$Baboonpdiff.FDR <- fdr

#only consider significant if FDR < 0.1 in all populations
Significant <- All %>% filter(FDR.in.LPS. < .1 & PopDE_LPS_FDR < .1 & Baboonpdiff.FDR < .1)
NonSignificant <- All %>% filter(FDR.in.LPS. >= .1 & PopDE_LPS_FDR >= .1 & Baboonpdiff.FDR >= .1)

df = read.table("MetaAnalysisAllGenes02062023", header=TRUE)

#get ensembl gene name for all genes
library(org.Hs.eg.db)

ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=df$Gene, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
res <- inner_join(df, ens2symbol, by=c("Gene"="ENSEMBL"))

#identify genes from the meta-analysis that were LPS DE or not
LPSGenes <- res[(res$SYMBOL %in% Significant$X),]
NonLPSGenes <- res[(res$SYMBOL %in% NonSignificant$X),]

wilcox.test(LPSGenes$Beta, NonLPSGenes$Beta, alternative = "two.sided")

res$Disease = NA

for(i in 1:nrow(res)){
  if(res$SYMBOL[i] %in% Significant$X){
    res$Disease[i] = "LPS"
  } else{res$Disease[i] = "No"}
}

length(which(res$Disease=="LPS"))
length(which(res$Disease=="No"))

p=ggplot(res, aes(x=Disease,y=Beta,  fill = Disease, color = Disease)) + 
  geom_violin(linewidth = 2) + 
  dark_theme_classic() +
  scale_color_manual(values=c("#55bae9", "#e98455")) +
  scale_fill_manual(values=c("#89cff0", "#f0aa89")) +
  geom_boxplot(width=0.2, linewidth = 2) + 
  theme(legend.position = "none", text = element_text(size = 30)) 

ggsave(p, file="plot_MetaLPSCatalogViolin.pdf", height=6, width =6, useDingbats=FALSE)


dat <- data.frame(
  "Signif_no" = c(length(which(LPSGenes$pdiff.FDR > 0.1)), length(which(NonLPSGenes$pdiff.FDR > 0.1))),
  "Signif_yes" = c(length(which(LPSGenes$pdiff.FDR < 0.1)), length(which(NonLPSGenes$pdiff.FDR < 0.1))),
  row.names = c("LPS", "Non-LPS"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("NS", "Significant")

dat

fisher.test(dat)

dat <- data.frame(
  "Signif_no" = c(length(which(LPSGenes$pdiff.FDR > 0.1)) + length(which(LPSGenes$pdiff.FDR < 0.1 & LPSGenes$Beta < 0)), length(which(NonLPSGenes$pdiff.FDR < 0.1 & NonLPSGenes$Beta < 0)) + length(which(NonLPSGenes$pdiff.FDR > 0.1))),
  "Signif_Urban" = c(length(which(LPSGenes$pdiff.FDR < 0.1 & LPSGenes$Beta >= 0)), length(which(NonLPSGenes$pdiff.FDR < 0.1 & NonLPSGenes$Beta >= 0))),
  row.names = c("LPS", "Non-LPS"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("NS", "Significant")

dat

fisher.test(dat)
