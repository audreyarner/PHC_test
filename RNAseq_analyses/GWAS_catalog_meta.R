library(tidyr)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject/")

#read in entire GWAS catalog database
GWASCatalog <- read.delim("gwas_catalog_v1.0.2-associations_e108_r2023-01-30.tsv", header =T, sep="\t", quote = "")
GWASDisease <- GWASCatalog %>% dplyr::select(DISEASE.TRAIT)
write.table(GWASDisease, "GWASDiseases.txt", quote = F, row.names = F)

#autoimmune diseases
v2 <- c("SjÃ¶gren's syndrome", "Hashimoto thyroiditis", "Primary biliary cirrhosis",
        "Systemic lupus erythematosus", "Graves' disease", "Rheumatoid arthritis", "Multiple sclerosis",
        "Alzheimer's disease", "Celiac disease", "Myasthenia gravis", "Crohn's disease", "Psoriasis",
        "Psoriatic arthritis", "Type 1 diabetes", "Ulcerative colitis", "Ankylosing spondylitis")

#identify genes 
GWASCatalog2 <- GWASCatalog[(GWASCatalog$DISEASE.TRAIT %in% v2),]
GWASCatalog3 <- separate_rows(GWASCatalog2, REPORTED.GENE.S., sep=", ")
GWASCatalog4 <- unique(GWASCatalog3[c("REPORTED.GENE.S.")])

#genes from meta-analysis
df = read.table("MetaAnalysisAllGenes02062023", header=TRUE)

#get ensembl gene IDs for all genes
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=df$Gene, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
res <- inner_join(df, ens2symbol, by=c("Gene"="ENSEMBL"))

#identify autoimmune and non-autoimmune genes
DiseaseGenes <- res[(res$SYMBOL %in% GWASCatalog4$REPORTED.GENE.S.),]
NonDiseaseGenes <- res[(!res$SYMBOL %in% GWASCatalog4$REPORTED.GENE.S.),]

wilcox.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, alternative = "greater")

res$Disease = NA

for(i in 1:nrow(res)){
  if(res$SYMBOL[i] %in% GWASCatalog4$REPORTED.GENE.S.){
    res$Disease[i] = "Auto"
  } else{res$Disease[i] = "No"}
}

#plot of effect sizes
p=ggplot(res, aes(x=Disease,y=Beta,  fill = Disease, color = Disease)) + 
  geom_violin(linewidth = 2) + 
  theme_classic() +
  ylab("Effect size") +
  xlab("Gene category") +
  scale_x_discrete(labels=c("Auto" = "Autoimmune\ndisease genes", "No" = "Non-autoimmune\ndisease genes")) +
  scale_color_manual(values=c("#55bae9", "#dddddd")) +
  scale_fill_manual(values=c("#89cff0", "#f0f0f0")) +
  geom_boxplot(width=0.2, linewidth = 2) + 
  theme(legend.position = "none", text = element_text(size = 30), axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black")) 

p

ggsave(p, file="plot_MetaAutoimmuneCatalogViolin3Oct23.pdf", height=8, width =6, useDingbats=FALSE)

#analysis with difference between genes that are and are not significant and are or are not associated with disease
dat <- data.frame(
  "Signif_no" = c(length(which(DiseaseGenes$pdiff.FDR > 0.1)), length(which(NonDiseaseGenes$pdiff.FDR > 0.1))),
  "Signif_yes" = c(length(which(DiseaseGenes$pdiff.FDR < 0.1)), length(which(NonDiseaseGenes$pdiff.FDR < 0.1))),
  row.names = c("Disease", "Non-Disease"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("NS", "Significant")

dat

#fisher test to see if there is a significant difference
fisher.test(dat)

#analysis also including whether or not there was a larger sex difference in urban vs rural
dat <- data.frame(
  "Signif_no" = c(length(which(DiseaseGenes$pdiff.FDR > 0.1)) + length(which(DiseaseGenes$pdiff.FDR < 0.1 & DiseaseGenes$Beta < 0)), length(which(NonDiseaseGenes$pdiff.FDR < 0.1 & NonDiseaseGenes$Beta < 0)) + length(which(NonDiseaseGenes$pdiff.FDR > 0.1))),
  "Signif_Urban" = c(length(which(DiseaseGenes$pdiff.FDR < 0.1 & DiseaseGenes$Beta >= 0)), length(which(NonDiseaseGenes$pdiff.FDR < 0.1 & NonDiseaseGenes$Beta >= 0))),
  row.names = c("Disease", "Non-Disease"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("NS", "Significant")

dat

fisher.test(dat)


#cancer analysis
v3 <- c("Thyroid cancer", "Colorectal cancer", "Lung cancer", "Multiple myeloma", "Melanoma",
        "Kidney cancer", "Bladder cancer", "Esophageal cancer")

GWASCatalog2 <- GWASCatalog[(GWASCatalog$DISEASE.TRAIT %in% v3),]
GWASCatalog3 <- separate_rows(GWASCatalog2, REPORTED.GENE.S., sep=", ")
GWASCatalog4 <- unique(GWASCatalog3[c("REPORTED.GENE.S.")])

#get ensemble gene ids
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=df$Gene, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
res <- inner_join(df, ens2symbol, by=c("Gene"="ENSEMBL"))

#genes that are and are not associated with cancer
CancerGenes <- res[(res$SYMBOL %in% GWASCatalog4$REPORTED.GENE.S.),]
NonCancerGenes <- res[(!res$SYMBOL %in% GWASCatalog4$REPORTED.GENE.S.),]

wilcox.test(CancerGenes$Beta, NonCancerGenes$Beta, alternative = "less")

res$Disease = NA

for(i in 1:nrow(res)){
  if(res$SYMBOL[i] %in% GWASCatalog4$REPORTED.GENE.S.){
    res$Disease[i] = "Cancer"
  } else{res$Disease[i] = "No"}
}

p=ggplot(res, aes(x=Disease,y=Beta,  fill = Disease, color = Disease)) + 
  geom_violin(linewidth = 2) + 
  theme_classic() +
  ylab("Effect size") +
  xlab("Gene category") +
  scale_x_discrete(labels=c("Cancer" = "Cancer\ndisease genes", "No" = "Non-cancer\ndisease genes")) +
  scale_color_manual(values=c("#f0aa89", "#dddddd")) +
  scale_fill_manual(values=c("#f7d0bd", "#f0f0f0")) +
  geom_boxplot(width=0.2, linewidth = 2) + 
  theme(legend.position = "none", text = element_text(size = 30), axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black")) 

p

ggsave(p, file="plot_CancerCatalogViolin3Oct23.pdf", height=8, width =6, useDingbats=FALSE)

#analysis with genes that are and are not significant in meta analysis and are or are not associated with disease
dat <- data.frame(
  "Signif_no" = c(length(which(CancerGenes$pdiff.FDR > 0.1)), length(which(NonCancerGenes$pdiff.FDR > 0.1))),
  "Signif_yes" = c(length(which(CancerGenes$pdiff.FDR < 0.1)), length(which(NonCancerGenes$pdiff.FDR < 0.1))),
  row.names = c("Disease", "Non-Disease"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("NS", "Significant")

dat

fisher.test(dat)

#analysis also including whether or not there was a larger sex difference in urban vs rural
dat <- data.frame(
  "Signif_no" = c(length(which(CancerGenes$pdiff.FDR > 0.1)) + length(which(CancerGenes$pdiff.FDR < 0.1 & CancerGenes$Beta < 0)), length(which(NonCancerGenes$pdiff.FDR < 0.1 & NonCancerGenes$Beta < 0)) + length(which(NonCancerGenes$pdiff.FDR > 0.1))),
  "Signif_Urban" = c(length(which(CancerGenes$pdiff.FDR < 0.1 & CancerGenes$Beta >= 0)), length(which(NonCancerGenes$pdiff.FDR < 0.1 & NonCancerGenes$Beta >= 0))),
  row.names = c("Cancer", "Non-Cancer"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("NS", "Significant")

dat

fisher.test(dat)