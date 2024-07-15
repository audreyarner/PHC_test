#using fgsea package
library(data.table)
library(fgsea)
library(ggplot2)
library(org.Hs.eg.db)
library(tibble)

df <- read.table("MetaAnalysisAllGenes02062023", header = T)

#get gene symbols from ensembl ID
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=df$Gene, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
ens2symbol

res <- inner_join(df, ens2symbol, by=c("Gene"="ENSEMBL"))

res2 <- res %>% arrange(Beta) %>%
  dplyr::select(SYMBOL, Beta)

ranks <- deframe(res2)

#use gene ontology gene sets
pathways.hallmark <- gmtPathways("c5.go.v2022.1.Hs.symbols.gmt.txt")

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks) %>% 
  as_tibble() %>% 
  arrange(padj)

#remove list before writing table
fgseaRes <- fgseaRes[-c(8)]

write.table(fgseaRes, "GSEAAllGenes.txt", quote=F, row.names = F)

#use clusterProfiler
library(clusterProfiler)

sample_gene = res %>% filter(pdiff.FDR <= 0.1) %>% 
  dplyr::select(Gene) 

sample_gene2 = sample_gene[,1]

sample_test <- enrichGO(sample_gene2, OrgDb=org.Hs.eg.db, ont = "BP", keyType = "ENSEMBL", pvalueCutoff=1, qvalueCutoff=1)
Ex = as.data.frame(head((sample_test), 120))

sample_test <- enrichGO(sample_gene2, OrgDb=org.Hs.eg.db, ont = "MF", keyType = "ENSEMBL", pvalueCutoff=1, qvalueCutoff=1)
Ex = as.data.frame(head((sample_test), 120))

sample_gene = res %>% filter(pdiff.FDR <= 0.1 & Beta < 0) %>% 
  dplyr::select(Gene) 

sample_gene2 = sample_gene[,1]

sample_test <- enrichGO(sample_gene2, OrgDb=org.Hs.eg.db, ont = "BP", keyType = "ENSEMBL", pvalueCutoff=1, qvalueCutoff=1)
head(summary(sample_test))

sample_test <- enrichGO(sample_gene2, OrgDb=org.Hs.eg.db, ont = "MF", keyType = "ENSEMBL", pvalueCutoff=1, qvalueCutoff=1)
head(summary(sample_test))

sample_gene = res %>% filter(pdiff.FDR <= 0.1 & Beta >= 0) %>% 
  dplyr::select(Gene) 

sample_gene2 = sample_gene[,1]

sample_test <- enrichGO(sample_gene2, OrgDb=org.Hs.eg.db, ont = "BP", keyType = "ENSEMBL", pvalueCutoff=1, qvalueCutoff=1)
head(summary(sample_test))

sample_test <- enrichGO(sample_gene2, OrgDb=org.Hs.eg.db, ont = "MF", keyType = "ENSEMBL", pvalueCutoff=1, qvalueCutoff=1)
head(summary(sample_test))


