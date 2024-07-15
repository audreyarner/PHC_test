setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")

library(dplyr)
library(qqman)
library(tidyverse)
library(sva)
library(reshape2)
library(data.table)
library(doParallel)
library(edgeR)
library(biomaRt)

#macaque metadata and gene count matrix
meta <- readRDS("metadata_age.rds")
genes1 <- read.delim("counts_matrix_dedupAGATPrimary (1).csv", sep=",")
Parity <- read.delim("parity_females_flow_data_MSR.txt")

meta2019=meta%>%filter(str_detect(trapid, 'T19'))
meta2020=meta%>%filter(str_detect(trapid, 'T20'))
meta2020$year=paste('2020')
meta2019$year=paste('2019')
meta=rbind(meta2019,meta2020)

#clean up parity data
Parity2019=Parity%>%filter(str_detect(blood_sample_date, '2019'))
Parity2020=Parity%>%filter(str_detect(blood_sample_date, '2020'))
Parity2020$year=paste('2020')
Parity2019$year=paste('2019')
Parity=rbind(Parity2019,Parity2020)

colnames(Parity)[1] <- "cayoid"

meta <- merge(meta, Parity, by=c("cayoid", "year"))
#make rownames the genes
rownames(genes1)=genes1$X
genes1 = genes1[,-1]

#select only individuals with lid in metadata
lid <- meta$lid
genes1 <- genes1 %>% dplyr::select(all_of(lid))

abundance <- genes1

abundance$sums=rowSums(abundance)
dim(abundance)
abundance=abundance %>% filter(sums > 0)

genes <- read.delim("genes.txt")
all(genes$converted == genes$name)
genes=genes[,-2]

#genes to remove
remove <- c('AHSP', 'HBA', 'HBE1', 'HBQ1', 'HBM', 'RHEX', 'MB', 'CYGB', 'B2M', 'HPX')
genes=genes %>% filter(!gene %in% remove)
genes<- genes %>% filter (!str_detect(desc, "mitochon"))
genes<- genes %>% filter (!str_detect(desc, "Y-linked"))
genes<- genes %>% filter (!str_detect(gene, "Y_RNA"))
genes<- genes %>% filter (!str_detect(gene, "ENSMMUG0000"))

# optional genes to remove
genes<- genes %>% filter (!str_detect(desc, "spliceoso"))
genes<- genes %>% filter (!str_detect(desc, "ribosom"))
genes<- genes %>% filter (!str_detect(desc, "dehydroge"))


#find any genes that are on the Y chromosome
ensembl=useMart("ensembl", dataset="mmulatta_gene_ensembl")

YAnalysis<-getBM(attributes=c('ensembl_gene_id',"hgnc_symbol",'chromosome_name','start_position','end_position'), mart = ensembl, filters="ensembl_gene_id", values = genes$name)                 

#identify and remove genes on the Y chromosome
NoY <- YAnalysis %>% dplyr::filter(chromosome_name %in% c(1:22) | chromosome_name=='X')

abundance2 <- abundance %>% filter(rownames(abundance) %in% genes$name)
abundance3 <- abundance %>% filter(rownames(abundance) %in% genes$gene)

abundance <- rbind(abundance2, abundance3)

abundance <- distinct(abundance)

all(rownames(abundance) == genes$alias)

# work on original dataset and normalize it
abundance=apply(abundance,2,function(x){x/sum(x)*1e6})
colSums(abundance)
abundance=as.data.frame(abundance)
dim(abundance)

# filter by median expression
gene_sum=apply(abundance,1,sum)
head(gene_sum)
hist(gene_sum)
hist(log10(gene_sum))

# genes with median abundance higher than 7 (10k)
median_abundance<-apply(abundance, 1, median)
hist(median_abundance)
hist(log10(median_abundance))
sum(median_abundance > 7)
abundance10=abundance[median_abundance>7,]
dim(abundance10)

# eliminate all but control
metafilter=meta %>% filter (stimulation %in% c("cont"))
lid=metafilter$lid
abundancefilter <- abundance10 %>% dplyr::select(lid)
all(metafilter$name==colnames(abundancefilter))
dim(abundancefilter)

SequencingDepth <- colSums(abundancefilter)

#add sequencing depth information to covariate table
metafilter$SeqDepth <- SequencingDepth

metafilter <- subset(metafilter, sex == "f")

#functions to mean center data
pretty_up_cols=function(cols){
  cols=refactorize(cols)
  cols=mean_center_clean(cols,"n_infants")
  cols=mean_center_clean(cols,"SeqDepth")
  cols=mean_center_clean(cols,"age")
  
  
  return(cols)
}

refactorize=function(cols_local){
  
  cols_local$sex=factor(cols_local$sex)
  cols_local$year=factor(cols_local$year)
  
  return(cols_local)
}

mean_center_clean=function(cols_local,column){
  index=which(colnames(cols_local)==column)
  
  cols_local[,index]=cols_local[,index]-mean(as.numeric(unlist(cols_local[,index])))
  return(cols_local)
}

control_data = pretty_up_cols(metafilter)

dge <- DGEList(abundancefilter)

#calculate normalization factors for use downstream
dge <- calcNormFactors(dge)

#specify model
design=model.matrix(~n_infants+age+SeqDepth, data = control_data)

#normalize with voom
v <- voom(dge,design, plot = FALSE)

v_combat = ComBat(dat=as.matrix(v$E), batch=control_data$year, mod=design, par.prior=TRUE)

v$E=v_combat

#using duplicate correlations because many individuals were sampled twice
dupcor <- duplicateCorrelation(v,design,block=metafilter$cayoid)

v <- voom(dge,design, plot = FALSE, block=metafilter$cayoid, correlation = dupcor$consensus.correlation)

v_combat = ComBat(dat=as.matrix(v$E), batch=control_data$year, mod=design, par.prior=TRUE)

v$E=v_combat

dupcor <- duplicateCorrelation(v, design, block=metafilter$cayoid)

#fits a linear model using weighted least squares for each gene
fit <-lmFit(v,design, block=metafilter$cayoid, correlation=dupcor$consensus.correlation)

#Empirical Bayes smoothing of standard errors
fit <- eBayes(fit)

#summary of what is significant/insignificant for each covariate
summary(decideTests(fit))

#table of results (with LFC, p-val, adj.P-val, etc for sex)
results=topTable(fit,coef="n_infants", n=Inf)

#get standard errors needed for meta analysis
SEs2 <- as.data.frame(fit$stdev.unscaled[,2]*fit$sigma)

colnames(SEs2)[1] <- "SE"

results <- merge(results, SEs2, by = "row.names")

Betas <- as.data.frame(fit$coefficients[,2])
colnames(Betas)[1] <- "Beta"

# reassigning row names
rownames(results) <- results$Row.names

results2 <- merge(results, Betas, by = "row.names")

sum(results$adj.P.Val < 0.05)
sum(results$adj.P.Val < 0.1)

#look at results
qq(results$P.Value)
hist(results$P.Value)

write.table(results2, "MacaqueParityDEAnalysis02062023.txt", quote=F, row.names = F)

TopGenesInfo <- subset(genes1, rownames(genes1) == "ENSMMUG00000022004")

t <- t(TopGenesInfo)

t <- as.data.frame(t)

t$names <- rownames(t)

Sex <- metafilter %>% dplyr::select(lid, n_infants)
Sex$Lifestyle = "Macaque"
t2 <- merge(Sex, t, by.x = "lid", by.y = "names")

long3 <- melt(setDT(t2), id.vars = c("lid","n_infants", "Lifestyle"), variable.name = "Gene")

long3$value = as.numeric(long3$value)

long3$n_infants = as.numeric(long3$n_infants)

colnames(long3)[2] ="NumberPreg"

p=ggplot(long3, aes(x=NumberPreg, y=log10(long3$value))) +
  dark_theme_bw() +
  geom_point(aes(color = Lifestyle), size = 3) +
  geom_smooth(method = lm, se = F) + 
  scale_color_manual(values=c("#9595ea")) +
  theme(text = element_text(size = 30), legend.position = "none")
p





