setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")
library(edgeR)
library("readxl")
library(sva)
library(dplyr)
library(biomaRt)
library(qqman)
library(ggplot2)

#Harrison data with covariate information
cols <- read.delim("Inputs/1_DE_analyses/RNAseq_metadata.txt", sep = " ")

### Get condition-specific metacols tables
cols_CTL=cols[which(cols$Condition=="CTL"),]

#mean center covariates and make sure others are factors with following three functions
pretty_up_cols=function(cols){
  
  cols=refactorize(cols)
  
  cols=mean_center_clean(cols,"CD3")
  cols=mean_center_clean(cols,"CD4")
  cols=mean_center_clean(cols,"CD8")
  cols=mean_center_clean(cols,"CD14")
  cols=mean_center_clean(cols,"CD20")
  cols=mean_center_clean(cols,"CD56")
  cols=mean_center_clean(cols,"fraction_assigned")
  
  return(cols)
}

mean_center_clean=function(cols_local,column){
  index=which(colnames(cols_local)==column)
  
  cols_local[,index]=cols_local[,index]-mean(as.numeric(unlist(cols_local[,index])))
  return(cols_local)
}

refactorize=function(cols_local){
  
  cols_local$Genotyping_ID=factor(cols_local$Genotyping_ID,levels=unique(as.character(cols_local$Genotyping_ID)))
  cols_local$Condition=factor(cols_local$Condition)
  cols_local$Flowcell=factor(cols_local$Flowcell)
  cols_local$Sex=factor(cols_local$Sex)
  
  return(cols_local)
}

reads <- read.delim("Inputs/1_DE_analyses/RNAseq_reads_matrix.txt", sep = " ")

cols=cols[which(cols$PopDE_set==1),]
rownames(cols)=cols$Sample
cols=cols[order(rownames(cols)),]
reads=reads[,which(colnames(reads) %in% rownames(cols))]
reads=reads[,order(colnames(reads))]

length(which(rownames(cols)!=colnames(reads)))
# 0

### Get condition-specific reads matrices
reads_CTL=reads[,which(cols$Condition=="CTL")]

### Get condition-specific metadata tables
cols_CTL=cols[which(cols$Condition=="CTL"),]

### Drop absent levels & Mean-center technical covariates in condition-specific metadata
cols_CTL=pretty_up_cols(cols_CTL)

#drop Y chrom genes
ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")

YAnalysis<-getBM(attributes=c('ensembl_gene_id',"hgnc_symbol",'chromosome_name','start_position','end_position'), mart = ensembl, filters="ensembl_gene_id", values = rownames(reads_CTL))                 

#identify and remove genes on the Y chromosome
NoY1 <- YAnalysis %>% dplyr::filter(chromosome_name %in% c(1:22) | chromosome_name=='X')

YAnalysis<-getBM(attributes=c('ensembl_gene_id',"hgnc_symbol",'chromosome_name','start_position','end_position'), mart = ensembl, filters="hgnc_symbol", values = rownames(reads_CTL))                 

#identify and remove genes on the Y chromosome
NoY2 <- YAnalysis %>% dplyr::filter(chromosome_name %in% c(1:22) | chromosome_name=='X')

NoY <- rbind(NoY1, NoY2)

reads_CTL1 <- subset(reads_CTL, rownames(reads_CTL) %in% NoY$ensembl_gene_id)
reads_CTL2 <- subset(reads_CTL, rownames(reads_CTL) %in% NoY$hgnc_symbol)

reads_CTL <- rbind(reads_CTL1, reads_CTL2)

#Form a DGEList object combining the counts and associated annotation
dge <- DGEList(reads_CTL)

#calculate cpm
dge <- calcNormFactors(dge)

cpm <- cpm(dge, log = TRUE, normalized.lib.sizes=TRUE)

cpm1 <- as.data.frame(cpm)

cpm1$Median <-apply(cpm1,1,median)

cpm2 <- cpm1 %>% filter(Median>=.1)

reads_CTL <- subset(reads_CTL, rownames(reads_CTL) %in% rownames(cpm2))

dge <- DGEList(reads_CTL)

#calculate normalization factors for use downstream
dge <- calcNormFactors(dge)

#specify model
design=model.matrix(~Sex+fraction_assigned+Admixture+CD14+CD4+CD20, data = cols_CTL)

#normalize with voom
v <- voom(dge,design, plot = F)

#control for sequencing flow cell batch effects
v_combat = ComBat(dat=as.matrix(v$E), batch=cols_CTL$Flowcell, mod=design, par.prior=TRUE)

v$E=v_combat

#fits a linear model using weighted least squares for each gene
fit <-lmFit(v,design)

#Empirical Bayes smoothing of standard errors
fit <- eBayes(fit)

#summary of what is significant/insignificant for each covariate
summary(decideTests(fit))

#table of results (with LFC, p-val, adj.P-val, etc for sex)
results=topTable(fit,coef="SexMale", n=Inf)

#get standard errors
SEs2 <- as.data.frame(fit$stdev.unscaled[,2]*fit$sigma)
colnames(SEs2)[1] <- "SE"

results <- merge(results, SEs2, by = "row.names")

Betas <- as.data.frame(fit$coefficients[,2])
colnames(Betas)[1] <- "Beta"

# reassigning row names
rownames(results) <- results$Row.names

results2 <- merge(results, Betas, by = "row.names")

sum(results$adj.P.Val < 0.05)

qq(results$P.Value)

hist(results$P.Value)

write.table(results2, "TwaKigaRawDEAnalysis1-31-22.txt", quote=F, row.names = F)

#sex permutation
GeomPolyTable <- data.frame(Pvals=numeric(),
                            Permutation=character(), 
                            stringsAsFactors=FALSE) 

n=1

#from here on permutation
for(j in 1:20){
  cols_CTL$SexPermute <- sample(cols_CTL$Sex, nrow(cols_CTL), replace = FALSE, prob = NULL)
  dge <- DGEList(reads_CTL)
  
  #calculate normalization factors for use downstream
  dge <- calcNormFactors(dge)
  
  #specify model
  design=model.matrix(~SexPermute+fraction_assigned+Admixture+CD14+CD4+CD20, data = cols_CTL)
  
  #normalize with voom
  v <- voom(dge,design, plot = F)
  
  #control for sequencing flow cell batch effects
  v_combat = ComBat(dat=as.matrix(v$E), batch=cols_CTL$Flowcell, mod=design, par.prior=TRUE)
  
  v$E=v_combat
  
  #fits a linear model using weighted least squares for each gene
  fit <-lmFit(v,design)
  
  #Empirical Bayes smoothing of standard errors
  fit <- eBayes(fit)
  
  results=topTable(fit,coef="SexPermuteMale", n=Inf)
  
  PTable <- results %>% dplyr::select(P.Value) 
  
  PTable <- mutate(PTable, Permutation=n)
  
  GeomPolyTable <- rbind(GeomPolyTable,PTable)
  
  n=n+1
  
}


GeomPolyTable$Permutation <- as.character(GeomPolyTable$Permutation)
FigOverall <- ggplot(GeomPolyTable, aes(P.Value)) +
  geom_freqpoly(bins=100)

write.table(GeomPolyTable, "TwaKigaSexLMempiricalPerm.txt", quote=F, sep = "\t")
