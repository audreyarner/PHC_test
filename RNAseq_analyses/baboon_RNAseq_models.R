setwd("~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/RotationProject")

library(edgeR)
library("readxl")
library(dplyr)
library(ggplot2)
library(plyr)
library(EMMREML)
library(doParallel)
library(qqman)
library(reshape2)
library(data.table)
library(biomaRt)

#Baboon data with covariate information
control_data <- read.delim("BaboonCovariates.txt", header=T)

#Take only those individs who didn't have a treatment
control_data <- control_data %>% filter(treatment=="NULL")

#baboon count matrix & protein coding genes
GeneMatrix <- read.delim("BaboonRNASeqcounts.txt", header =T)
GTFFile <- read.delim("Panubis1.0ProteinCoding.txt", header=F)

#functions to get gene names from GTF files
MiddleName = function(information){
  return(gsub(".*gene_id ", "", information))
}

GeneName = function(information){
  return(gsub(";.*","", information))
}

#extract gene names
ProteinCoding <- mutate(GTFFile, GeneID=GeneName(MiddleName(GTFFile$V9)))

#only use genes that are protein coding
CodingMatrix <-GeneMatrix[(row.names(GeneMatrix) %in% ProteinCoding$GeneID),]

#find any genes that are on the Y chromosome
ensembl=useMart("ensembl", dataset="panubis_gene_ensembl")

YAnalysis<-getBM(attributes=c('ensembl_gene_id',"hgnc_symbol",'chromosome_name','start_position','end_position'), mart = ensembl, filters="hgnc_symbol", values = row.names(CodingMatrix))                 

#identify and remove genes on the Y chromosome
NoY <- YAnalysis %>% dplyr::filter(chromosome_name %in% c(1:22) | chromosome_name=='X')

CodingMatrix <- subset(CodingMatrix, rownames(CodingMatrix) %in% NoY$hgnc_symbol)

#only use the null treatment samples
CodingMatrix2 <- CodingMatrix[, (colnames(CodingMatrix) %in% control_data$sra_sname), drop = FALSE]

#only use genes with our threshold cpm
dge <- DGEList(CodingMatrix2)

dge <- calcNormFactors(dge)

cpm <- cpm(dge, log = TRUE, normalized.lib.sizes=TRUE)

cpm1 <- as.data.frame(cpm)

cpm1$Median <-apply(cpm1,1,median)

cpm2 <- cpm1 %>% filter(Median>=.1)

#subset only those genes that hit our cpm threshold
OnlyUsed <- subset(CodingMatrix2, rownames(CodingMatrix2) %in% rownames(cpm2))

#find sequencing depth (sum of each column) for each sample
SequencingDepth <- colSums(CodingMatrix2)

#add sequencing depth information to covariate table
control_data$SeqDepth <- SequencingDepth

#functions for normalizing data
pretty_up_cols=function(cols){
  
  cols=refactorize(cols)
  
  cols=mean_center_clean(cols,"age")
  cols=mean_center_clean(cols,"SeqDepth")
  cols=mean_center_clean(cols,"Flow_PC1")
  cols=mean_center_clean(cols,"Flow_PC2")
  cols=mean_center_clean(cols,"Flow_PC3")
  
  return(cols)
}

mean_center_clean=function(cols_local,column){
  index=which(colnames(cols_local)==column)
  
  cols_local[,index]=cols_local[,index]-mean(as.numeric(unlist(cols_local[,index])))
  return(cols_local)
}

refactorize=function(cols_local){
  
  cols_local$sex=factor(cols_local$sex)
  cols_local$year=factor(cols_local$year)
  
  return(cols_local)
}

#normalize relevant factors
control_data2 = pretty_up_cols(control_data)

#get names listed in GRM from SI tables
Names <- read.delim("FINAL_SI_tables.txt", header = F)

Together <- merge(control_data2, Names, by.x = "sra_sname", by.y="V2")

#kinship matrix
GRM <- read.delim("13Jun17_kmatrix.txt")

#take out the duplicated sample names (there are two of each for the control vs treatment)
GRM2 <- unique(GRM)

GRM2 <- GRM2[!duplicated(as.list(GRM2))]

#make the rownames the first col list of individual IDs, then remove that first row
rownames(GRM2) <- GRM2[,1]
GRM2[,1] <- NULL

#only use individuals who are used in the RNA-seq
GRM2 <- GRM2[,(row.names(GRM2) %in% Together$V1)]
GRM2 <- GRM2[(row.names(GRM2) %in% Together$V1),]

#Only have RNASeq data from individs used
OnlyUsed <- OnlyUsed[, (colnames(OnlyUsed) %in% Together$sra_sname), drop = FALSE]

#order the columns of gene count matrix, other necessary matrices already ordered
Ordering <- colnames(OnlyUsed)
Together2 <- Together[match(Ordering, Together$sra_sname),]

#make z matrix for EMMA
mat=matrix(nrow=nrow(Together2),ncol=length(unique(Together2[,"sra_sname"]))); colnames(mat)=unique(Together2[,"sra_sname"]); rownames(mat)=Together2[,"sra_sname"]
for (r in 1:nrow(mat)) {
  for (c in 1:ncol(mat)) {
    if (rownames(mat)[r]==colnames(mat)[c]) {mat[r,c]=1}else {mat[r,c]=0}} }

#make year character, not number
Together2$year <- as.character(Together2$year)

#get residuals using limma
design <- model.matrix(~0+year+SeqDepth, data = Together2)
dge <- DGEList(OnlyUsed)
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=FALSE)
fit <-lmFit(v,design)
resids <- residuals.MArrayLM(object=fit, v)

#mixed effect modeling with emmreml (following Jordan & Noah's code)
ncores <- detectCores(logical = TRUE)  # find number of available cores
clus <- makeCluster(ncores, setup_strategy = "sequential")  # create cluster
registerDoParallel(cores = ncores)
clusterExport(clus,
              varlist = c("resids", "Together2", "mat", "GRM2"),
              envir = environment())  # initiate local cluster

# Model residual gene expression as a function of sex
emma_Lifestyle <- data.frame(t(parApply(clus, resids, 1, function(y){
  library(EMMREML)
  emma <- emmreml(y = y,  # model each gene
                  X = model.matrix(~Together2$sex + Together2$age + Together2$Flow_PC1 + Together2$Flow_PC2 + Together2$Flow_PC3),  # design matrix
                  Z = as.matrix(mat),  # identity matrix
                  K = as.matrix(GRM2),  # relatedness matrix
                  varbetahat = T, varuhat = T, PEVuhat = T, test = T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
})))

stopCluster(clus)

#get FDR for each p-value
LifestyleSexFDR <- unname(unlist(emma_Lifestyle[,"V14"]))
fdr <- p.adjust(LifestyleSexFDR, method = "fdr")
emma_Lifestyle$pdiff.FDR <- fdr

write.table(emma_Lifestyle, "BaboonLMEOutput2-1-2023.txt", quote=F, sep = "\t")

#figures from raw data
qq(emma_Lifestyle$V14)

hist(emma_Lifestyle$V14)

sum(emma_Lifestyle$pdiff.FDR < 0.05)

#Sex permutation
GeomPolyTable <- data.frame(Pvals=numeric(),
                            Permutation=character(), 
                            stringsAsFactors=FALSE) 

n=1

#from here on permutation
for(j in 1:10){
  ncores <- detectCores(logical = TRUE)  # find number of available cores
  clus <- makeCluster(ncores, setup_strategy = "sequential")  # create cluster
  registerDoParallel(cores = ncores)
  clusterExport(clus,
                varlist = c("resids", "Together2", "mat", "GRM2"),
                envir = environment())  # initiate local cluster
  
  Together2$SexPermute <- sample(Together2$sex, nrow(Together2), replace = FALSE, prob = NULL)
  
  emma_Lifestyle2 <- data.frame(t(parApply(clus, resids, 1, function(y){
    library(EMMREML)
    emma <- emmreml(y = y,  # model each gene
                    X = model.matrix(~Together2$SexPermute + Together2$age + Together2$Flow_PC1 + Together2$Flow_PC2 + Together2$Flow_PC3),  # design matrix
                    Z = as.matrix(mat),  # identity matrix
                    K = as.matrix(GRM2),  # relatedness matrix
                    varbetahat = T, varuhat = T, PEVuhat = T, test = T)
    p=emma$pvalbeta
    varb=emma$varbetahat
    b=emma$betahat
    return(c(b,varb,p[,"none"]))
  })))
  
  PTable <- emma_Lifestyle2 %>% dplyr::select(V14) 
  
  PTable <- mutate(PTable, Permutation=n)
  
  GeomPolyTable <- rbind(GeomPolyTable,PTable)
  
  n=n+1
  stopCluster(clus)  
  
}
  

GeomPolyTable$Permutation <- as.character(GeomPolyTable$Permutation)
FigOverall <- ggplot(GeomPolyTable, aes(V14)) +
  geom_freqpoly(bins=100)

write.table(GeomPolyTable, "BaboonLMMEmpiricalFDR.txt", quote=F, sep = "\t")
