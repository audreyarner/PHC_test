library(edgeR)
library("readxl")
library(sva)
library(dplyr)
library(ggplot2)
library(plyr)
library(EMMREML)
library(doParallel)
library(biomaRt)
library(qqman)

#Harrison data with covariate information
control_data <- read.delim("RNASeqIndividDatawithBloodCount.txt", header=T, sep="\t")

cpm2 <- read.delim("GenesUsedTurkanaNoCellCount.txt", sep="\t")

#Take only those individs who are strictly Urban or urban
control_data <- control_data %>% filter(Lifestyle=="Urban")

#make sure everything is a factor
control_data$Lifestyle <- as.factor(as.character(control_data$Lifestyle))
control_data$NumberPreg <- as.numeric(control_data$NumberPreg)

GeneMatrix <- read.delim("24Aug21_feature_counts.txt", header =T)

#take out any duplicated genes
GeneMatrix <- GeneMatrix[!duplicated(GeneMatrix[,c("Gene")]),]

#make gene column the row names
rownames(GeneMatrix) <- GeneMatrix[,1]

#only use genes that are protein coding
OnlyUsed <- subset(GeneMatrix, rownames(GeneMatrix) %in% rownames(cpm2))

ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")

YAnalysis<-getBM(attributes=c('ensembl_gene_id',"hgnc_symbol",'chromosome_name','start_position','end_position'), mart = ensembl, filters="ensembl_gene_id", values = row.names(OnlyUsed))                 

#identify and remove genes on the Y chromosome
NoY <- YAnalysis %>% dplyr::filter(chromosome_name %in% c(1:22) | chromosome_name=='X')

OnlyUsed <- subset(OnlyUsed, rownames(OnlyUsed) %in% NoY$ensembl_gene_id)

TechCovariates1 <- read_excel("24Aug21_summaryLaneBatch.xlsx")
TechCovariates2 <- read_excel("24Aug21_summaryTechnicalCovariates.xlsx")

MergedCovariates <- merge(TechCovariates1, TechCovariates2, by.x="Sample", by.y="Barcode")

MergeAll <- merge(control_data, MergedCovariates)

TurkanaRepeats <- ddply(MergeAll,.(SampleNumber),nrow)

MergeAll2 <- MergeAll

MergeAll2[is.na(MergeAll2)] <- 0


MergeAll2[c(9,12,30)] <- scale(MergeAll2[c(9,12,30)])

MergeAll2$Lane <- as.factor(MergeAll2$Lane)

pretty_up_cols=function(cols){
  
  cols=refactorize(cols)
  
  cols=mean_center_clean(cols,"Age")
  cols=mean_center_clean(cols,"Count")

  return(cols)
}

mean_center_clean=function(cols_local,column){
  index=which(colnames(cols_local)==column)
  
  cols_local[,index]=cols_local[,index]-mean(as.numeric(unlist(cols_local[,index])))
  return(cols_local)
}

refactorize=function(cols_local){
  
  cols_local$Gender=factor(cols_local$Gender)
  
  return(cols_local)
}

MergeAll2 = pretty_up_cols(MergeAll2)

GRM <- read.delim("FullTurkanaGRM.txt", check.names = F)

MergedTurkanaData <- read.delim("TurkanaGRMandRNASeqData.txt", header=T)

MergedTurkanaData <- MergedTurkanaData[((MergedTurkanaData$Sample) %in% MergeAll2$Sample),]

GRM2 <- GRM[,(row.names(GRM) %in% MergedTurkanaData$IID)]

GRM2 <- GRM2[(row.names(GRM2) %in% MergedTurkanaData$IID),]

MergeAll2 <- MergeAll2[((MergeAll2$Sample) %in% MergedTurkanaData$Sample),]

OnlyUsed <- OnlyUsed[, (colnames(OnlyUsed) %in% MergeAll2$Sample), drop = FALSE]

Ordering <- MergedTurkanaData$IID
GRM3 <- GRM2[match(Ordering, colnames(GRM2)),]
GRM3 <- GRM3[,match(Ordering, colnames(GRM3))]

MergeAll2 <- MergeAll2[!duplicated(MergeAll2$Sample), ]

mat=matrix(nrow=nrow(MergedTurkanaData),ncol=length(unique(MergedTurkanaData[,"IID"]))); colnames(mat)=unique(MergedTurkanaData[,"IID"]); rownames(mat)=MergedTurkanaData[,"IID"]
for (r in 1:nrow(mat)) {
  for (c in 1:ncol(mat)) {
    if (rownames(mat)[r]==colnames(mat)[c]) {mat[r,c]=1}else {mat[r,c]=0}} }

#get residuals from limma
design <- model.matrix(~0+Lane, data = MergeAll2)
dge <- DGEList(OnlyUsed)
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=FALSE)
fit <-lmFit(v,design)
resids <- residuals.MArrayLM(object=fit, v)

#mixed effect modeling with emmreml
ncores <- detectCores(logical = TRUE)  # find number of available cores
clus <- makeCluster(ncores, setup_strategy = "sequential")  # create cluster
registerDoParallel(cores = ncores)
clusterExport(clus,
              varlist = c("resids", "MergeAll2", "mat", "GRM3"),
              envir = environment())  # initiate local cluster

# Model residual gene expression as a function of diet in genewise models
emma_Lifestyle <- data.frame(t(parApply(clus, resids, 1, function(y){
  library(EMMREML)
  emma <- emmreml(y = y,  # model each gene
                  X = model.matrix(~MergeAll2$Age + MergeAll2$Gender + MergeAll2$Count + MergeAll2$PC1),  # design matrix
                  Z = as.matrix(mat),  # identity matrix
                  K = as.matrix(GRM3),  # relatedness matrix
                  varbetahat = T, varuhat = T, PEVuhat = T, test = T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
})))

stopCluster(clus)

LifestyleSexFDR <- unname(unlist(emma_Lifestyle[,"V13"]))
fdr <- p.adjust(LifestyleSexFDR, method = "fdr")
emma_Lifestyle$pdiff.FDR <- fdr

qq(emma_Lifestyle$V13)

hist(emma_Lifestyle$V13)

sum(emma_Lifestyle$pdiff.FDR < 0.05)

write.table(emma_Lifestyle, "TurkanaUrbanLMEOutput.txt", quote=F, sep = "\t")


#get empirical FDR
#get empirical FDR
GeomPolyTable <- data.frame(Pvals=numeric(),
                            Permutation=character(), 
                            stringsAsFactors=FALSE) 

n=1

#from here on permutation
for(j in 1:5){
  
  MergeAll2$LifestylePermute <- sample(MergeAll2$Gender, nrow(MergeAll2), replace = FALSE, prob = NULL)
  
  MergeAll2$LifestylePermute <- as.character(MergeAll2$LifestylePermute)
  
  ncores <- detectCores(logical = TRUE)  # find number of available cores
  clus <- makeCluster(ncores, setup_strategy = "sequential")  # create cluster
  registerDoParallel(cores = ncores)
  clusterExport(clus,
                varlist = c("resids", "MergeAll2", "mat", "GRM3"),
                envir = environment())  # initiate local cluster
  
  
  emma_Lifestyle <- data.frame(t(parApply(clus, resids, 1, function(y){
    library(EMMREML)
    emma <- emmreml(y = y,  # model each gene
                    X = model.matrix(~MergeAll2$Age + MergeAll2$Gender + MergeAll2$Count + MergeAll2$PC1),  # design matrix
                    Z = as.matrix(mat),  # identity matrix
                    K = as.matrix(GRM3),  # relatedness matrix
                    varbetahat = T, varuhat = T, PEVuhat = T, test = T)
    p=emma$pvalbeta
    varb=emma$varbetahat
    b=emma$betahat
    return(c(b,varb,p[,"none"]))
  })))
  
  PTable <- emma_Lifestyle %>% dplyr::select(V13) 
  
  PTable <- mutate(PTable, Permutation=n)
  
  GeomPolyTable <- rbind(GeomPolyTable,PTable)
  
  n=n+1
  stopCluster(clus)  
  
}


GeomPolyTable$Permutation <- as.character(GeomPolyTable$Permutation)
FigOverall <- ggplot(GeomPolyTable, aes(V13)) +
  geom_freqpoly(bins=100)

write.table(GeomPolyTable, "TurkanaUrbanLMMempiricalPerm.txt", quote=F, sep = "\t")


