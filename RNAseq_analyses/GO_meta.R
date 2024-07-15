library(biomaRt)

#meta analysis genes
df = read.table("MetaAnalysisAllGenes12-8", header=TRUE)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO value
#GO term for adaptive immune response
gene.data <- getBM(attributes=c('ensembl_gene_id'),
                   filters = 'go', values = 'GO:0002250', mart = ensembl)

#genes associated
DiseaseGenes <- df[(df$Gene %in% gene.data$ensembl_gene_id),]

#genes not associated
NonDiseaseGenes <- df[(!df$Gene %in% gene.data$ensembl_gene_id),]

t.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, var.equal=FALSE)

#GO term for regulation of adaptive immune system
gene.data <- getBM(attributes=c('ensembl_gene_id'),
                   filters = 'go', values = 'GO:0002819', mart = ensembl)

#genes associated
DiseaseGenes <- df[(df$Gene %in% gene.data$ensembl_gene_id),]
nrow(DiseaseGenes)

#genes not associated
NonDiseaseGenes <- df[(!df$Gene %in% gene.data$ensembl_gene_id),]

t.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, var.equal=FALSE)

#GO term for regulation of the immune system
gene.data <- getBM(attributes=c('ensembl_gene_id'),
                   filters = 'go', values = 'GO:0050776', mart = ensembl)

DiseaseGenes <- df[(df$Gene %in% gene.data$ensembl_gene_id),]
nrow(DiseaseGenes)

NonDiseaseGenes <- df[(!df$Gene %in% gene.data$ensembl_gene_id),]

t.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, var.equal=FALSE)

#GO term for parturition
gene.data <- getBM(attributes=c('ensembl_gene_id'),
                   filters = 'go', values = 'GO:0007567', mart = ensembl)

DiseaseGenes <- df[(df$Gene %in% gene.data$ensembl_gene_id),]
nrow(DiseaseGenes)

NonDiseaseGenes <- df[(!df$Gene %in% gene.data$ensembl_gene_id),]

t.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, var.equal=FALSE)

#GO term for female pregnancy
gene.data <- getBM(attributes=c('ensembl_gene_id'),
                   filters = 'go', values = 'GO:0007565', mart = ensembl)

DiseaseGenes <- df[(df$Gene %in% gene.data$ensembl_gene_id),]
nrow(DiseaseGenes)

NonDiseaseGenes <- df[(!df$Gene %in% gene.data$ensembl_gene_id),]

t.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, var.equal=FALSE)

#analysis for all of the above GO terms
gene.data <- getBM(attributes=c('ensembl_gene_id'),
                   filters = 'go', values = c('GO:0007565', 'GO:0007567', 'GO:0050776', 'GO:0002819', 'GO:0002250'), mart = ensembl)

DiseaseGenes <- df[(df$Gene %in% gene.data$ensembl_gene_id),]
nrow(DiseaseGenes)

NonDiseaseGenes <- df[(!df$Gene %in% gene.data$ensembl_gene_id),]

t.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, var.equal=FALSE)

#analysis for pregnancy related GO terms
gene.data <- getBM(attributes=c('ensembl_gene_id'),
                   filters = 'go', values = c('GO:0007565', 'GO:0007567'), mart = ensembl)

DiseaseGenes <- df[(df$Gene %in% gene.data$ensembl_gene_id),]
nrow(DiseaseGenes)

NonDiseaseGenes <- df[(!df$Gene %in% gene.data$ensembl_gene_id),]

t.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, var.equal=FALSE)

#Analysis for immune-system related GO terms
gene.data <- getBM(attributes=c('ensembl_gene_id'),
                   filters = 'go', values = c('GO:0050776', 'GO:0002819', 'GO:0002250'), mart = ensembl)

DiseaseGenes <- df[(df$Gene %in% gene.data$ensembl_gene_id),]
nrow(DiseaseGenes)

NonDiseaseGenes <- df[(!df$Gene %in% gene.data$ensembl_gene_id),]

t.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, var.equal=FALSE)

#GO term for immune response
gene.data <- getBM(attributes=c('ensembl_gene_id'),
                   filters = 'go', values = c('GO:0006955', mart = ensembl))

DiseaseGenes <- df[(df$Gene %in% gene.data$ensembl_gene_id),]
nrow(DiseaseGenes)

NonDiseaseGenes <- df[(!df$Gene %in% gene.data$ensembl_gene_id),]

t.test(DiseaseGenes$Beta, NonDiseaseGenes$Beta, var.equal=FALSE)