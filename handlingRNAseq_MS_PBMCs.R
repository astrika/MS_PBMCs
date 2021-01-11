# Astrid M Manuel
# 01/09/2021
# Handling the MS PBMCs to aquire differentially expressed genes

geneEx2 = read.table("C:/Users/astri/Documents/ZhaoLab/MS_PBMCs/GSE136411_Matrix-merged-non-normalized-raw.txt", as.is = T, header = T)

illumina_annotations = read.delim("C:/Users/astri/Documents/ZhaoLab/MS_PBMCs/GPL6104-11576.txt", as.is = T)
idx <- match(geneEx[,1],illumina_annotations[,1])
geneEx[,1] <- (illumina_annotations[idx,12])
genes <- geneEx[,1]

geneEx <- geneEx[,-1]

for (row in 1:nrow(geneEx)){
  geneEx[row, ] <- as.numeric(gsub(",", "", geneEx[row, ]))
}

geneEx_num <- geneEx

#removing duplicated genes
geneEx$Symbol = genes
dupgenes <- geneEx$Symbol[which(duplicated(geneEx$Symbol))]

removeduplicates <- geneEx[!duplicated(geneEx$Symbol),]
nodupGeneRecords <- removeduplicates[,-337]
rownames(nodupGeneRecords) <- removeduplicates$Symbol

for(gene in dupgenes){
  duplicateRecords <- geneEx[geneEx$Symbol == gene, ]
  
  duplicateRecords <- duplicateRecords[,-337]
  duplicateRecords <- apply(duplicateRecords, 2,function(x) as.numeric(x))
  geneProbeAverage <- colMeans(duplicateRecords)
  
  nodupGeneRecords[gene,] <- geneProbeAverage

}

nodupGeneRecords <- apply(nodupGeneRecords, 2,function(x) as.numeric(x))

#quantile normalization by Dave Tang, source: https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#Log 2 transform, this is common for RNA seq analysis:
RNA_log2T <- (log2(rpkm+1))

rpkm <- nodupGeneRecords/1000000

#performing quantile normalization, this is not common for RNA seq analysis, but may help for our analysis:
RNA_normal <- quantile_normalisation(RNA_log2T)

par(mfrow=c(2,1))
boxplot(RNA_log2T[,1:50], main="After log2 Tranformation")
boxplot(RNA_normal[,1:20], main= "After log2 Transform & Quantile Normalization")

ms <- RNA_log2T[,61:242]
controls <- RNA_log2T[, 243:309]


# # Running T-test and p-value adjustments:

# allResults will store p-values from t-test
allResults <- c()

# runs through each GEmatrix and runs t-test for each gene to consider differential gene expression
for(i in 1:nrow(ms)){
  
  # conducting t-test for i'th gene
  result <- t.test(ms[i,], controls[i,], var.eq=TRUE, alternative = "two.sided")
  
  #appending the p-values to allResults varaibale
  allResults[i] <- result$p.value
}

#Aquiring FDR p-values with alpha < 0.05

fdr_pvalues <- p.adjust(allResults, method="fdr")

genes[which(fdr_pvalues < 0.05)]




######### Getting the deferentially expressed genes
# Following the tutorial: https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html#Overview 
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

namesOfsamples_MS <- colnames(ms)
namesOfsamples_controls <- colnames(controls)
all_names <- c(namesOfsamples_MS, namesOfsamples_controls)
status <- c(rep("MS", 182), rep("Control", 67))
sampleInfo <- cbind(all_names, status)
countdata <- data.frame(cbind(ms,controls))

# add the groups (ms vs controls levels)
y <- DGEList(nodupGeneRecords)
group <- x
group <- factor(group)
y$samples$group <- group

# get annotation for genes
ann <- removeduplicates$Symbol
y$genes <- ann

# creating the design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)

#testing for differential expression
fit <- lmFit(v)
cont.matrix <- makeContrasts(RRMSvsHC=RRMS - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

