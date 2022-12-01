##############################################
# This script runs a Differential Expression #
# Gene Analysis with DESeq2 package          #
##############################################

## CLEAN ENVIRONMENT AND INSTALL/SET LIBRARIES
#####################
rm(list=ls())

repos = "http://cran.us.r-project.org"
if ("optparse" %in% row.names(installed.packages())  == FALSE) install.packages("optparse", repos = repos)
if ("gplots" %in% row.names(installed.packages())  == FALSE) install.packages("gplots", repos = repos)
if ("ggplot2" %in% row.names(installed.packages())  == FALSE) install.packages("ggplot2", repos = repos)
if ("RColorBrewer" %in% row.names(installed.packages())  == FALSE) install.packages("RColorBrewer", repos = repos)
if ("cluster" %in% row.names(installed.packages())  == FALSE) install.packages("cluster", repos = repos)
if ("pheatmap" %in% row.names(installed.packages())  == FALSE) install.packages("pheatmap", repos = repos)
if ("grid" %in% row.names(installed.packages())  == FALSE) install.packages("grid", repos = repos)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)
if ("DESeq2" %in% row.names(installed.packages()) == FALSE) BiocManager::install("DESeq2")
if ("AnnotationDbi" %in% row.names(installed.packages()) == FALSE) BiocManager::install("AnnotationDbi")
if ("vsn" %in% row.names(installed.packages())  == FALSE) BiocManager::install("vsn")
if ("EnhancedVolcano" %in% row.names(installed.packages())  == FALSE) BiocManager::install("EnhancedVolcano")
if ("tidyr" %in% row.names(installed.packages())  == FALSE) BiocManager::install("tidyr")
if ("org.Hs.eg.db" %in% row.names(installed.packages()) == FALSE) BiocManager::install("org.Hs.eg.db")
if ("apeglm" %in% row.names(installed.packages()) == FALSE) BiocManager::install("apeglm")

suppressPackageStartupMessages({
  library(BiocManager, quietly = TRUE)
  library(optparse, quietly = TRUE)
  library(vsn, quietly = TRUE)
  library(gplots, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(RColorBrewer, quietly = TRUE)
  library(cluster, quietly = TRUE)
  library(pheatmap, quietly = TRUE)
  library(grid, quietly = TRUE)
  library(BiocManager, quietly = TRUE)
  library(DESeq2, quietly = TRUE)
  library(AnnotationDbi, quietly = TRUE)
  library(EnhancedVolcano, quietly = TRUE)
  library(tidyr, quietly = TRUE)
  library(ggrepel, quietly = TRUE)
  library(cp4p, quietly = TRUE)
  library(FDRestimation, quietly = TRUE)
  orgdb <- "org.Hs.eg.db"
  library(orgdb, quietly = TRUE, character.only = TRUE)
})


## CONFIGURATION
#####################

#Load config file
config <- "C:/Users/CBM/Desktop/BWH_counts/configfile_def.txt" 


#Load count data
sampleTable <- read.table(config, header=TRUE)

cutoff <- 0.05


#Transform data variables to factor
sampleTable$sample <- factor(sampleTable$sample)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$gender <- factor(sampleTable$gender)
sampleTable$age <- factor(sampleTable$age)
sampleTable$PED <- factor(sampleTable$PED)

## DESEQ2
#####################

## Data

#Build model using condition as main variable
data <- DESeqDataSetFromHTSeqCount(sampleTable, directory=".", 
                                   design = ~ age + gender + PED + condition)


#Set reference condition
data$condition <- relevel(data$condition, ref = "Unaffected")


#Pre-filtering (keep genes with counts along all samples >=10)
keep <- rowSums(counts(data)) >= 10
data <- data[keep,]
rm(keep)
#With this filter, the object goes from 61806 elements to 28525 elements


#Run deseq
dds <- DESeq(data)
#Alternative deseq run
#dds <- estimateSizeFactors(data)
#dds <- estimateDispersions(dds)
#dds <- nbinomWaldTest(dds, maxit = 10000)


#Obtain and write a file with normalized counts
ncounts <- counts(dds, normalized=TRUE)
write.table (ncounts, file="definitive_results/tsv/deseq2_normalized_counts_.tsv", quote=FALSE, sep = "\t", col.names=NA) 


## Visualize data (transformations of data)

#Tikhonov Regularized Log Transformation (RLT)
rlt <- rlog(dds, blind = FALSE) #calculates RLT
rltMat <- assay(rlt) #takes the data layer on which RLT is saved


#Variance stabilizing transformation (VST)
vst <- varianceStabilizingTransformation(dds, blind = FALSE) #calculates VST
vstMat <- assay(vst) #takes the data layer on which VST is saved
mm <- model.matrix(~condition, colData(vst))
mat <- limma::removeBatchEffect(vstMat,batch=vst$PED,batch2=vst$gender,batch3=vst$age,design=mm)
assay(vst) <- mat


#PCA plot_labelled
tiff(filename = "definitive_results/plots/deseq2_pca.tiff", units="in", width=5, height=5, res=300)
pca <- plotPCA(vst)
pca + ggtitle("Principal Components Plot") + geom_text_repel(aes(label=colnames(vst)),  size=2)
invisible(dev.off())


#Save levels
levels <- unique(sampleTable$condition)
ll <- length(levels)
l1 <- toString(levels[1])
l2 <- toString(levels[2])


#Execute DEGs analysis
res <- results(dds, contrast=c("condition", l1, l2))
res
res$FoldChange <- 2^res$log2FoldChange  
res <- res[colnames(res)[c(1,7,2:6)]] # order columns


#Add symbol and description to results
symbol <- mapIds(get('org.Hs.eg.db'), keys=row.names(res), column="SYMBOL", 
                 keytype="ENSEMBL", multiVals="first") #to obtain gene symbols
description <- mapIds(get('org.Hs.eg.db'), keys=row.names(res),
                      column="GENENAME", keytype="ENSEMBL", 
                      multiVals="first") #to obtain description
res <- cbind(symbol, res) #to add the symbols to the results file
res$description <- description #to add the descriptions to the results file


## Results

my_colour <- list(df=c(l1="orange", l2="skyblue"))
suffix <- paste(l1, l2, sep="_vs_")

conditions <- c(l1, l2)
conds <- subset(sampleTable, sampleTable$condition %in% conditions)
samples <- conds$sample
df <- data.frame(condition=conds$condition)
rownames(df) <- samples
significant <- subset(res, res$padj < cutoff)
significant <- significant[order(significant$padj),]
subcounts <- subset(ncounts, rownames(ncounts) %in% rownames(significant))
subcounts <- subcounts[,rownames(df)]
lsubcounts <- log2(subcounts+1)


#Write all genes
file <- paste("definitive_results/tsv/all_genes_", suffix, ".tsv", sep="")
write.table (res, file=file, quote=FALSE, sep="\t", col.names=NA) 


#Write significant genes (padj<0.05)
file <- paste("definitive_results/tsv/0.05_sig_padj_", suffix, ".tsv", sep="")
write.table (significant, file=file, quote=FALSE, sep="\t", col.names=NA)


#Write significant genes (pvalue instead of padj)
sig_pval <- subset(res, res$pvalue < cutoff)
file <- paste("definitive_results/tsv/deseq2_sig_padj_", suffix, ".tsv", sep="")
write.table (sig_pval, file=file, quote=FALSE, sep="\t", col.names=NA)


## PLOT RESULTS
#####################

## MA-plot

file <- paste("definitive_results/plots/deseq2_maplot.tiff", sep="")
main <- "MA-plot"

jpeg(filename = file, units="in", width=5, height=5, res=300)
DESeq2::plotMA(res, alpha= 0.05, main=main)
invisible(dev.off())


## Volcano plot

#Pre-filtering padj:
keep_noNA <- !is.na(res$padj)
res_2 <- res[keep_noNA,]
rm(keep_noNA)

#Pre-filtering (remove outliers):
keep_FC <- (res_2$log2FoldChange <= 4.5) & (res_2$log2FoldChange >= -4.5)
res_2 <- res_2[keep_FC,]
rm(keep_FC)

file <- paste("definitive_results/plots/volcano_plot_", suffix, ".tiff", sep="")
jpeg(filename = file, units="in", width=10, height=9, res=300)   
EnhancedVolcano(res_2,
                lab = res_2$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = expression(~Log[2]~Fold~Change),
                ylab = expression(~-Log[10]~italic(pvalue)),
                ylim = c(0, 11),
                xlim = c(-3.5, 3.5),
                title = '',
                subtitle = '',
                pCutoff = 0.0002,
                FCcutoff = 0.3,
                pointSize = 1.5,
                labSize = 5,
                colAlpha = 0.5,
                legendLabels=c('Not sig', expression(~Log[2]~Fold~Change), 'pvalue', expression(pvalue~and~Log[2]~FC)),
                legendPosition = "right",
                drawConnectors = TRUE,
                arrowheads = FALSE,
                widthConnectors = 0.7,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
                )
invisible(dev.off())  


## Heatmap 
hmcol <- colorRampPalette(c("red", "yellow", "blue"))(299)

genes <- 17 #number of significant genes (padj < 0.01)
significant_symbol <- as.character(significant$symbol)
my_colour = list(df=c(l1="skyblue", l2="orange"))

file <- paste("definitive_results/plots/heatmap_", suffix, ".jpeg", sep="")
main <- paste("Heatmap of genes with p value adjusted < 0.01", sep=" ")
jpeg(filename = file, units="in", width=8, height=5, res=300) 
pheatmap(mat=lsubcounts[1:genes,], 
         scale="row", 
         cluster_cols=T, 
         labels_row = significant_symbol,
         cluster_rows=T, 
         legend = T, 
         drop_levels = T, 
         fontsize_row=9, 
         annotation_col=df, 
         annotation_colors=my_colour, 
         show_rownames=T,
         show_colnames=T, 
         annotation_names_col=F, 
         annotation_names_row=F, 
         main=main
         )  
invisible(dev.off())  


##################################################### TESTING ###############################################################









