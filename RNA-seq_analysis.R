#Set Working Directory And Load Packages

getwd()
setwd("C:/Users/Laven/Documents/Data_Analysis/Single-cell_transcriptomics_identifies_an")


BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install("vsn")
library(vsn)

install.packages("hexbin")

BiocManager::install("EnhancedVolcano")


library(hexbin)

library(factoextra)
library(ggplot2)

list.files()

#Create and Inspect Matricies From Data Tables
gene_exp_matrix <- read.table('NCOMMS-19-7936188_bulk_RNAseq_raw_counts.txt', 
                        header = TRUE, sep = '\t')

dim(gene_exp_matrix)

gene_exp_matrix[1:4, 1:4]

pheno_matrix <- read.table('NCOMMS-19-7936188_bulk_RNAseq_metadata.txt', 
                           header = TRUE, sep = '\t', stringsAsFactors = TRUE)

pheno_matrix[1:4, 1:4]

#Rename rows of pheno matrix to sample_id values
rownames(pheno_matrix) <- pheno_matrix$sample_id

#Inspect pheno matrix
dim(pheno_matrix)
head(pheno_matrix)

#Check that the matricies are properly aligned
all(rownames(pheno_matrix) == colnames(gene_exp_matrix))

#Setting Sample Parameters To Focus on to inspect

stimtime1 <- '5d'
conditions1 <- c('Th2', 'Th0')
celltype1 <- 'CD4_Memory'

stimtime2 <- '5d'
conditions2 <- c('Th2', 'Th0')
celltype2 <- 'CD4_Naive'

# Create indices for subsetting based on conditions

toSelect1 <- pheno_matrix$stimulation_time == stimtime1 & 
  pheno_matrix$cytokine_condition %in% conditions1 &
  pheno_matrix$cell_type == celltype1

toSelect2 <- pheno_matrix$stimulation_time == stimtime2 & 
  pheno_matrix$cytokine_condition %in% conditions2 &
  pheno_matrix$cell_type == celltype2

#Check that these indices work

toSelect1 #Yes
toSelect2 #Yes

#Create Subsets of the phenotype matrix and gene expression matrix for the CD4_Memory Cells

pheno_matrix.subset <- pheno_matrix[toSelect1, ]
gene_exp_matrix.subset <- gene_exp_matrix[ , toSelect1]

#Create Subsets of the phenotype matrix and gene expression matrix for the CD4_Naive Cells

pheno_matrix.subset2 <- pheno_matrix[toSelect2, ]
gene_exp_matrix.subset2 <- gene_exp_matrix[ , toSelect2]

# Create a DESeq2 Object for differential expression analysis. This is a key data 
#structure used in the DESeq2 package for RNA-Seq data analysis. 
#It encapsulates all the information necessary for differential expression analysis, 
#including the raw count data, sample information (metadata), 
#and the experimental design. 

#This one is for the CD4_Memory Cells
dds <- DESeqDataSetFromMatrix(countData = gene_exp_matrix.subset, 
                              colData = , pheno_matrix.subset,
                              design = ~ cytokine_condition)

#What's going on with the CD4_Naive Cells
dds2 <- DESeqDataSetFromMatrix(countData = gene_exp_matrix.subset2, 
                               colData = , pheno_matrix.subset2,
                               design = ~ cytokine_condition)

#Inspecting the components of the DESeqDataSet objects

#CD4_Memory
print(dds)
counts(dds)
counts(dds, normalized = TRUE)
sizeFactors(dds)
colData(dds)
design(dds)
resultsNames(dds)

#CD4_Naive
print(dds2)
counts(dds2)
counts(dds2, normalized = TRUE)
sizeFactors(dds2)
colData(dds2)
design(dds2)
resultsNames(dds2)

#Filter out genes with less than 10 reads (memory)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#Filter out genes with less than 10 reads (naive)

keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]
dds2

# Apply a pseudocount of 1 and apply log2 (memory)

#Adding a pseudocount of 1 to each raw count ensures that there are no zero counts 
#in the data. This is important because the log of zero is undefined. 
#By adding a small pseudocount, you avoid taking the log of zero.

#Doing the log2 transformation to the counts stabilizes 
#the variance across the range of counts.

normtfd <- normTransform(dds)

# Apply a pseudocount of 1 and apply log2 (naive)

normtfd2 <- normTransform(dds2)

#This creates a plot that shows the relationship between the mean and 
#the standard deviation of the transformed data.
#It helps to assess if the variance is stabilized after transformation.

# Compare mean to sd (memory)
p1 <- meanSdPlot(assay(normtfd))
p1$gg + ggtitle("Mean-SD Plot For CD4 Memory Cells") + 
  xlab("Rank(Mean)") + 
  ylab("Standard Deviation") +
  theme(plot.title = element_text(hjust = 0.5))

# Compare mean to sd (naive)

p2 <- meanSdPlot(assay(normtfd2))
p2$gg + ggtitle("Mean-SD Plot For CD4 Naive Cells") + 
  xlab("Rank(Mean)") + 
  ylab("Standard Deviation") +
  theme(plot.title = element_text(hjust = 0.5))

# Let's calculate rlog values and take another look. (memory)
#This should provide more effective stabilization

rltfd <- rlog(dds, blind=FALSE)

p3 <- meanSdPlot(assay(rltfd))
p3$gg + ggtitle("Mean-SD Plot for Memory Cells (rlog method)") + 
  xlab("Rank(Mean)") + 
  ylab("Standard Deviation") +
  theme(plot.title = element_text(hjust = 0.5))

# Let's calculate rlog values and take another look. (naive)

rltfd2 <- rlog(dds2, blind=FALSE)

p4 <- meanSdPlot(assay(rltfd))
p4$gg + ggtitle("Mean-SD Plot for Naive Cells (rlog method)") + 
  xlab("Rank(Mean)") + 
  ylab("Standard Deviation") +
  theme(plot.title = element_text(hjust = 0.5))

#DESeq2 uses the median of ratios method where the counts are
#divided by sample-specific size factors determined by median ratio of
#gene counts relative to geometric mean per gene. (memory)

#Size factors in RNA sequencing  help us adjust the data so that 
#we can compare the number of gene counts each sample 
#has, even though the samples might have been processed or sequenced 
#differently.

dds <- estimateSizeFactors(dds)

sizeFactors(dds)


plot(sizeFactors(dds), 
     colSums(counts(dds, normalized=F)), 
     xlab = 'Size Factor',
     ylab = 'Total Number Of Reads',
     main = 'Comparing Size Factors To Total Number Of Reads For Different Samples Of CD4 Memory Cells',
     pch = 19)

#gene counts relative to geometric mean per gene. (naive)

dds2 <- estimateSizeFactors(dds2)

sizeFactors(dds2)


plot(sizeFactors(dds2), 
     colSums(counts(dds2, normalized=F)), 
     xlab = 'Size Factor',
     ylab = 'Total Number Of Reads',
     main = 'Comparing Size Factors To Total Number Of Reads For Different Samples Of CD4 Naive Cells',
     pch = 19)

#Samples with more reads have larger size factors so that
#dividing the counts by the size factors accounts for the differences in
#sequencing depth between samples. (memory)

rltfd.pca <- prcomp(t(assay(rltfd)), scale = TRUE)

class(rltfd.pca)

#Samples with more reads have larger size factors so that
#dividing the counts by the size factors accounts for the differences in
#sequencing depth between samples. (naive)

rltfd2.pca <- prcomp(t(assay(rltfd2)), scale = TRUE) #Error... too many genes with zero variance

zero_variance_genes <- which(apply(assay(rltfd2), 1, var) == 0) #Identifying the zero variance genes

filtered_assay <- assay(rltfd2)[-zero_variance_genes, ] #removing the zero variance genes
class(filtered_assay)

rltfd2.pca <- prcomp(t(filtered_assay), scale = TRUE) #performing PCA on the new matix

class(rltfd2.pca) #testing that it is the proper class

# Extract the principal component scores (memory)
pca_scores <- rltfd.pca$x

# Extract the principal component scores (naive)
pca_scores2 <- rltfd2.pca$x

# Convert to a data frame and inspect the head (memory)
pca_df <- as.data.frame(pca_scores)
head(pca_df)

# Convert to a data frame and inspect the head (naive)
pca_df2 <- as.data.frame(pca_scores2)
head(pca_df2)


require(factoextra)

p5 <- fviz_eig(rltfd.pca) #Create Scree Plot (First PC is 31.9% of the variance and second is 22.7% of the variance)
p5 + ggtitle("Scree Plot for CD4 Memory Cells") +
  xlab("Principle Component") + 
  ylab("Percentage Of Explained Variances") +
  theme(plot.title = element_text(hjust = 0.5))

p6 <- fviz_pca_ind(rltfd.pca) #Create PCA Plot
p6 + ggtitle("Individuals - PCA - CD4 Memory Cells") +
  xlab("PC1 (31.9%)") +
  ylab("PC2 (22.7%)") +
  theme(plot.title = element_text(hjust = 0.5))
p6

# Extract the principal component scores (naive)
pca_scores2 <- rltfd2.pca$x


p7 <- fviz_eig(rltfd2.pca) #Create Scree Plot (First PC is 30.2% of the variance, second PC is 19% of the variance)
p7 + ggtitle("Scree Plot for CD4 Naive Cells") +
  xlab("Princple Component") +
  ylab("Percentage Of Explained Variances") +
  theme(plot.title = element_text(hjust = 0.5))


p8 <- fviz_pca_ind(rltfd2.pca) #Create PCA Plot 
p8 + ggtitle("Individuals - PCA - CD4 Naive Cells") +
  xlab("PC1 (30.2%)") +
  ylab("PC2 (19%)") +
  theme(plot.title = element_text(hjust = 0.5))

#Check if our samples for C4 Memory Cells group by sequencing batch which would indicate a technical confounding factor.
colData(rltfd)$sequencing_batch <- factor(colData(rltfd)$sequencing_batch)
#This does not seem to be the case

# Generate PCA plot
pca_plot7 <- plotPCA(rltfd, intgroup = 'sequencing_batch', ntop = 26656, returnData = TRUE)
pca_plot7
percentVar7 <- round(attr(pca_plot7, "percentVar") * 100, 1)

# Create the PCA plot with ggplot2
p7 <- ggplot(pca_plot7, aes(PC1, PC2, color = sequencing_batch)) +
  geom_point(size = 3) +
  ggtitle("Checking If PCA Results Vary Between Sequencing Batch For CD4 Memory Cells") +
  xlab(paste0("PC1: 31.9%")) +
  ylab(paste0("PC2: 22.7%")) +
  labs(color = "Sequencing Batch") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot
print(p7)
#They do not cluster by sequencing batch

#Check if our samples for C4 Naive Cells group by sequencing batch which would indicate a technical confounding factor.
colData(rltfd2)$sequencing_batch <- factor(colData(rltfd2)$sequencing_batch); plotPCA(rltfd2, intgroup = 'sequencing_batch', ntop = 26656)
#This does not seem to be the case

# Generate PCA plot
pca_plot8 <- plotPCA(rltfd2, intgroup = 'sequencing_batch', ntop = 26656, returnData = TRUE)
percentVar8 <- round(attr(pca_plot8, "percentVar") * 100, 1)

# Create the PCA plot with ggplot2
p8 <- ggplot(pca_plot8, aes(PC1, PC2, color = sequencing_batch)) +
  geom_point(size = 3) +
  ggtitle("Checking If PCA Results Vary Between Sequencing Batch For CD4 Naive Cells") +
  xlab(paste0("PC1: 30.2%")) +
  ylab(paste0("PC2: 19%")) +
  labs(color = "Sequencing Batch") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot
print(p8)


#This command is used to open the help page for the plotPCA function. 
#It provides details about how to use the function, its parameters, and examples.
?plotPCA #It will display in the help pane.

#This command retrieves the specific method implementation of the plotPCA function 
#for objects of the class DESeqTransform. #It the actual code that defines how plotPCA works 
#for DESeqTransform objects, 
#which can help you understand the function's internals or modify it if needed.
getMethod("plotPCA","DESeqTransform")


#Set up the parameters for the plotPCA function
object=rltfd
object2=rltfd2
intgroup = 'sequencing_batch'
ntop=26656  
returnData = FALSE

rv <- rowVars(assay(object)) #Computes the row variances of the assay data
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))] #Selects the top ntop genes with highest variance

pca <- prcomp(t(assay(object)[, ]), scale = TRUE) #Performs PCA on the selected genes for CD4 Memory

pca$sdev #check pca standard deviation
sum(pca$sdev) #check sum of standard deviation

percentVar <- pca$sdev^2 / sum(pca$sdev^2)

percentVar

pca2 <- prcomp(t(assay(object2)[, ]), scale = TRUE) #Performs PCA on the selected genes for CD4 Naive

pca2$sdev #check pca standard deviation
sum(pca2$sdev) #check sum of standard deviation

percentVar2 <- pca2$sdev^2 / sum(pca2$sdev^2)

percentVar2

if (!all(intgroup %in% names(colData(object)))) {
  stop("the argument 'intgroup' should specify columns of colData(dds)")
}
intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE]) #Creates a data frame of the grouping factors.

group <- if (length(intgroup) > 1) {
  factor(apply(intgroup.df, 1, paste, collapse = ":"))
} else {
  colData(object)[[intgroup]]
}

#Create data frame of the principle components (CD4 memory cells)
d <- data.frame(
  PC1 = pca$x[, 1], 
  PC2 = pca$x[, 2], 
  sequencing_batch = factor(colData(object)[[intgroup]]), 
  cytokine_condition = colData(object)$cytokine_condition,
  name = colnames(object)
)

#Create data frame of the principle components (CD4 naive cells)
d2 <- data.frame(
  PC1 = pca2$x[, 1], 
  PC2 = pca2$x[, 2], 
  sequencing_batch = factor(colData(object2)[[intgroup]]),
  cytokine_condition = colData(object2)$cytokine_condition, 
  name = colnames(object2)
)


if (returnData) {      #if returnData is set to TRUE add the percent variance to the data frame which will be plotted with ggplot
  attr(d, "percentVar") <- percentVar[1:2]
  return(d)
}



#Both the plot created in ggplot and created with plotPCA do show the same results.
#With both CD4 Memory and CD4 Naive the sequencing reads is not a confounding factor.

#Check if our samples group by treatment
# Generate PCA plot data
pca_plot9 <- plotPCA(rltfd, intgroup = 'cytokine_condition', ntop = 26656, returnData = TRUE)

# Create the PCA plot with ggplot2 and add a title
p9 <- ggplot(pca_plot9, aes(PC1, PC2, color = cytokine_condition)) +
  geom_point(size = 3) +
  ggtitle("Examining if CD4 Memory Cells Cluster By Cytokine Treatment") +
  xlab(paste0("PC1: 31.9%")) +
  ylab(paste0("PC2: 22.7%")) +
  labs(color = "Cytokine Condition") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot
print(p9)
#No they do not for CD4 Memory Cells.

pca_plot10 <- plotPCA(rltfd2, intgroup = 'cytokine_condition',  ntop = 26656, returnData = TRUE)

p10 <- ggplot(pca_plot10, aes(PC1, PC2, color = cytokine_condition)) +
  geom_point(size = 3) +
  ggtitle("Examining if CD4 Naive Cells Cluster By Cytokine Treatment") +
  xlab(paste0("PC1: 30.2%")) +
  ylab(paste0("PC2: 19%")) +
  labs(color = "Cytokine Condition") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot
p10
#Yes they do for CD4 Naive Cells.

# Differential Expression.

## Load the libraries

library(EnhancedVolcano) #Contains the enhanced volcano function for making volcano plots
library(pheatmap)

## DESeq object for CD4 Memory

dds <- DESeq(dds)

class(dds)

res <- results(dds)
dim(res)

res ##This data frame shows the difference in expression between the Th2 condition 
#and the Th0 condition for the CD4 Memory Cells

##DESeq object for CD4 Naive

dds2 <- DESeq(dds2)

class(dds2)

res2 <- results(dds2)
dim(res2)

res2 ##This data frame shows the difference in expression between the Th2 condition 
#and the Th0 condition for the CD4 Naive Cells

## Analysis of the outcome

#How many significantly differentially expressed genes do we find for the
#current contrast Th2 vs Th0 in CD4+ memory cells?

sum(res$padj <= 0.01 & 
      abs(res$log2FoldChange) > 1, na.rm = TRUE)

## Visualization 1: Volcano Plot

#Volcano plots are a helpful tool to visualize the log-fold changes and
#corresponding differential expression p-values

p11 <- EnhancedVolcano(res, lab = rownames(res),
                x = 'log2FoldChange', y = 'padj', 
                subtitle = 'Th2 vs Th0 for CD4 Memory Cells', labSize = 3, 
                pCutoff = 0.01,
                FCcutoff = 1,
                drawConnectors = TRUE) +
                theme(plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5))

p11


p12 <- EnhancedVolcano(res2, lab = rownames(res2),
                x = 'log2FoldChange', y = 'padj', 
                subtitle = 'Th2 vs Th0 for CD4 Naive Cells', labSize = 3, 
                pCutoff = 0.01,
                FCcutoff = 1,
                drawConnectors = TRUE) +
                theme(plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5))

p12

#Counting significantly differentially expressed genes
significant_genes_CD4_Memory <- res[which(res$padj < 0.01 & abs(res$log2FoldChange) > 1), ]
num_significant_genes_CD4_Memory <- nrow(significant_genes_CD4_Memory)
print(num_significant_genes_CD4_Memory)

significant_genes_CD4_Naive <- res2[which(res2$padj < 0.01 & abs(res2$log2FoldChange) > 1), ]
num_significant_genes_CD4_Naive <- nrow(significant_genes_CD4_Naive)
print(num_significant_genes_CD4_Naive)

## Visualization 2: Heatmap

#Lastly, it can be helpful to visualize the individual expression values
#for a set of genes of interest over the different samples. This is
#commonly done using heatmaps.

#First, we select our genes of interest, here the differentially
#expressed genes for CD4 Memory

DEG.idx <- which(res$padj <= 0.01 & 
                   abs(res$log2FoldChange) > 1)
res[DEG.idx,]
df <- as.data.frame(colData(dds)[,c("cytokine_condition","donor_id", "sequencing_batch")])

df

#And for CD4 Naive

DEG.idx2 <- which(res2$padj <= 0.01 & 
                   abs(res2$log2FoldChange) > 1)
res2[DEG.idx2,]
df2 <- as.data.frame(colData(dds2)[,c("cytokine_condition","donor_id", "sequencing_batch")])

subset_data2 <- assay(rltfd2)[DEG.idx2, ]
print(subset_data2)

df2

print(DEG.idx2)


#Secondly, we use the pheatmap function. Importantly, this function can
#perform a clustering of rows to group genes with similar expression
#patterns, as well as clustering of columns to group samples with similar
#patterns. Here, we see that samples cluster nicely by treatment
#(cytokine_condition). Also note, that the expression values are scaled
#by row, i.e.gene to compensate for differences in based expression and
#focus on expression changes between samples.

#CD4 Memory Cell Heatmap
p13 <- pheatmap(assay(rltfd)[DEG.idx,], annotation_col=df,
         treeheight_row = 0, treeheight_col = 0, scale = "row",
         main = "CD4 Memory Cell Differential Expression Heatmap")

p13

#CD4 Naive Cell Heatmap

p14 <- pheatmap(assay(rltfd2)[DEG.idx2,], annotation_col=df2,
         treeheight_row = 0, treeheight_col = 0, scale = "row")

#This heatmap is too messy. I will have to make a higher minimum effect size so that I get a cleaner heatmap

DEG.idx2.2 <- which(res2$padj <= 0.01 & 
                    abs(res2$log2FoldChange) > 4)
res2[DEG.idx2.2,]
df2.2 <- as.data.frame(colData(dds2)[,c("cytokine_condition","donor_id", "sequencing_batch")])

p15 <- pheatmap(assay(rltfd2)[DEG.idx2.2,], annotation_col=df2,
         treeheight_row = 0, treeheight_col = 0, scale = "row",
         main = "CD4 Naive Cell Differential Expression Heatmap")

p15






