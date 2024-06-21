# Transcriptomics_Of_Bulk_RNA_And_Single_Cell_RNA

This is a project where I analyzed Bulk RNA and Single-Cell RNA sequencing data from the study “Single-cell transcriptomics identifies an effectorness gradient shaping the response of CD4+ T cells to cytokines” published in Nature Communications in 2020.

This study analyzes gene expression using Bulk RNA and Single-Cell RNA to quantify gene expression of human naïve and memory CD4+ T cells in response to cytokines. In this study, Th0 refers to the cytokine unstimulated condition, and Th2 refers to a cytokine-stimulated condition.
Using R the Bulk RNA sequencing data was processed using the DESeq2 package. The data was visualized using the vsn, ggplot, EnhancedVolcano, and pheatmap packages.

Using Python the Single-Cell RNA sequencing data was processed using the pandas and scanpy packages and visualized using the matplotlib and scanpy packages.

Principle Component Analysis in R revealed that neither the memory cells nor the naïve cells clustered by sequencing batch, meaning that there was no technical confounding error in sequencing that differentiated the reads. The memory cells did not cluster by cytokine treatment, but the naïve cells did, meaning that in the naïve cell samples there was an overall difference in the transcriptome but not in the memory cells. This is in line with the results of the study which found that the naïve cells but not the stimulated cells can be induced towards the Th2 phenotype. Creating volcano plots and heatmaps confirmed that cytokine treatment induced a much stronger gene expression change in naïve cells than memory cells.

Using Python, the single-cell RNA sequencing data was analyzed to find the number of genes expressed per cell, number of RNA reads per cell, and percentage of mitochondrial reads per cell. Differentially expressed genes were identified and the principle components of gene expressions were calculated. UMAP (Uniform Manifold Approximation and Projection) is a dimension reduction technique used to visualize differences in gene expression between cells. UMAPs were prepared that visualized the cells split into cells clustered by the Leiden algorithm and cells clustered by cell identity.
