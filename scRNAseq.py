# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 18:38:00 2024

@author: Lavender
"""

#In this code I am using the anndata package.

"""
#anndata stands for Annotated Data.
#It is a Python package that allows for represention of gene
#expression values, sample/gene metadata, and even results from a single-cell 
#RNA-seq experiment within a single Python object. The main advantage of this 
#approach is a uniform, coherent representaton of all different types of information. 
#It is useful for h5ad files.

#When you load an .h5ad file in your working directory, 
#you're typically using the anndata library in Python, 
#which is designed for handling annotated data matrices, 
#such as those used in single-cell RNA sequencing.


# This project also uses Scanpy.
# Scanpy is a flexible toolkit for examining single-cell gene expression data.
# It's packed with tools for preprocessing, visualization, and clustering,
# all designed with single-cell RNA sequencing (scRNA-seq) in mind.
# Scanpy works seamlessly with anndata, giving you a solid setup
# for everything you need in scRNA-seq data analysis.
"""
#%% Cell 1: Importing Libraries
# importing anndata,
import anndata as ad

# importing numpy and pandas and os as well
import numpy as np
import pandas as pd
import os
import sys

#We will now import a subset of the single cell data 
#from the T-cell use case. The data are stored in a compressed 
#format that is optimized for the scRNA-seq experiments.

#%% Cell 2: Inspecting Libraries

# List all loaded modules
loaded_modules = sys.modules.keys()
print(loaded_modules)

#%% Cell 3: Set the working directory to the desired path

desired_path = 'C:/Users/Laven/Documents/Data_Analysis/Single-cell_transcriptomics_identifies_an'
os.chdir(desired_path)

#%% Cell 4: importing scRNA-seq data from a .h5ad file

scdr = ad.read_h5ad('scdr.h5ad')
scdr

#Let's dissect this first batch of information. 
#The standard printout of the anndata object shows that 
#the transcriptomics profiles of 5016 cells (denoted as \"observations\", `obs`) 
#are present. Each profile is composed by 20953 gene expression values, or 
#variables in anndata's terminology. We also have the metadata data frame `obs`, 
#which contains several colums: 'cell.type', 'cytokine.condition', 'donor.id', 
#'batch.10X', 'nGene', etc.

#%% Cell 5: checking data in scdr

# matrix containing the gene expression values
scdr.X

# Accessing the observation annotations (metadata for cells)
cell_metadata = scdr.obs

# Accessing the variable annotations (metadata for genes)
gene_metadata = scdr.var

# Accessing a specific column in the cell metadata
cell_types = scdr.obs['cell.type']

#Not much to see here, at a first glance at least. 
#The key point is that expression values in single cell experiments are very sparse, 
#i.e., most of the genes have no expression in most of the cells. 
#In these cases it is much more efficient to store only the values that are 
#different from zeros (approx. 5 millions values in our case),
#rather than the whole matrix with 5016 x 20953 = 105100248

#%% Cell 6: Let's check what are the highest levels of expressions:

# first 10 highest expression values.
np.flip(np.sort(scdr.X.data))[range(10)]

#Let's now check the 'cell_metadata` data frame

#cell_metadata data frame
cell_metadata

# Display the first few rows of the observation metadata
print(cell_metadata.head())

#We see that the highest expression values across all genes and genes are between 800 / 1100 reads.


#The `cell_metadata` data frame contains one row for each cell, and each column provides 
#information on a different aspects of the cells. For example, the 
#\"cell.type\" column allows to distinguish between Naive and Memory T-cells, 
#while \"cluster.id\" identifies the specific cell type more in detail. Usually, 
#these types of information are not readily available when scRNA-seq data are 
#produced; in this case, they were derived in the original study.

#On a side note, the anndata package using pandas for representing 
#tabular data makes a lot of sense, given the versatility of pandas data frames. 
#In general, one advantage of Python and of programming in general is that well-written, 
#useful modules can then be reused as building-block for future packages.

#%% Cell 7: #### Using an anndata object

#The information in `obs` can be used for narrowing down the scope of our 
#analysis by subsetting the anndata object. We will now retain only the Naive cells:
    
# selecting only Naive T-cells",
ids = scdr.obs['cell.type'] == 'Naive'

# Subset the AnnData object
scdr = scdr[ids, :]


scdr

#We see that `scdr` now contains only the 2141 Naive cells. Let's check the `obs` data frame:
    
# how many rows on obs are left now?
scdr.obs

#set naive cell metadata
cell_metadata_naive = scdr.obs

#Here an important point: *subsetting the anndata object with the command 
#`scdr = scdr[ids, :]` has modified both the gene expression matrix X and 
#the metadata data frame `obs`*. This is the essence of having a single, 
#unified object for representing all information produced in a scRNA-seq 
#experiments: changes are automatically propagated through all the 
#relevant data and metadata in a coherent way.

#Let's now add another type of metadata to the `scdr` object: 
#the `var` data frame that contains information on each single gene. 
#For this, we will first compute the total number of reads for each gene across all cells

# summing over the first axis / dimension",
num_reads = scdr.X.sum(axis=0)
num_reads

#Let's now create the `var` data frame:

# creating the vars data frame",
scdr.var = pd.DataFrame({'num_reads':np.array(num_reads).flatten()}, index=scdr.var_names)
scdr.var

#We now have one more metadata information, this time focusing on 
#characterizing genes, not cells. This change is of course also reflected 
#in the anndata object standard printout, with `var` appearing right below `obs`.

# showing scdr content
scdr

#Many more types of information can be stored within an anndata object. 
#We will see most of them during the next coding lessons, when we will start preprocessing the data.

#### Reading / writing anndata objects

#We already used the `read_h5ad` function for reading an anndata object 
#saved in a binary file based on the hdf5 format. 
#You can also use the `write_h5ad` function for creating a similar 
#file for our current `scdr` object:
    
# creating a directory for storing the new data",
if not os.path.exists('output_data'):
    os.makedirs('output_data')

# writing h5ad objects
scdr.write_h5ad('output_data/scdr_naive_cells.h5ad')

# check the content of the directory
os.listdir('output_data')

#The hdf5 format is unparalled in terms of reading / writing speed, 
#especially for very large files. However, if you want to store your 
#anndata object in comma separate values (csv) files, you can use the `write_csvs` function:

#%% Cell 8: Writing and listing all the written data

# using the write_csv function:
scdr.write_csvs('output_data/csvs', skip_data=True)
# listing all written files
os.listdir('output_data/csvs')

#This function has created a new folder, \"csvs\", 
#which contain one csv file for each element contained 
#within the `scdr` object. *The X matrix is not written, though, 
#unless you explicitly set the `skip_data` argument to `False`*. 
#This is because writing the gene expression matrix in a csv format 
#would usually take very large amount of disk space. In our case `X` 
#would take ~180 MB in csv format, but larger experiments would easily 
#scale up to several GB.

## Preprocessing scRNA-seq data

#### Why preprocessing the data?

#%% Cell 9: Importing and data and pacakages for preprocessing

#Before starting the analysis, it is important to ensure that (a) no low-quality samples 
#or measurements will negatively affect the results, and (b) the data have been 
#transformed in a way that is appropriate for the analysis to perform.
#We will perform the preprocessing of the T-cell use-case data using the scanpy package. 
#Preprocessing is such an important step in the analysis that scanpy has a dedicate module for it, 
#namely `pp`
#First, let's import the relevant packages and load the data

# packages
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# loading the T-cell data
scdr = sc.read_h5ad('scdr.h5ad')
(scdr.n_obs, scdr.n_vars)


#### Basic quality control

scdr

#%% Cell 10: Filtering Cells

#We will start by removing cells where fewer than 200 expressed genes, 
#as these samples may probably be debris or dead cells. 
#This threshold is purely empirical, and it must be adjusted 
#depending on the specific experiment.

# filtering cells
sc.pp.filter_cells(scdr, min_genes=200)
(scdr.n_obs, scdr.n_vars)

#No cells were eliminated this time!

#Similarly, we will remove genes that are expressed only in 3 or fewer cells. 
#Most likely such genes would not contribute significantly to the analyses, 
#and would only aggravate the computational requirements

#%% Cell 11: Filtering Genes

# filtering genes
sc.pp.filter_genes(scdr, min_cells=3)
(scdr.n_obs, scdr.n_vars)

#Circa 6000 genes were filtered out, for having being detected in very few cells.

#We will now compute the percentage of reads that originates from mitochondrial genes. 
#An excess of mitochondrial reads usually indicate RNA contamination:

# the names of mitochondrial genes starts with 'MT-'. Let's create a boolean 
#column in the var data frame to identify them"
scdr.var['mt'] = scdr.var_names.str.startswith('MT-')
# The following function will compute the percentage of mithocondrial reads for each cell
sc.pp.calculate_qc_metrics(scdr, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#%% Cell 12: Creating violin plots
#Let's plot the distribution of relevant quality control quantities:
    
# quality control violin plots
sc.pl.violin(scdr, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

#The first plot reports how many genes are expressed for each cell. 
#The second plot provides the distribution of the total number of reads across cells, 
#while the third plot illustrate the percentage of mitochondrial reads.
#Notice that all distribution are skewed torwards higher values. 
#In the case of the first two plots, these might be doublets, meaning 
#two or more cells mistakenly sequenced as a single cell. In the case of mitochondrial reads percentage, 
#dots above 5% might correspond to droplets contaminated with external RNA.
#We may want to filter out these anomalous samples. Once again we will define some thresholds, heuristically.

#%% Cell 13: Filtering out doublets:
# filtering out possible doublets
scdr = scdr[scdr.obs.n_genes < 2000, :].copy()


#%% Cell 14: Fitlering out cells affected by contaimination
# filtering out cells possibly affected by contamination

scdr = scdr[scdr.obs.pct_counts_mt < 5, :].copy()
(scdr.n_obs, scdr.n_vars)

# plotting again the quality control violin plots
sc.pl.violin(scdr, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

#"The three distributions appear less skewed, and we can now assume that we discarded the most problematic samples.


#### Standard preprocessing pipeline

#After quality control, preprocessing can start. 
#We usually refer to a preprocessing \"pipeline\" because the 
#whole task is broken down in a number of steps: normalization, 
#variance-stabilization transformation, regressing out confounders,
# scaling. Depending on the specific study at hand, this standard 
#pipeline can be altered in different ways. Let's see how these main steps operates on our data.

#The first step is correcting for the different sequencing depth 
#characterizing each cell. To assess the extent of this issue, 
#let's compute the total number of reads for each cell:
    
#%% Cell 15: Making sum of reads per cell into a histogram
    
# summing reads for each cell,
library_size = scdr.X.sum(axis = 1)

# transforming in a pandas data frame
library_size = pd.DataFrame({'library size':np.array(library_size).flatten()}, index=scdr.obs_names)

# plotting the library size distribution
tmp = library_size.hist(bins=20)

#Most cells have around four thousand reads. But there are cells with over 800. The histogram is slightly skewed right.The disparity
#does not allow for direct comparison of the same gene across different samples: a gene which has the same relative level of expression
#between two cells would appear more expressed in the cell with more reads
#Use the Normalize function to bring down these disparities

#%% Cell 16: Make backup copy

scdr_orig = scdr.copy()

#%% Cell 17: Using the Normalize Total Function to make each cell have approixmately 10,000 reads

# normalization. After normalization each cell will have approximately 1e4 normalized reads",
sc.pp.normalize_total(scdr, target_sum=1e4)

# new library size

library_size = scdr.X.sum(axis = 1)
(library_size.min(), library_size.max())
#We can see that all the cell have a normalized library size equal to 10000, 
#up to rounding error.

#Let's now have a look to what the distribution of the new values look like. 
#We remove the zero values, which would other wise dominate the plot.



#%% Cell 18: Histogram of all values

tmp = np.array(scdr.X.todense()).flatten()
tmp = tmp[tmp > 0]
tmp = pd.DataFrame({'values': tmp}).hist(bins=100)

##Most of the values are very low, quite close to zero, while a small fraction of the values can be quite high. 
##Such a skewed distribution can create issues during the analysis, especially for methods assuming close-to-normalty distributions.

##One possible solution is to log-transform the values, so that to reduce the higher values to a more manageable range. 
##A small offset is added to the values before log-transforming so that to avoid taking the logarith of zero.

#%% Cell 18: Log Transforming the values

sc.pp.log1p(scdr)

#%% Cell 19: Histogram of all values after log transformation

tmp = np.array(scdr.X.todense()).flatten()
tmp = tmp[tmp > 0]
tmp = pd.DataFrame({'values': tmp}).hist(bins=1000)

#This is significantly better. The values are much less dispersed. It is still not a normal distribution.

#Before proceeding further, we will now identify the most variable genes, and we will operate from now on only on this subset. 
#The reason behind this is that these genes are likely to retain most of the signal present in the data, 
#while at the same time this will speed up the following computations.

#%% Cell 20: Identifying genes with the highest variation

sc.pp.highly_variable_genes(scdr)
sc.pl.highly_variable_genes(scdr)

#These two plots convey the same information: when the average expression level is contrasted against expession values' dispersions, 
#some gene are cleary on the top side, meaning that their dispersion is greater than the average one. These are the genes we want to focus on.

#%% Cell 21: What is the percentage of highly variable genes

print('%.1f' % np.round(scdr.var.highly_variable.sum()/scdr.n_vars * 100, 3) + '%')
#10.3%

#%% Cell 22: Current status of the data stored in raw

scdr.raw = scdr.copy()


#%% Cell 23: Restrict the X matrix to highly variable genes

scdr = scdr[:, scdr.var.highly_variable].copy()
scdr.n_vars #Output is "1517" showing that genes are highly variable.
#Much less than the original 20953.

#%% Cell 24: Scaling

sc.pp.scale(scdr, max_value = 10)

#%% Cell 25: Inspecting what the gene PASK looks like after scaling.

tmp = pd.DataFrame({'PASK': scdr[:, 'PASK'].X.flatten()}).hist(bins=20)

#printing the average value and standard deviation of PSK
print(np.round(scdr[:, 'PASK'].X.mean()), np.round(scdr[:, 'PASK'].X.std()))

#%% Cell 26: Saving the data

scdr.write_h5ad('scdr_preprocessed.h5ad')


#%% Cell 27: Loading the data

# Reading the saved .h5ad file
scdr_preprocessed = sc.read_h5ad('scdr_preprocessed.h5ad')

#%% Cell 28: Using Principle Component Analysis (PCA) for identifying a lower dimensional space to adequately represent our data

# PCA. Notice that we now use the tools module, abbreviated as "tl"
#This code will store the results in the AnnData object
sc.tl.pca(scdr, svd_solver='arpack')

#%% Cell 29: Creating a plot to examine how much data is explained with each component

# PCA plot. Notice that we now use the plot module, abbreviated as "pl"

sc.pl.pca_variance_ratio(scdr, log=True)

#%% Cell 30: Clustering cells with scanpy package

#clustering. identification of the closest neighbours. Using the first 30 components
sc.pp.neighbors(scdr, n_neighbors=10, n_pcs=30)

#clustering: applying the Leiden algorithm
sc.tl.leiden(scdr)

#the clusters are recorded in the 'leiden' column of the obs data frame
scdr.obs['leiden']

#We have now grouped are cells in different clusters and this information is 
#recorded in `obs`.
#We will now visualize it using UMAP

#%% Cell 31: Computing and visualizing the UMAP

sc.tl.umap(scdr)
sc.pl.umap(scdr, color = 'leiden')

#The dots in the UMAP represent our single cells. 
#Dots close to each other indicate that the transcriptomics profile of the 
#corresponding cells are similar. Dots are colored according to the clusters 
#found by the Leiden algorithm. Notice how the 7 clusters gracefully separate 
#our cells, with minimal or no overlapping.

#We now need to find the differentially expressed genes characterizing each cluster

#%% Cell 32: Finding differentially expressed genes

# find the differentially expressed genes characterizing each cluster
sc.tl.rank_genes_groups(scdr, 'leiden', method='t-test')

# for each cluster, plot the 25 most differentially expressed genes
sc.pl.rank_genes_groups(scdr, n_genes=25, sharey=False)

#%% Cell 33: Markdown cell of information about cell types

#Markers                            | Cell Type
#---                                |---
#LRRN3, CCR7, SELL                  | Naive T-cells (TN)
#PASK                               | Central Memory T-cells (TCM)\
#IL7R, KLRB1, TNFSF13B              | Effector Memory T-cells (TEM)Accessing
#CL4, GZMH, GZMA, GNLY, NKG7, CST7  | Effector Memory T-cells re-expressing CD45RA (TEMRA)

#%% Cell 34: Create a new column in order to store the found identities
    
# let's create a new column in obs to store the found identities
scdr.obs['cell_identity'] = ''
# assigning cluster 6 cells to the "EMRA T-cell" type"
scdr.obs.loc[scdr.obs.leiden == '6', 'cell_identity'] = 'TEMRA Cells'


#Further crossing known markers and differentially expressed genes allows us to realize the following matching:",
#"- Clusters 0, 3, and 5 are characterized by high levels of LRRN3 or CCR7, and thus are recognized as Naive T-cells",
#"- Clusters 2 is characterized by high levels of PASK, a known marker for Central Memory T-cells",
#"- Cluster 1 is more complex to assign. It shares most of its differentially expressed genes with Cluster 2, and thus it is likely to be formed by Central Memory T-cells \n",
#"- All three TEM markers are highly expressed in Clusters 4.",
#"- Clusters 6 was already matched to TEMRA

#%% Cell 35: Storing the different cell types in the cell identity column

# storing the identified cell types in the \"cell_identity\" column\n",
ids = (scdr.obs.leiden == '0') | (scdr.obs.leiden == '3') | (scdr.obs.leiden == '5')
scdr.obs.loc[ids, 'cell_identity'] = 'Naive T-Cells'
ids = (scdr.obs.leiden == '1') | (scdr.obs.leiden == '2')
scdr.obs.loc[ids, 'cell_identity'] = 'Central Memory T-Cells'
ids = scdr.obs.leiden == '4'
scdr.obs.loc[ids, 'cell_identity'] = 'Effector Memory T-Cells'

# plotting the UMAP using colors for marking cell identities\n",
sc.pl.umap(scdr, color = 'cell_identity')

#%%  Cell 36: Checking what is inside the scdr at the end of the analysis

scdr

#%% Cell 37: Saving Post-processed data

scdr.write_h5ad('scdr_post-processed2.h5ad')



