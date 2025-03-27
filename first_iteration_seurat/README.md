I'll create a README markdown document that explains the Seurat analysis in a clear, layman-friendly format.
markdownCopy# Single-Cell RNA Sequencing Analysis with Seurat

This repository contains R code for analyzing single-cell RNA sequencing data from mouse embryonic forelimb samples using the Seurat package. The analysis compares wildtype (WT) and knockout (KO) samples from embryonic day 10.5 mice.

## What Does This Code Do?

This analysis workflow takes the output from Cell Ranger (10x Genomics) and performs a comprehensive analysis to identify cell types and gene expression differences between normal and knockout mice. Here's a breakdown of the main steps:

### 1. Data Import and Quality Control

- Loads gene expression data from both wildtype and knockout mice
- Calculates quality metrics for each cell (number of genes detected, total RNA counts, percentage of mitochondrial genes)
- Creates visualizations to assess data quality
- Filters out low-quality cells based on customized thresholds

### 2. Doublet Detection and Removal

- Uses the DoubletFinder algorithm to identify "doublets" (droplets that captured two cells instead of one)
- Removes these doublets to ensure clean, single-cell data

### 3. Data Integration and Normalization

- Normalizes gene expression data to account for differences in sequencing depth
- Identifies highly variable genes that can distinguish between cell types
- Integrates the wildtype and knockout datasets to allow direct comparison
- Scales the data to prepare for dimensional reduction

### 4. Clustering and Visualization

- Performs principal component analysis (PCA) to reduce data dimensions
- Identifies cell clusters using a graph-based clustering approach
- Creates UMAP visualizations to show how cells group together
- Generates plots showing:
  - Cells colored by their sample origin (WT or KO)
  - Cells colored by their cluster identity

### 5. Cell Type Analysis

- Identifies marker genes for each cell cluster
- Creates heatmaps and dot plots to visualize gene expression patterns
- Compares cell type proportions between wildtype and knockout samples
- Performs statistical tests to determine if cell type distributions differ significantly between conditions

### 6. Differential Expression Analysis

- Identifies genes that are expressed differently between wildtype and knockout cells
- Creates volcano plots to visualize differential expression results
- Annotates differentially expressed genes with functional information from public databases
- Exports results to tab-separated files for further analysis

## Outputs

The analysis generates multiple visualizations and data files:
- Quality control plots before and after filtering
- UMAP visualizations showing cell clusters
- Heatmaps of marker gene expression
- Dot plots showing expression of selected genes across cell types
- Stacked bar plots comparing cell type proportions between conditions
- Volcano plots highlighting differentially expressed genes
- Annotated lists of differentially expressed genes with statistical metrics

## Usage

This script is designed to be run in R after Cell Ranger processing has been completed. The input files should be organized in the following directories:
- `E10.5_WT_inputfiles/`: Cell Ranger output for wildtype samples
- `E10.5_cKO_inputfiles/`: Cell Ranger output for knockout samples

The final Seurat object is saved as "merged_integrated_obg_10.5.RDS" for future analysis.
