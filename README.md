# scRNAseq analysis of mouse forelimbs at 10.5 and 11.5 days
This repository contains the code and documentation for the analysis of single-cell RNA sequencing (scRNA-seq) data from mouse embryonic forelimb samples. The analysis compares wild-type (WT) and conditional knockout (cKO) samples at two developmental stages: embryonic day 10.5 (E10.5) and embryonic day 11.5 (E11.5).

## Analysis Overview

We performed two iterations of data analysis using different computational frameworks:

### Iteration 1: Seurat-based Analysis (R)

Our initial analysis utilized the Seurat framework in R, which provided:

- Robust quality control and data filtering
- Doublet detection and removal using DoubletFinder
- Integration of samples for batch correction
- Cell clustering and visualization with UMAP and t-SNE
- Differential expression analysis between genotypes
- Cell type proportion comparison between WT and cKO

This analysis established the foundational characterization of our dataset and identified key gene expression differences between genotypes.

### Iteration 2: Scanpy/scVelo Analysis (Python)

We performed a second analysis iteration using the Scanpy framework in Python, which offered:

- Similar QC, filtering, and clustering approaches to Seurat
- Improved integration with RNA velocity analysis
- Compatibility with the scVelo package for developmental trajectory inference

The Python-based workflow was particularly valuable for RNA velocity analysis, which uses spliced and unspliced mRNA counts processed from Cell Ranger output using velocyto to infer cellular dynamics and differentiation trajectories.

## Key Components

The repository includes:

- **Data Processing**: Scripts for preparing reference genomes, running Cell Ranger, and processing raw sequencing data
- **Quality Control**: Filtering cells based on gene counts, mitochondrial percentage, and doublet detection
- **Cluster Analysis**: Identification and characterization of cell clusters
- **Differential Expression**: Comparison of gene expression between conditions
- **RNA Velocity**: Analysis of cellular dynamics and developmental trajectories using velocyto and scVelo
- **Visualization**: UMAP plots, heatmaps, and RNA velocity stream plots

## Directories

- [/first_iteration_seurat/](./first_iteration_seurat/): Contains R scripts and outputs from the Seurat-based analysis
- [/Scanpy_Analysis/](./Scanpy_Analysis/): Contains Python scripts and outputs from the Scanpy/scVelo analysis

## Why Two Analysis Approaches?

While Seurat is a powerful and user-friendly tool for scRNA-seq analysis, we complemented it with Scanpy/scVelo for several reasons:

1. **RNA Velocity Integration**: The Scanpy ecosystem integrates seamlessly with scVelo for RNA velocity analysis, which helps infer the future state of cells based on spliced/unspliced mRNA ratios
2. **Computational Reproducibility**: Python-based workflows can be more easily integrated into computational pipelines
3. **Cross-Validation**: Using two independent analytical frameworks provides validation of key findings

Both approaches yielded consistent results regarding cell clustering and differential expression, with the Scanpy/scVelo analysis providing additional insights into developmental trajectories and cellular dynamics.

## Data Processing Workflow

1. **Reference Genome Preparation**: Custom mm10 reference with GENCODE M23 annotations
2. **Read Alignment**: Cell Ranger was used to process raw sequencing data
3. **RNA Velocity Preparation**: Velocyto was used to quantify spliced/unspliced counts from Cell Ranger BAM files
4. **Quality Control**: Filtering based on gene counts, UMI counts, and mitochondrial percentage
5. **Analysis**: Dimensionality reduction, clustering, and differential expression
6. **Trajectory Analysis**: RNA velocity to infer cellular dynamics


## Getting Started

Detailed instructions for each analysis framework are provided in their respective directories:
- See [/first_iteration_seurat/README.md](./first_iteration_seurat/README.md): for the R-based workflow
- See [/second_iteration_scanpy/README.md](./second_iteration_scanpy/README.md): for the Python-based workflow
