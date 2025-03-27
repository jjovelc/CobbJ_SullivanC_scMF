# Analysis of Single-Cell RNA Sequencing Data

This document explains in simple terms how we processed and analyzed single-cell RNA sequencing data from mouse embryonic forelimb samples.

## Data Preparation

We began with raw sequencing data from mouse embryonic forelimbs at two stages of development (days 10.5 and 11.5), comparing normal (wild-type) and gene-knockout samples. Here's what we did:

1. **Created a reference map**: We prepared a customized mouse genome reference using the standard mouse genome (mm10) with annotations from GENCODE. This gave us a foundation to accurately map our sequencing reads.

2. **Processed raw data**: We used Cell Ranger software to align the raw sequencing data to our reference genome and count the number of RNA molecules for each gene in each cell.

3. **Special processing for RNA velocity**: We set up a specialized workflow to track RNA splicing, which helps us understand how gene expression changes over time. This involved creating a reference that includes both spliced and unspliced transcripts.

## Quality Control

Once we had our basic data, we cleaned it up to ensure accurate analysis:

1. **Filtered out low-quality cells**: We removed cells with:
  - Fewer than 200 genes (likely empty droplets)
  - More than 5% mitochondrial genes (likely dying cells)
  - More than 6,000 genes (likely multiple cells captured together)

2. **Detected and removed doublets**: We used software called Scrublet to identify and remove cases where two cells were accidentally captured together.

3. **Normalized the data**: We adjusted the data so that differences in sequencing depth between cells wouldn't affect our comparisons.

## Analysis Techniques

With clean data in hand, we performed multiple analyses:

1. **Dimensionality reduction**: We used PCA, UMAP, and t-SNE techniques to visualize the high-dimensional data in 2D space, making patterns easier to see.

2. **Cell clustering**: We used the Leiden algorithm to group similar cells together, helping us identify different cell types.

3. **Differential expression**: We identified genes that were expressed differently between:
  - Different cell clusters
  - Wild-type vs. knockout samples

4. **Cluster composition analysis**: We statistically compared how cell types were distributed between wild-type and knockout samples.

5. **RNA velocity analysis**: We examined the balance between newly-made (unspliced) and mature (spliced) RNA to predict the future state of cells, helping us understand developmental trajectories.

## Visualization and Integration

Finally, we created various visualizations to understand the data:

1. **UMAP plots**: Showed how cells group together based on their gene expression profiles.

2. **Stream plots**: Visualized the direction of cellular development using RNA velocity vectors.

3. **Phase portraits**: Examined the splicing dynamics of specific genes.

4. **Stacked bar plots**: Compared the distribution of cell types across different conditions.

We also integrated data from all four samples (two time points × two conditions) to compare cellular dynamics across development and between wild-type and knockout samples.

All analyses were performed using Python with specialized bioinformatics packages including Scanpy, scVelo, and others.
