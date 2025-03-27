library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)
library(EnhancedVolcano)
library(pals)
library(DESeq2)
library(biomaRt)

setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/carlySullivan')

prefix <- "E10.5"

# Set path to data
E10.5.wt.data.dir <- file.path(getwd(), "E10.5_WT_inputfiles")
E10.5.ko.data.dir <- file.path(getwd(), "E10.5_cKO_inputfiles")

# Read the data
data.wt <- Read10X(data.dir = E10.5.wt.data.dir)
data.ko <- Read10X(data.dir = E10.5.ko.data.dir)

# Create the Seurat objects
E10.5_wt.obj <- CreateSeuratObject(counts = data.wt)
E10.5_ko.obj <- CreateSeuratObject(counts = data.ko)

# Add 'orig.ident' metadata
E10.5_wt.obj$orig.ident <- "E10.5_WT"
E10.5_ko.obj$orig.ident <- "E10.5_KO"

# Determine percentage of mitochondrial reads
E10.5_wt.obj[["percent.mt"]] <- PercentageFeatureSet(E10.5_wt.obj, pattern = "^mt-")
E10.5_ko.obj[["percent.mt"]] <- PercentageFeatureSet(E10.5_ko.obj, pattern = "^mt-")

# Generate the data layer
E10.5_wt.obj <- NormalizeData(E10.5_wt.obj)
E10.5_ko.obj <- NormalizeData(E10.5_ko.obj)

# Inspect quality prior to merging and filtering
# Three features WT
file.name <- paste0(prefix, "_wt_unfiltered_features.png")
png(file.name)
VlnPlot(E10.5_wt.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Three features KO
file.name <- paste0(prefix, "_ko_unfiltered_features.png")
png(file.name)
VlnPlot(E10.5_ko.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Feature relationship unfiltered
# WT counts vs mito
file.name <- paste0(prefix, "_wt_unfiltered_counts-vs-mito.png")
png(file.name)
plot1 <- FeatureScatter(E10.5_wt.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
print(plot1)
dev.off()

# KO counts vs mito
file.name <- paste0(prefix, "_ko_unfiltered_counts-vs-mito.png")
png(file.name)
plot1 <- FeatureScatter(E10.5_ko.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
print(plot1)
dev.off()

# WT counts vs features
file.name <- paste0(prefix, "_wt_unfiltered_counts-vs-features.png")
png(file.name)
plot2 <- FeatureScatter(E10.5_wt.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot2)
dev.off()

# KO counts vs features
file.name <- paste0(prefix, "_ko_unfiltered_counts-vs-features.png")
png(file.name)
plot2 <- FeatureScatter(E10.5_ko.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot2)
dev.off()

# 2. Pre-filtering
# Define thresholds
nFeature_lower_threshold <- 500
nFeature_upper_threshold <- 4000
nCount_lower_threshold <- 500
nCount_upper_threshold <- 10000
percent_mt_threshold <- 5

# WT dataset filtering
E10.5_wt.obj <- subset(E10.5_wt.obj, subset = nFeature_RNA > nFeature_lower_threshold & nFeature_RNA < nFeature_upper_threshold & nCount_RNA > nCount_lower_threshold & nCount_RNA < nCount_upper_threshold & percent.mt < percent_mt_threshold)

# KO dataset filtering
E10.5_ko.obj <- subset(E10.5_ko.obj, subset = nFeature_RNA > nFeature_lower_threshold & nFeature_RNA < nFeature_upper_threshold & nCount_RNA > nCount_lower_threshold & nCount_RNA < nCount_upper_threshold & percent.mt < percent_mt_threshold)

# 3. Preprocessing and Integration
# Normalize, find variable features, and scale data separately
E10.5_wt.obj <- NormalizeData(E10.5_wt.obj)
E10.5_wt.obj <- FindVariableFeatures(E10.5_wt.obj)
E10.5_wt.obj <- ScaleData(E10.5_wt.obj)

E10.5_ko.obj <- NormalizeData(E10.5_ko.obj)
E10.5_ko.obj <- FindVariableFeatures(E10.5_ko.obj)
E10.5_ko.obj <- ScaleData(E10.5_ko.obj)

# Inspect quality prior to merging but after filtering
# Three features WT
file.name <- paste0(prefix, "_wt_filtered_features.png")
png(file.name)
VlnPlot(E10.5_wt.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Three features WT
file.name <- paste0(prefix, "_ko_filtered_features.png")
png(file.name)
VlnPlot(E10.5_ko.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Feature relationship unfiltered
# WT counts vs mito
file.name <- paste0(prefix, "_wt_filtered_counts-vs-mito.png")
png(file.name)
plot1 <- FeatureScatter(E10.5_wt.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
print(plot1)
dev.off()

# KO counts vs mito
file.name <- paste0(prefix, "_ko_filtered_counts-vs-mito.png")
png(file.name)
plot1 <- FeatureScatter(E10.5_ko.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
print(plot1)
dev.off()

# WT counts vs features
file.name <- paste0(prefix, "_wt_filtered_counts-vs-features.png")
png(file.name)
plot2 <- FeatureScatter(E10.5_wt.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot2)
dev.off()

# KO counts vs features
file.name <- paste0(prefix, "_ko_filtered_counts-vs-features.png")
png(file.name)
plot2 <- FeatureScatter(E10.5_ko.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot2)
dev.off()

# Merge the datasets
merged_obj <- merge(E10.5_wt.obj, y = E10.5_ko.obj, add.cell.ids = c("WT", "KO"), project = "E10.5")
merged_obj <- JoinLayers(merged_obj)
# Preprocess the merged data
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj, npcs = 30)

# propose the number of doublets to find
nExp_poi <- round(0.075 * ncol(merged_obj))  

# Run DoubletFinder
merged_obj <- doubletFinder(merged_obj, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi)

# Visualize doublet scores
file.name <- paste0(prefix, "_doublets_scores.png")
png(file.name)
VlnPlot(merged_obj, features = "pANN_0.25_0.09_235", pt.size = 0.1)
dev.off()

# Set a threshold for doublet identification (e.g., 0.3)
doublet_threshold <- 0.3

# Filter out doublets directly using the column name in the subset function
merged_obj <- subset(merged_obj, subset = pANN_0.25_0.09_235 < doublet_threshold)

# Quality control before removing doublets
file.name <- paste0(prefix, "_filtered_features_withoutDoublets.png")
png(file.name)
VlnPlot(merged_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
dev.off()

# Feature relationship filtered
file.name <- paste0(prefix, "_filtered_counts-vs-mito_withoutDoublets.png")
png(file.name)
plot1 <- FeatureScatter(merged_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
print(plot1)
dev.off()

file.name <- paste0(prefix, "_filtered_counts-vs-features_withoutDoublets.png")
png(file.name)
plot2 <- FeatureScatter(merged_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot2)
dev.off()

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = SplitObject(merged_obj, split.by = "orig.ident"), dims = 1:30)

# Integrate data
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:30)

# Proceed with downstream analysis on the integrated data
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data, npcs = 30)

VizDimLoadings(integrated_data, dims = 1:2, reduction = "pca")

DimHeatmap(integrated_data, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(integrated_data)


FeaturePlot(integrated_data, features = c("Id1", "Tubb2b", "Tcf15", "Tubb3", "S100a6"))

# Find neighbors and clusters
integrated_data <- FindNeighbors(integrated_data, dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Run UMAP
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)

# Visualize UMAP
file.name <- paste0(prefix, "_UMAP_byOrigIdent.png")
png(file.name)
DimPlot(integrated_data, reduction = "umap", group.by = "orig.ident", label = TRUE)
dev.off()

# Plot UMAP with clusters
file.name <- paste0(prefix, "_UMAP_byCluster.png")
png(file.name)
DimPlot(integrated_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
dev.off()

# Set identities for DE analysis
#Idents(integrated_data) <- "orig.ident"

# This step requires manual inspection of marker genes
markers <- FindAllMarkers(integrated_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Visualize top markers using DoHeatmap
file.name <- paste0(prefix, "_markers_heatmap.png")
png(file.name, width = 600, height = 900)
DoHeatmap(integrated_data, features = top_markers$gene)
dev.off()

# Example input variables (replace these with actual input handling in Shiny app)
input <- list(
  dotGenes = c(
    "Pib2", "Emid1", "Gal", "Enpp3", "Tbx20", "Adamts15", "Rbp4",
    "Ube2v2", "Adar", "Rnf145", "Shisa5", "Aspa", "Ndp", "Cdc25c", "Cdkn1a", 
    "Ccl25", "Gm13403", "Hdc", "Plek2", "Nnat", "Nrn1", "St7", "Prrx2", 
    "Id1", "Tubb2b", "Tcf15", "Tubb3", "S100a6", "Spink1", "Gm12689",
    "Eif3b", "Eif4a2", "Eef2", "Hist1h2bn", "Hist1h3c", "Hist1h4d",
    "Ptk2b", "Ppbp", "Cldn5", "Hba-a1", "Hbb-bs", "Hbb-bt", "Hba-a2",
    "Hbb-y", "4933411F08Rik", "2600014E01Rik", "Tef", "Tek", "Fap", 
    "Nnat", "Mmp3", "Gm2049", "Gm16291", "Gm12942", "Gm13030", 
    "A230057D06Rik", "Gm15636", "Gm1992", "Gm20536", "Gm28043", 
    "Gm4609", "Gm5334", "Gm5546", "Gm5736", "Gm5866"
  ),
  geneList = "",  # Additional genes as a comma or space-separated string
  dotGroupBy = "seurat_clusters"  # Grouping variable for the dot plot
)

# Combine and clean up the gene list
dotGenes <- unique(c(input$dotGenes, unlist(strsplit(input$geneList, "[,\\s]+"))))
dotGenes <- dotGenes[dotGenes != ""]  # Remove any empty strings

# Grouping variable
group.by <- input$dotGroupBy

# Ensure clusters are assigned
if (!"seurat_clusters" %in% colnames(integrated_data@meta.data)) {
  integrated_data <- FindNeighbors(integrated_data, dims = 1:30)
  integrated_data <- FindClusters(integrated_data, resolution = 0.5)
}

# Set identities to clusters
Idents(integrated_data) <- "seurat_clusters"

# Generate Dot Plot with proper grouping
plot <- DotPlot(integrated_data, features = dotGenes, group.by = group.by) + RotatedAxis()

# Display the plot
print(plot)

# Save Dot Plot to a file (optional)
png("DotPlot_Markers_Per_Cluster.png", width = 800, height = 600)
print(plot)
dev.off()


# Compare cell type proportions
# Ensure clusters are assigned
if (!"seurat_clusters" %in% colnames(integrated_data@meta.data)) {
  integrated_data <- FindNeighbors(integrated_data, dims = 1:30)
  integrated_data <- FindClusters(integrated_data, resolution = 0.5)
}

# Set identities to clusters
Idents(integrated_data) <- "seurat_clusters"

# Compute the cell type proportions for each orig.ident
cell_proportions <- prop.table(table(Idents(integrated_data), integrated_data$orig.ident), margin = 2)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Origin", "Proportion")

# Visualize cell type proportions
p <- ggplot(cell_proportions_df, aes(x = Origin, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = as.vector(polychrome(26))) +
  theme_minimal() +
  labs(title = "Cell Type Proportions in WT vs KO", x = "Origin", y = "Proportion")
png("cellProportions_WT_vs_KO.png", width = 12, height = 8)
print(p)
dev.off()

# Statistical test for differences in proportions
cell_counts <- table(Idents(merged_obj), merged_obj$orig.ident)
chi_sq_test <- chisq.test(cell_counts)
print(chi_sq_test)


# Volcano Plot Function
makeVolcanoPlot <- function(df, vp_file) {
  keyvals <- ifelse(
    (df$avg_log2FC < -0.25 & df$p_val_adj < 0.05), 'forestgreen',
    ifelse((df$avg_log2FC > 0.25 & df$p_val_adj < 0.05), 'firebrick1',
           ifelse((abs(df$avg_log2FC) < 0.25 & df$p_val_adj < 0.05), 'dodgerblue1',
                  'black')))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'firebrick1'] <- 'High'
  names(keyvals)[keyvals == 'dodgerblue1'] <- 'small FC'
  names(keyvals)[keyvals == 'forestgreen'] <- 'Low'
  
  evp <- EnhancedVolcano(
    df,
    lab = rownames(df),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    pCutoff = 0.05,
    FCcutoff = 0.25,
    colCustom = keyvals,
    title = NULL,
    subtitle = NULL,
    colAlpha = 0.8,
    shape = 20,
    pointSize = 3.5,
    labSize = 6
  )
  png(vp_file, width = 800, height = 800)
  print(evp)
  dev.off()
}


# Create remote connection
my_mart <- useEnsembl('ensembl', dataset = "mmusculus_gene_ensembl")

# Make list of attributes to retrieve
attributes <- c("ensembl_transcript_id",
                "ensembl_gene_id",
                "entrezgene_id",
                "external_gene_name",
                "wikigene_description",
                "name_1006",
                "definition_1006",
                "namespace_1003")

# Function to extract records
pullRecords <- function(attributes, mart, filter_values){
  records <- getBM(attributes = attributes, filters = "external_gene_name", values = filter_values, mart = mart)
  return(records)
}

# Extract first hit only
getFirstMatch <- function(records, genes) {
  first_annotation_df <- data.frame()
  for (gene in genes) {
    hits <- which(records$external_gene_name == gene)
    if (length(hits) > 0) {
      first <- hits[1]
      first_annotation <- records[first,]
      first_annotation_df <- rbind(first_annotation_df, first_annotation)
    } else {
      first_annotation <- c(gene, "_", "-", "-", "-", "-", "-", "-")
      first_annotation_df <- rbind(first_annotation_df, first_annotation)
    }
  }
  return(first_annotation_df)
}

renameColumns <- function(df){
  colnames(df)[6] <- "GO_group"
  colnames(df)[7] <- "GO_definition"
  colnames(df)[8] <- "Ontology"
  return(df)
}

# Extract list of genes
genes <- rownames(merged_obj[["RNA"]]$counts)
genes_annotations <- pullRecords(attributes, my_mart, genes)
annotation_df <- getFirstMatch(genes_annotations, genes)
annotation_df <- renameColumns(annotation_df)

# Join layers
merged_obj <- JoinLayers(merged_obj)

# Set the tissue type as identity class
Idents(merged_obj) <- merged_obj$orig.ident


# Define the tissue types for comparison
group1 <- "E10.5_WT"
group2 <- "E10.5_KO"

# Run FindMarkers
de_results <- FindMarkers(merged_obj, ident.1 = group1, ident.2 = group2)
vp_file_name <- paste0(prefix, "_volcanoPlot_Wilcox.png")
makeVolcanoPlot(de_results, vp_file_name)

# Define the annotation function
get_annotation <- function(de_res, annotation_df) {
  annotated_df <- data.frame()
  for (gene in rownames(de_res)) {
    hits <- which(annotation_df$external_gene_name == gene)
    if (length(hits) > 0) {
      first <- hits[1]
      annotation <- annotation_df[first,]
      annotated_row <- cbind(de_res[gene,], annotation)
    } else {
      annotation <- data.frame(
        ensembl_transcript_id = NA, 
        ensembl_gene_id = NA, 
        entrezgene_id = NA, 
        external_gene_name = gene, 
        wikigene_description = NA, 
        GO_group = NA, 
        GO_definition = NA, 
        Ontology = NA
      )
      annotated_row <- cbind(de_res[gene,], annotation)
    }
    annotated_df <- rbind(annotated_df, annotated_row)
  }
  return(annotated_df)
}



new_annotations = get_annotation(de_results, annotation_df)
file.name <- paste0(prefix, "_DE_results_Wilcox.tsv")
write.table(new_annotations, file.name, sep = "\t", quote = F)

saveRDS(merged_obj, "merged_integrated_obg_10.5.RDS")
