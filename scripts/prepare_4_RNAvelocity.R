library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)
library(EnhancedVolcano)
library(pals)
library(DESeq2)
library(biomaRt)
library(SeuratDisk)
library(Azimuth)

setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/carlySullivan')

prefix <- "allDatasets"

# Set path to data
E10.5.wt.data.dir <- file.path(getwd(), "E10.5_WT_inputfiles")
E10.5.ko.data.dir <- file.path(getwd(), "E10.5_cKO_inputfiles")
E11.5.wt.data.dir <- file.path(getwd(), "E11.5_WT_inputfiles")
E11.5.ko.data.dir <- file.path(getwd(), "E11.5_cKO_inputfiles")


# Read the data
data.10.5.wt <- Read10X(data.dir = E10.5.wt.data.dir)
data.10.5.ko <- Read10X(data.dir = E10.5.ko.data.dir)
data.11.5.wt <- Read10X(data.dir = E11.5.wt.data.dir)
data.11.5.ko <- Read10X(data.dir = E11.5.ko.data.dir)

# Create the Seurat objects
E10.5_wt.obj <- CreateSeuratObject(counts = data.10.5.wt)
E10.5_ko.obj <- CreateSeuratObject(counts = data.10.5.ko)
E11.5_wt.obj <- CreateSeuratObject(counts = data.11.5.wt)
E11.5_ko.obj <- CreateSeuratObject(counts = data.11.5.ko)

# Add 'orig.ident' metadata
E10.5_wt.obj$orig.ident <- "E10.5_WT"
E10.5_ko.obj$orig.ident <- "E10.5_KO"
E11.5_wt.obj$orig.ident <- "E11.5_WT"
E11.5_ko.obj$orig.ident <- "E11.5_KO"

# Determine percentage of mitochondrial reads
E10.5_wt.obj[["percent.mt"]] <- PercentageFeatureSet(E10.5_wt.obj, pattern = "^mt-")
E10.5_ko.obj[["percent.mt"]] <- PercentageFeatureSet(E10.5_ko.obj, pattern = "^mt-")
E11.5_wt.obj[["percent.mt"]] <- PercentageFeatureSet(E11.5_wt.obj, pattern = "^mt-")
E11.5_ko.obj[["percent.mt"]] <- PercentageFeatureSet(E11.5_ko.obj, pattern = "^mt-")

# 3. Preprocessing and Integration
# Normalize, find variable features, and scale data separately
E10.5_wt.obj <- NormalizeData(E10.5_wt.obj)
E10.5_wt.obj <- FindVariableFeatures(E10.5_wt.obj)
E10.5_wt.obj <- ScaleData(E10.5_wt.obj)

E10.5_ko.obj <- NormalizeData(E10.5_ko.obj)
E10.5_ko.obj <- FindVariableFeatures(E10.5_ko.obj)
E10.5_ko.obj <- ScaleData(E10.5_ko.obj)

E11.5_wt.obj <- NormalizeData(E11.5_wt.obj)
E11.5_wt.obj <- FindVariableFeatures(E11.5_wt.obj)
E11.5_wt.obj <- ScaleData(E11.5_wt.obj)

E11.5_ko.obj <- NormalizeData(E11.5_ko.obj)
E11.5_ko.obj <- FindVariableFeatures(E11.5_ko.obj)
E11.5_ko.obj <- ScaleData(E11.5_ko.obj)

# Merge the datasets
merged_obj <- merge(E10.5_wt.obj, y = c(E10.5_ko.obj, E11.5_wt.obj, E11.5_ko.obj), add.cell.ids = c("WT-10", "KO-10","WT-11", "KO-11"), project = "E11.5")
merged_obj <- JoinLayers(merged_obj)

# 2. Pre-filtering
# Define thresholds
nFeature_lower_threshold <- 500
nFeature_upper_threshold <- 4000
nCount_lower_threshold <- 500
nCount_upper_threshold <- 10000
percent_mt_threshold <- 5

# datasets filtering
merged_obj   <- subset(merged_obj, subset = nFeature_RNA > nFeature_lower_threshold & nFeature_RNA < nFeature_upper_threshold & nCount_RNA > nCount_lower_threshold & nCount_RNA < nCount_upper_threshold & percent.mt < percent_mt_threshold)

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
file.name <- paste0("allDatasets", "_doublets_scores.png")
png(file.name)
VlnPlot(merged_obj, features = "pANN_0.25_0.09_429", pt.size = 0.1)
dev.off()

# Set a threshold for doublet identification (e.g., 0.3)
doublet_threshold <- 0.3

# Filter out doublets directly using the column name in the subset function
merged_obj <- subset(merged_obj, subset = pANN_0.25_0.09_429 < doublet_threshold)

# Quality control before removing doublets
file.name <- paste0("allDatasets", "_filtered_features_withoutDoublets.png")
png(file.name)
VlnPlot(merged_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 3)
dev.off()

# Feature relationship filtered
file.name <- paste0("allDatasets", "_filtered_counts-vs-mito_withoutDoublets.png")
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
anchors <- FindIntegrationAnchors(object.list = SplitObject(merged_obj, split.by = "orig.ident"), dims = 1:20)

# Integrate data
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:20, k.weight = 50)

# Proceed with downstream analysis on the integrated data
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data, npcs = 20)

VizDimLoadings(integrated_data, dims = 1:2, reduction = "pca")

DimHeatmap(integrated_data, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(integrated_data)

# Find neighbors and clusters
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Run UMAP
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:20)

# Visualize UMAP
file.name <- paste0("allDatasets", "_UMAP_byOrigIdent.png")
png(file.name)
DimPlot(integrated_data, reduction = "umap", group.by = "orig.ident", label = TRUE)
dev.off()

# Plot UMAP with clusters
file.name <- paste0(prefix, "_UMAP_byCluster.png")
png(file.name)
DimPlot(integrated_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
dev.off()

saveRDS(integrated_data, "merged_integrated_obj_all.RDS")
SaveH5Seurat(integrated_data, filename = "merged_integrated_obj_all.h5Seurat", overwrite = TRUE)
as.loom(integrated_data, filename = "merged_integrated_obj_all.loom", overwrite = TRUE)

