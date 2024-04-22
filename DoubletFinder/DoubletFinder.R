library(remotes)
install_github('chris-mcginnis-ucsf/DoubletFinder')

library(DoubletFinder)
library(ggplot2)
library(tidyverse)
library(Seurat)

cts <- ReadMtx(mtx = '../Data/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix/matrix.mtx.gz',
        features = '../Data/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix/features.tsv.gz',
        cells = '../Data/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix/barcodes.tsv.gz')
str(cts)
nrow(cts)
ncol(cts)
cts[1:20,1:20]

Seurat_obj <- CreateSeuratObject(counts = cts , project = "PBMC")
str(Seurat_obj)
nrow(Seurat_obj)
ncol(Seurat_obj)
view(Seurat_obj@meta.data)

# QC
Seurat_obj[["MT_percent"]] <- PercentageFeatureSet(Seurat_obj , patter = "^MT-")
view(Seurat_obj@meta.data)

# Filtering
Seurat_obj <- subset(Seurat_obj, subset = nFeature_RNA > 500 & nCount_RNA > 800 & MT_percent < 10)
str(Seurat_obj)
nrow(Seurat_obj)
ncol(Seurat_obj)
view(Seurat_obj@meta.data)

### Preprocess workflow
# Normalization
Seurat_obj <- NormalizeData(object = Seurat_obj)
# Collect highly variable features
Seurat_obj <- FindVariableFeatures(object = Seurat_obj)
# Scaling
Seurat_obj <- ScaleData(object = Seurat_obj)
# PCA
Seurat_obj <- RunPCA(object = Seurat_obj)
elbowplot <- ElbowPlot(Seurat_obj)
ggsave("elbow_plot.pdf", plot = elbowplot, device = "pdf")
# Clustering
Seurat_obj <- FindNeighbors(object = Seurat_obj, dims = 1:20)
Seurat_obj <- FindClusters(object = Seurat_obj)
# Non-Linear
Seurat_obj <- RunUMAP(Seurat_obj , dims = 1:20)
