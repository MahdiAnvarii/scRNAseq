library(Seurat)
library(tidyverse)
library(ggplot2)
library(harmony)
library(devtools)
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
AvailableData()
install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source")
library(ifnb.SeuratData)
InstallData("ifnb")
data("ifnb")
ifnb.updated = UpdateSeuratObject(object = ifnb)
str(ifnb.updated)
View(ifnb.updated@meta.data)
dim(ifnb.updated)

# QC
ifnb.updated[["MT_percent"]] <- PercentageFeatureSet(ifnb.updated , patter = "^MT-")
view(ifnb.updated@meta.data)

# Filtering
ifnb.updated <- subset(ifnb.updated, subset = nFeature_RNA > 200 & nCount_RNA > 800 & MT_percent < 5)
str(ifnb.updated)
dim(ifnb.updated)
view(ifnb.updated@meta.data)

### Preprocess workflow
# Normalization
ifnb.updated <- NormalizeData(object = ifnb.updated)
# Collect highly variable features
ifnb.updated <- FindVariableFeatures(object = ifnb.updated)
# Scaling
ifnb.updated <- ScaleData(object = ifnb.updated)
# PCA
ifnb.updated <- RunPCA(object = ifnb.updated)
elbowplot <- ElbowPlot(ifnb.updated)
ggsave("elbow_plot.pdf", plot = elbowplot, device = "pdf")
# Non-Linear
ifnb.updated <- RunUMAP(ifnb.updated , dims = 1:20 , reduction = 'pca')
dimplot <- DimPlot(ifnb.updated, reduction = "umap", group.by = 'stim')
ggsave("before_dim_plot.pdf", plot = dimplot, device = "pdf")

