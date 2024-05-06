library(Seurat)
library(DESeq2)
library(tidyverse)
library(devtools)
setRepositories()
install.packages("ExperimentHub")
install.packages("AnnotationHub")
library(AnnotationHub)
library(ExperimentHub)

experimenthub = ExperimentHub(localHub=FALSE)
query(experimenthub,"Kang")
install.packages("SingleCellExperiment")
library(SingleCellExperiment)
install.packages("muscData")
library(muscData)
sce_obj <- experimenthub[['EH2259']]
sce_obj

Seurat_obj <- as.Seurat(sce_obj , data = NULL)
str(Seurat_obj)
dim(Seurat_obj)
view(Seurat_obj@meta.data)

# QC
Seurat_obj[["MT_percent"]] <- PercentageFeatureSet(Seurat_obj , patter = "^MT-")
view(Seurat_obj@meta.data)

# Filtering
Seurat_obj <- subset(Seurat_obj, subset = nFeature_originalexp > 200 & nCount_originalexp > 800 
                     & MT_percent < 5 & multiplets == 'singlet')
str(Seurat_obj)
dim(Seurat_obj)
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
# Non-Linear
Seurat_obj <- RunUMAP(Seurat_obj , dims = 1:20)
dimplot1 <- DimPlot(Seurat_obj, reduction = "umap", group.by = 'cell', label = TRUE)
ggsave("cell_types_dim_plot.pdf", plot = dimplot1, device = "pdf", width = 16 , height = 10)
dimplot2 <- DimPlot(Seurat_obj, reduction = "umap", group.by = 'stim')
ggsave("conditions_dim_plot.pdf", plot = dimplot2, device = "pdf", width = 16 , height = 10)
cc <- dimplot1|dimplot2
ggsave("conditions_clusters_dim_plot.pdf", plot = cc, device = "pdf", width = 16 , height = 10)
