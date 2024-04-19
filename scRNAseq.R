library(Seurat)
library(tidyverse)

Human_DTC <- Read10X_h5(filename = '../Data/10k_Human_DTC_Melanoma_3p_gemx_Multiplex_count_raw_feature_bc_matrix.h5')
str(Human_DTC)
nrow(Human_DTC)
ncol(Human_DTC)
#view(Human_DTC)

Seurat_obj <- CreateSeuratObject(counts = Human_DTC , project = "DTC_Melanoma" , min.cells = 5 , min.features = 200)
str(Seurat_obj)
nrow(Seurat_obj)
ncol(Seurat_obj)
view(Seurat_obj@meta.data)

# QC
Seurat_obj[["MT_percent"]] <- PercentageFeatureSet(Seurat_obj , patter = "^MT-")
view(Seurat_obj@meta.data)

vln <- VlnPlot(Seurat_obj, features = c("nCount_RNA","nFeature_RNA","MT_percent"), ncol=3)
ggsave("vln_plot.pdf", plot = vln, device = "pdf")
FS <- FeatureScatter(Seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
ggsave("fs_plot.pdf", plot = FS, device = "pdf")

# Filtering
Seurat_obj <- subset(Seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & MT_percent < 5)


