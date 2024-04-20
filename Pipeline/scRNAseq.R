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

# Normalization
Seurat_obj <- NormalizeData(Seurat_obj)

# Collect highly variable features
Seurat_obj <- FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
topfeatures <- head(VariableFeatures(Seurat_obj), 10)
topfeatures

featuresplot <- VariableFeaturePlot(Seurat_obj)
topfeaturesplot <- LabelPoints(plot = featuresplot, points = topfeatures, repel = TRUE)
ggsave("features_plot.pdf", plot = topfeaturesplot, device = "pdf", height = 6, width = 18)

# Scaling
Genes <- rownames(Seurat_obj)
Genes
Seurat_obj <- ScaleData(Seurat_obj, features = Genes)

view(Seurat_obj@assays$RNA$counts)
view(Seurat_obj@assays$RNA$data)
view(Seurat_obj@assays$RNA$scale.data)

# PCA
Seurat_obj <- RunPCA(Seurat_obj, features = VariableFeatures(Seurat_obj))

print(Seurat_obj[["pca"]], dims = 1:10 , nfeatures = 5)
elbowplot <- ElbowPlot(Seurat_obj)
ggsave("elbow_plot.pdf", plot = elbowplot, device = "pdf")

# Clustering
Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:20)
Seurat_obj <- FindClusters(Seurat_obj,resolution = c(0.01,0.1,0.3,0.5))
view(Seurat_obj@meta.data)

DimPlot(Seurat_obj, group.by = "RNA_snn_res.0.1", label = TRUE)
Idents(Seurat_obj) <- "RNA_snn_res.0.1"

# Non-Linear
Seurat_obj <- RunUMAP(Seurat_obj , dims = 1:20)
dimplot <- DimPlot(Seurat_obj, reduction = "umap")
ggsave("dim_plot.pdf", plot = dimplot, device = "pdf")