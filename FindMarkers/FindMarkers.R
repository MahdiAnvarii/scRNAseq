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
ggsave("before_dim_plot.pdf", plot = dimplot, device = "pdf", width = 16 , height = 10)

# Harmony
ifnb.harmony <- ifnb.updated %>%
  RunHarmony(group.by.vars = 'stim' , plot_convergence = FALSE)
ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony , "harmony")
ifnb.harmony.embed[1:20,1:20]

ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony' , dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony' , dims = 1:20) %>%
  FindClusters(resolution = 0.5)

dimplot2 <- DimPlot(ifnb.harmony, reduction = "umap", group.by = 'stim')
ggsave("after_dim_plot.pdf", plot = dimplot2, device = "pdf", width = 16 , height = 10)
dimplot|dimplot2

view(ifnb.harmony@meta.data)
clusters <- DimPlot(ifnb.harmony, reduction = 'umap' , group.by = 'seurat_clusters', label = TRUE)
ggsave("clusters_dim_plot.pdf", plot = clusters, device = "pdf")
conditions <- DimPlot(ifnb.harmony, reduction = 'umap' , group.by = 'stim')
ggsave("conditions_dim_plot.pdf", plot = conditions, device = "pdf")
cc <- conditions|clusters
ggsave("conditions_clusters_dim_plot.pdf", plot = cc, device = "pdf", width = 16 , height = 10)

# to identify cell types that form each cluster
# find all markers
AllMarkers <- FindAllMarkers(ifnb.harmony,
               logfc.threshold = 0.25,
               min.pct = 0.3,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')
str(AllMarkers)
view(AllMarkers)

install.packages('metap')
library(BiocManager)
BiocManager::install('multtest')
library(metap)
library(multtest)
# find conserved markers
markers_cluster3 <- FindConservedMarkers(ifnb.harmony,
                                         ident.1 = 3,
                                         grouping.var = 'stim')
head(markers_cluster3)
view(markers_cluster3)

FCGR3A <- FeaturePlot(ifnb.harmony, features = c('FCGR3A'), min.cutoff = 'q10')
ggsave("FCGR3A_dim_plot.pdf", plot = cc, device = "pdf")

head(Idents(ifnb.harmony))
ifnb.harmony <- RenameIdents(ifnb.harmony, '3' = 'CD16 Mono')
head(Idents(ifnb.harmony))
#single cell signature database, PanglaoDB, CellMarker
DimPlot(ifnb.harmony, reduction = 'umap' , label = T)

Idents(ifnb.harmony) <- ifnb.harmony@meta.data$seurat_annotations
head(Idents(ifnb.harmony))
CellTypes <- DimPlot(ifnb.harmony, reduction = 'umap' , label = T)
ggsave("cell_types_dim_plot.pdf", plot = CellTypes, device = "pdf", width = 16 , height = 10)