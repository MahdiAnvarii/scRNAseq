library(devtools)
devtools::install_github('satijalab/seurat-wrappers')
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)

markers <- read.delim('../Data/ABC_Marker.txt' , header = T) # gene metadata
metadata <- read.delim('../Data/ABC_Meta.txt' , header = T) # cells metadata
expressions <- read.delim('../Data/ABC_umi_matrix_7551_cells.csv' , header = T , sep = ',') # expression matrix
head(markers)
head(metadata)
expressions[1:10,1:10]

exp.t <- t(expressions)
exp.t[1:10,1:10]
Seurat_obj <- CreateSeuratObject(counts = exp.t)
view(Seurat_obj@meta.data)
Seurat_obj@meta.data <- merge(Seurat_obj@meta.data, metadata, by.x = 'row.names' , by.y = 'cell_id')
view(Seurat_obj@meta.data)
Seurat_obj@meta.data <- Seurat_obj@meta.data %>%
  column_to_rownames(var = 'Row.names')
str(Seurat_obj)
dim(Seurat_obj)

# QC
Seurat_obj@meta.data$MT_percent <- PercentageFeatureSet(Seurat_obj , patter = "^MT-")
view(Seurat_obj@meta.data)

# Filtering
Seurat_obj <- subset(Seurat_obj, subset = nFeature_RNA > 200 & nCount_RNA > 800 & MT_percent < 5)
str(Seurat_obj)
dim(Seurat_obj)
view(Seurat_obj@meta.data)

# Choosing B cells
unique(Seurat_obj@meta.data$population)
Idents(Seurat_obj) <- Seurat_obj$population
B_Seurat_obj <- subset(Seurat_obj, idents = 'b')
str(B_Seurat_obj)
dim(B_Seurat_obj)
view(B_Seurat_obj@meta.data)
unique(B_Seurat_obj@meta.data$redefined_cluster)

### Preprocess workflow
# Normalization
B_Seurat_obj <- NormalizeData(object = B_Seurat_obj)
# Collect highly variable features
B_Seurat_obj <- FindVariableFeatures(object = B_Seurat_obj)
# Scaling
B_Seurat_obj <- ScaleData(object = B_Seurat_obj)
# PCA
B_Seurat_obj <- RunPCA(object = B_Seurat_obj)
elbowplot <- ElbowPlot(B_Seurat_obj)
ggsave("elbow_plot.pdf", plot = elbowplot, device = "pdf")
# Clustering
B_Seurat_obj <- FindNeighbors(B_Seurat_obj, dims = 1:30)
B_Seurat_obj <- FindClusters(B_Seurat_obj,resolution = 0.9)
# Non-Linear
B_Seurat_obj <- RunUMAP(B_Seurat_obj , dims = 1:30, n.neighbors = 50)
dimplot1 <- DimPlot(B_Seurat_obj, reduction = "umap", group.by = 'redefined_cluster', label = TRUE)
ggsave("cell_types_dim_plot.pdf", plot = dimplot1, device = "pdf", width = 16 , height = 10)
dimplot2 <- DimPlot(B_Seurat_obj, reduction = "umap", group.by = 'seurat_clusters', label = TRUE)
ggsave("clusters_dim_plot.pdf", plot = dimplot2, device = "pdf", width = 16 , height = 10)
cc <- dimplot1|dimplot2
ggsave("cellTypes_clusters_dim_plot.pdf", plot = cc, device = "pdf", width = 16 , height = 10)

# Monocle3
B_CDS <- as.cell_data_set(B_Seurat_obj) # Convert to cell_data_set object
B_CDS
str(B_CDS)
colData(B_CDS) # Get cell metadata
fData(B_CDS) # Get gene metadata
fData(B_CDS)$gene_short_name <- rownames(fData(B_CDS))
fData(B_CDS)
counts(B_CDS) # Get counts

# Assign partitions
B_CDS@clusters$UMAP$partitions
B_CDS@colData@rownames
partition <- c(rep(1,length(B_CDS@colData@rownames)))
names(partition) <- B_CDS@colData@rownames
partition <- as.factor(partition)
partition
B_CDS@clusters$UMAP$partitions <- partition

# Assign cluster info
B_Seurat_obj@active.ident
B_CDS@clusters$UMAP$clusters <- B_Seurat_obj@active.ident

# Assign UMAP info
B_CDS@int_colData@listData$reducedDims$UMAP <- B_Seurat_obj@reductions$umap@cell.embeddings


clusters <- plot_cells(B_CDS, color_cells_by = 'cluster' , label_groups_by_cluster = F , group_label_size = 5) +
  theme(legend.position = 'right')

celltypes <- plot_cells(B_CDS, color_cells_by = 'redefined_cluster' , label_groups_by_cluster = F , group_label_size = 5) +
  theme(legend.position = 'right')

cc <- clusters|celltypes
ggsave("cellTypes_clusters_dim_plot_2.pdf", plot = cc, device = "pdf", width = 16 , height = 10)
