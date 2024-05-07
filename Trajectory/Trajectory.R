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
