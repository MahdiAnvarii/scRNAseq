library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

folders <- list.dirs(path = '../Data/', recursive = FALSE, full.names = FALSE)
folders
for(i in folders){
  cts <- ReadMtx(mtx = paste0('../Data/',i,'/matrix.mtx.gz'),
                 features = paste0('../Data/',i,'/features.tsv.gz'),
                 cells = paste0('../Data/',i,'/barcodes.tsv.gz'))
  
  name <- gsub('_filtered_feature_bc_matrix','',i)
  assign(name, CreateSeuratObject(counts = cts))
}

# Merge datassets
merged_seurat_obj <- merge(HB17_background, y=c(HB17_PDX,HB17_tumor,HB30_PDX,HB30_tumor,HB53_background,HB53_tumor),
      add.cell.ids = ls()[3:9],
      project = "HB")
merged_seurat_obj
view(merged_seurat_obj@meta.data)
merged_seurat_obj$sample <- rownames(merged_seurat_obj@meta.data)
merged_seurat_obj@meta.data <- separate(merged_seurat_obj@meta.data, col = 'sample', into = c('patient','type','barcode'), sep = '_')
view(merged_seurat_obj@meta.data)
unique(merged_seurat_obj@meta.data$patient)
unique(merged_seurat_obj@meta.data$type)

# QC
merged_seurat_obj[["MT_percent"]] <- PercentageFeatureSet(merged_seurat_obj , patter = "^MT-")
view(merged_seurat_obj@meta.data)










