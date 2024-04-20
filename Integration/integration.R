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
