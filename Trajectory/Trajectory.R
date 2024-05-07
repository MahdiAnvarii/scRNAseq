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
