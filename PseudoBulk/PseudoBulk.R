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
