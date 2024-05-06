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
