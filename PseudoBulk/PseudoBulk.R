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

### Preprocess workflow
# Normalization
Seurat_obj <- NormalizeData(object = Seurat_obj)
# Collect highly variable features
Seurat_obj <- FindVariableFeatures(object = Seurat_obj)
# Scaling
Seurat_obj <- ScaleData(object = Seurat_obj)
# PCA
Seurat_obj <- RunPCA(object = Seurat_obj)
elbowplot <- ElbowPlot(Seurat_obj)
ggsave("elbow_plot.pdf", plot = elbowplot, device = "pdf")
# Non-Linear
Seurat_obj <- RunUMAP(Seurat_obj , dims = 1:20)
dimplot1 <- DimPlot(Seurat_obj, reduction = "umap", group.by = 'cell', label = TRUE)
ggsave("cell_types_dim_plot.pdf", plot = dimplot1, device = "pdf", width = 16 , height = 10)
dimplot2 <- DimPlot(Seurat_obj, reduction = "umap", group.by = 'stim')
ggsave("conditions_dim_plot.pdf", plot = dimplot2, device = "pdf", width = 16 , height = 10)
cc <- dimplot1|dimplot2
ggsave("conditions_clusters_dim_plot.pdf", plot = cc, device = "pdf", width = 16 , height = 10)

#PseudoBulk workflow
view(Seurat_obj@meta.data)
Seurat_obj$sample <- paste0(Seurat_obj$ind, "_" , Seurat_obj$stim)
view(Seurat_obj@meta.data)
DefaultAssay(Seurat_obj)

#count aggregate to sample level
cts <- AggregateExpression(Seurat_obj,
                    group.by = c('cell' , 'sample'),
                    assays = "originalexp",
                    slot = "counts",
                    return.seurat = FALSE)
cts <- cts$originalexp
view(cts)

cts.t <- t(cts)
cts.t <- as.data.frame(cts.t)
cts.t[1:10,1:10]
SplitRows <- gsub('_.*','',rownames(cts.t))
SplitRows

# Split data frame
cts.split <- split.data.frame(cts.t, f = factor(SplitRows))
cts.split$`B cells`[1:10,1:10]

# Fix columns and retranspose
gsub('.*_(.*)','\\1','B cells_101-ctrl')
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)','\\1',rownames(x))
  t(x)
})
cts.split.modified$`B cells`[1:10,1:10]

# Run DESeq2 for B cells
BcellsCountMatrix <- cts.split.modified$`B cells`
view(BcellsCountMatrix)
BcellsMetadata <- data.frame(samples = colnames(BcellsCountMatrix))
view(BcellsMetadata)
BcellsMetadata <- BcellsMetadata %>%
  mutate(condition = ifelse(grepl('stim', samples), 'stim' , 'ctrl')) %>%
  column_to_rownames(var = 'samples')
view(BcellsMetadata)

dds <- DESeqDataSetFromMatrix(countData = BcellsCountMatrix,
                       colData = BcellsMetadata,
                       design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "condition_stim_vs_ctrl")
res
