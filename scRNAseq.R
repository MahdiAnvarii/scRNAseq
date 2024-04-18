library(Seurat)
library(tidyverse)

Human_DTC <- Read10X_h5(filename = '../Data/10k_Human_DTC_Melanoma_3p_gemx_Multiplex_count_raw_feature_bc_matrix.h5')
str(Human_DTC)
nrow(Human_DTC)
ncol(Human_DTC)