# scRNAseq
 Data for "Pipeline" folder:  
 https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_DTC_Melanoma_3p_gemx_Multiplex/10k_Human_DTC_Melanoma_3p_gemx_Multiplex_count_raw_feature_bc_matrix.h5  
   
 Data for "Integration" folder:  
 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180665  

 Data for "DoubletFinder" folder:  
 https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.tar.gz  
 
 Data for "FindMarkers" folder:  
 library(SeuratData)  
 install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source")  
 library(ifnb.SeuratData)  
 InstallData("ifnb")  
 data("ifnb")  

 Data for "PseudoBulk" folder:  
 library(ExperimentHub)  
 experimenthub = ExperimentHub(localHub=FALSE)  
 query(experimenthub,"Kang")  
 sce_obj <- experimenthub[['EH2259']]  

 Data for "Trajectory" folder:  
 http://scrna.sklehabc.com/session/70521e3f05de2335e17476d2de37e179/download/download.umi?w=  
 http://scrna.sklehabc.com/session/70521e3f05de2335e17476d2de37e179/download/download.marker?w=  
 http://scrna.sklehabc.com/session/70521e3f05de2335e17476d2de37e179/download/download.meta?w=  
