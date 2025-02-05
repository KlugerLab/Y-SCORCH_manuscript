library(Seurat)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(cowplot)
options(future.globals.maxSize = 1e12)

seu.obj.msn <- readRDS("./data/msn_subset.RDS")


msn_ref <- readRDS("/banach2/SCORCH/data/analysis/resources/primate_nacc_vst_msn_monkeyP.RDS")

# for each sample, run against the reference
Idents(seu.obj.msn) <- "sample_ids"
seurat_query_list <- SplitObject(seu.obj.msn, split.by = "sample_ids")

prediction.df <- lapply(seurat_query_list, function(seurat_obj){
  print(unique(Idents(seurat_obj)))
  if(ncol(seurat_obj) <= 50){NULL}
  else{
    seurat_obj <- SCTransform(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)
    
    transfer_anchors_vst_ref <- FindTransferAnchors(reference = msn_ref, query = seurat_obj, 
                                                    normalization.method = "SCT",
                                                    reduction = "pcaproject")
    predictions_msn <- TransferData(anchorset = transfer_anchors_vst_ref, k.weight = 30,
                                    refdata = msn_ref$MSN_type)
    #seurat_obj <- AddMetaData(seurat_obj, metadata = predictions_msn)
    return(predictions_msn)
  }
  
})
saveRDS(prediction.df,"./data/msn_prediction.RDS")
