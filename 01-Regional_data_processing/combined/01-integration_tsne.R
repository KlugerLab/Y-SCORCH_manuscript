# combined visualization

library(Seurat)

wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/01-Regional_data_processing/"
setwd(wd)

seurat.ins <- readRDS("./INS/data/03-cell_integration_vis/integrated_seurat_obj.RDS")
seurat.pfc <- readRDS("./PFC/data/03-cell_integration_vis/integrated_seurat_obj.RDS")
seurat.vst <- readRDS("./VST/data/04-VST-tsne_marker_plotting/integrated_seurat_obj.RDS")

count_data <- cbind(seurat.ins@assays$RNA@layers$counts,
                    seurat.pfc@assays$RNA@layers$counts,
                    seurat.vst@assays$RNA@layers$counts
)
cols2keep <- c("sample_name","sample_ids","region","condition_new",
               "consensus_cell_type_coarse","consensus_cell_type_coarse_l1","consensus_annotation_v2")
metadata <- rbind(seurat.ins@meta.data[,cols2keep],
                  seurat.pfc@meta.data[,cols2keep],
                  seurat.vst@meta.data[,cols2keep])

new_barcode_names <- c(paste0(colnames(seurat.ins),"-INS"),
                       paste0(colnames(seurat.pfc),"-PFC"),
                       paste0(colnames(seurat.vst),"-VST")
)
rownames(metadata) <- new_barcode_names
colnames(count_data) <- new_barcode_names
rownames(count_data) <- rownames(seurat.ins) 

seurat.combined <- CreateSeuratObject(counts = count_data,
                                      meta.data = metadata)

# integrate the data

seurat.combined[["RNA"]] <- split(seurat.combined[["RNA"]], f = seurat.combined$sample_name)
seurat.combined <- NormalizeData(seurat.combined)
seurat.combined <- FindVariableFeatures(seurat.combined)
seurat.combined <- ScaleData(seurat.combined)
seurat.combined <- RunPCA(seurat.combined)
seurat.combined <- IntegrateLayers(
  object = seurat.combined, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# visualization
seurat.combined <-  RunTSNE(
  seurat.combined,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path="/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)

saveRDS(seurat.combined,"./combined/data/seurat_combined.RDS")



