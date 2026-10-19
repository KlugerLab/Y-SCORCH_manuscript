# combined visualization
#
# REVISION: run on the CellBender-corrected data; paths (including fast_tsne_path) come
# from revision/config/paths.R and outputs go to OUT_ROOT. Changes from the published
# script are marked ## REVISION.

library(Seurat)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup("01-Regional_data_processing",
            subdirs = c("INS/data/03-cell_integration_vis", "PFC/data/03-cell_integration_vis", "VST/data/04-VST-tsne_marker_plotting", "combined/data"))
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
  fast_tsne_path=FAST_TSNE_PATH, ann_not_vptree=FALSE, nthreads=12
)

## REVISION (reviewer request): UMAP alongside the tSNE above.
##   Same harmony reduction and dims, so the two embeddings are directly
##   comparable. Stored on the object before it is saved. Not part of the
##   published analysis.
##
##   Threads: Seurat passes n_threads = future::nbrOfWorkers() to uwot, so a
##   direct n_threads argument to RunUMAP is silently dropped. run_umap_threaded()
##   sets the future plan for the call and restores it afterwards.
seurat.combined <- run_umap_threaded(seurat.combined, reduction = "harmony", dims = 1:50)

## REVISION: write named tSNE and UMAP panels. The published scripts call
##   DimPlot() with no device open, so under Rscript every plot lands unlabelled
##   in Rplots.pdf and the embedding figures cannot be found.
save_embedding_panels(seurat.combined, prefix = "01_integration_tsne")

saveRDS(seurat.combined,"./combined/data/seurat_combined.RDS")



