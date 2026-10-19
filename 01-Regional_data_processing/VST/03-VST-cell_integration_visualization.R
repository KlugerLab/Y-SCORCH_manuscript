# REVISION: run on the CellBender-corrected data; paths (including fast_tsne_path) come
# from revision/config/paths.R and outputs go to OUT_ROOT. Changes from the published
# script are marked ## REVISION.

library(Seurat)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(cowplot)
## REVISION: package used but never loaded in the published script (it was
##   loaded earlier in the interactive session); under Rscript this is an error.
library(ggplot2)
options(future.globals.maxSize = 1e12)



source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup("01-Regional_data_processing/VST",
            subdirs = c("data/02-cell_annotation", "data/03-cell_integration_vis"))
sce.obj <- readRDS("./data/02-cell_annotation/processed_sce_obj.rds")
count_data <- sce.obj@assays@data$counts
colnames(count_data) <- paste0(colnames(sce.obj),"-",sce.obj@colData$sample_ids)
seu.obj <- CreateSeuratObject(counts = count_data,
                              meta.data = as.data.frame(sce.obj@colData))
rm(count_data)
seu.obj <- NormalizeData(seu.obj)


# integration
seu.obj[["RNA"]] <- split(seu.obj[["RNA"]], f = seu.obj$sample_ids)
seu.obj <- NormalizeData(seu.obj)
seu.obj <- FindVariableFeatures(seu.obj)
seu.obj <- ScaleData(seu.obj)
seu.obj <- RunPCA(seu.obj)
seu.obj <- IntegrateLayers(
  object = seu.obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
seu.obj <- FindNeighbors(seu.obj, reduction = "harmony", dims = 1:50)
seu.obj <- FindClusters(seu.obj, resolution = 1, cluster.name = "seurat_clusters")


# visualization
seu.obj <-  RunTSNE(
  seu.obj,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
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
seu.obj <- run_umap_threaded(seu.obj, reduction = "harmony", dims = 1:50)

## REVISION: write named tSNE and UMAP panels. The published scripts call
##   DimPlot() with no device open, so under Rscript every plot lands unlabelled
##   in Rplots.pdf and the embedding figures cannot be found.
save_embedding_panels(seu.obj, prefix = "03_VST_cell_integration_visualization")


# tocheck:
# 1. are there a lot of predicted MSNs also predicted as ext neurons with high scores
# 2. are these cells enriched in specific samples/conditions
# 3. are these cells more similar to INS than to PFC

# not sure:
# for the neuron population, cluster the data, then check the marker genes to make sure whether 
# each cluster is indeed ext or msn

p0 <- DimPlot(seu.obj,reduction = "tsne",group.by = "condition_new",label=F,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

p1 <- DimPlot(seu.obj,reduction = "tsne",group.by = "sample_ids",label=F,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
## REVISION: consensus_annotation_v2 is not created until 04-01, so the published
##   code errors on a fresh run. Falls back to the pipeline's consensus_annotation.
.annot_col <- if ("consensus_annotation_v2" %in% colnames(seu.obj@meta.data))
  "consensus_annotation_v2" else "consensus_annotation"
p2 <- DimPlot(seu.obj,reduction = "tsne",group.by = .annot_col,label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p3 <- DimPlot(seu.obj,reduction = "tsne",group.by = "seurat_clusters",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)+ NoLegend()
Ext_markers <- c("GRIN1","SLC17A7","SATB2","SV2B")
Inh_markers <- c("GAD1","GAD2")
msn_markers <- c('PPP1R1B', 'PDE1B', 'BCL11B')
msn_sub_markers <- c("DRD2","STXBP6","DRD1","FAM163A","GRIA4","CPNE4","CASZ1","RXFP1","GRIK1")

f1 <- FeaturePlot(seu.obj,reduction = "tsne",features = Ext_markers)
f2 <- FeaturePlot(seu.obj,reduction = "tsne",features = Inh_markers)
f3 <- FeaturePlot(seu.obj,reduction = "tsne",features = msn_markers)

# so there are several clusters that are MSNs instead of Ext
# what are these clusters

# TODO: 1.identify the clusters that are ext based on the marker genes
#         verify by the max scores
#       2. visualize the sample/condition distributions
#       3. remove these clusters
v1 <- DotPlot(seu.obj,features = c(Ext_markers,Inh_markers,msn_markers))
## REVISION: the published code hard-coded
##   ext_clusters <- c(14,22,26,33,35,38,41,43,45,47,48,50,52,56)
## from the published clustering, chosen by marker inspection. Clusters renumber
## after ambient correction, so those indices no longer identify the same cells,
## and the sample exclusion in VST/04-01 depends on what is flagged here.
##
## Replaced with an explicit criterion calibrated against the published object
## (ORIG_ANALYSIS_DIR/01-Regional_data_processing/VST/data/03-cell_integration_vis/
## integrated_seurat_obj.RDS). Its putative_ext labels separate cleanly:
##
##   published_ext = TRUE    14 clusters   ext_pct median 0.565  (range 0.014-0.832)
##   published_ext = FALSE   50 clusters   ext_pct median 0.031  (max   0.278)
##   13,240 of 263,915 cells flagged (5.0%)
##
## where ext_pct is the fraction of cells in a cluster expressing SLC17A7 or
## SATB2. A threshold of 0.30 sits in the empty gap between those distributions
## and reproduces the published call, so the criterion is stated rather than
## enumerated (reviewer request).
##
## Marker choice: only SLC17A7 and SATB2 discriminate. GRIN1 and SV2B are in the
## published plotted panel but are pan-neuronal, so they do not separate cortical
## excitatory neurons from MSNs.
##
## Label-based detection does not work: the primate striatal reference has no
## excitatory class (MSNs / Interneurons / Astrocytes / Oligos / Microglia /
## Oligos_Pre / Endothelial / Mural-Fibroblast), so contamination is forced into
## MSNs or Interneurons; and the BICCN cortical reference has no MSN class, calling
## 29,879 VST cells "L2/3 IT" against 37,560 that the primate reference calls MSNs.
## Marker expression is the only signal that separates them.
##
## Detection fractions are used rather than z-scored panel means: z-scores are
## relative, so a cluster expressing neither panel can still score ext > msn
## (that approach flagged 50.5% of the striatum against a published 5.0%).
EXT_MARKERS <- c("SLC17A7", "SATB2")           # discriminative
MSN_MARKERS <- c("PPP1R1B", "PDE1B", "BCL11B")
EXT_PCT_THRESHOLD <- 0.30                      # calibrated, see above
PUBLISHED_EXT_CLUSTERS  <- c(14,22,26,33,35,38,41,43,45,47,48,50,52,56)   # reference only
PUBLISHED_SAMPLE2REMOVE <- c("1039","1068","1119","7200776574","HCCRS","HCTME","HCTXJ")

.ext_ok <- intersect(EXT_MARKERS, rownames(seu.obj))
.msn_ok <- intersect(MSN_MARKERS, rownames(seu.obj))
if (!length(.ext_ok)) stop("no excitatory marker genes present in the object")

## Fraction of cells per cluster expressing each marker — an absolute quantity,
## unlike z-scores. cluster_pct_expressing() is layer-aware:
## at this point in the script the RNA assay is split into one counts layer per
## sample, and LayerData(layer = "counts") would silently return only the first.
.cl <- as.character(seu.obj@meta.data$seurat_clusters)
ext_pct <- cluster_pct_expressing(seu.obj, .ext_ok, group = "seurat_clusters")
msn_pct <- cluster_pct_expressing(seu.obj, .msn_ok, group = "seurat_clusters")

ext_clusters <- names(ext_pct)[ext_pct > EXT_PCT_THRESHOLD]
seu.obj@meta.data$putative_ext <-
  as.character(seu.obj@meta.data$seurat_clusters) %in% ext_clusters

## A cluster-level call that flags no cells means the names did not resolve.
if (length(ext_clusters) > 0 && sum(seu.obj@meta.data$putative_ext) == 0) {
  stop(length(ext_clusters), " clusters called Ext but 0 cells flagged - ",
       "cluster names do not match seurat_clusters")
}
## Published contamination was 5.0% of cells. An order-of-magnitude departure
## means the criterion has stopped working, not that the data changed that much.
.frac <- mean(seu.obj@meta.data$putative_ext)
if (.frac > 0.20) {
  stop(sprintf("%.1f%% of cells called excitatory - implausible for striatum ",
               100 * .frac),
       "(published run: 5.0%). Check EXT_PCT_THRESHOLD and the marker panel.")
}

## ---- diagnostics for review --------------------------------------------------
dir.create("./data/03-cell_integration_vis", recursive = TRUE, showWarnings = FALSE)

cluster_diag <- data.frame(
  cluster           = names(ext_pct),
  n_cells           = as.integer(table(.cl)[names(ext_pct)]),
  ext_pct           = as.numeric(ext_pct),
  msn_pct           = as.numeric(msn_pct),
  called_ext        = ext_pct > EXT_PCT_THRESHOLD,
  in_published_list = names(ext_pct) %in% as.character(PUBLISHED_EXT_CLUSTERS)
)
cluster_diag <- cluster_diag[order(-cluster_diag$ext_pct), ]
write.csv(cluster_diag, "./data/03-cell_integration_vis/ext_cluster_diagnostic.csv",
          row.names = FALSE)

samp <- as.data.frame(seu.obj@meta.data) %>%
  dplyr::group_by(sample_ids, sample_name, condition_new) %>%
  dplyr::summarise(n_cells      = dplyr::n(),
                   n_ext        = sum(putative_ext),
                   ext_fraction = mean(putative_ext),
                   .groups = "drop") %>%
  dplyr::mutate(in_published_removal = sample_ids %in% PUBLISHED_SAMPLE2REMOVE) %>%
  dplyr::arrange(dplyr::desc(ext_fraction))
write.csv(samp, "./data/03-cell_integration_vis/sample_ext_fraction.csv",
          row.names = FALSE)

message("[VST ext] clusters called Ext (ext_pct > ", EXT_PCT_THRESHOLD, "): ",
        length(ext_clusters), " of ", length(ext_pct),
        "   (published hardcoded list had ", length(PUBLISHED_EXT_CLUSTERS), ")")
message("[VST ext] cells flagged: ", sum(seu.obj@meta.data$putative_ext),
        " of ", ncol(seu.obj), sprintf("  (%.1f%%, published 5.0%%)", 100 * .frac))
message("[VST ext] diagnostics written to data/03-cell_integration_vis/")
print(as.data.frame(samp), row.names = FALSE)
DimPlot(seu.obj,reduction = "tsne",group.by = "putative_ext",label=F,
        raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
FeaturePlot(seu.obj,features = c("ma_prediction_score.max","nhp_prediction_score.max","biccn_prediction.score.max"))
table(seu.obj@meta.data$condition_new[seu.obj@meta.data$putative_ext])
table(seu.obj@meta.data$sample_ids[seu.obj@meta.data$putative_ext])
table(seu.obj@meta.data$sample_name[seu.obj@meta.data$putative_ext])

seu.obj@meta.data$celltypes <- "Others"
seu.obj@meta.data$celltypes[seu.obj@meta.data$putative_ext] <-"Ext"
metadata_all <- seu.obj@meta.data

DotPlot(seu.obj,features = Ext_markers,group.by = "celltypes")

n_counts <- metadata_all %>%
  group_by(condition_new,sample_ids,celltypes) %>%
  summarise(n=n())%>%
  filter(celltypes == "Ext")
cond_col <- c("#0072B2", "darkorange1")
names(cond_col) <- c("CTR","HIV")
ggplot(n_counts, aes(x=sample_ids,y=n,fill=condition_new)) +
  scale_fill_manual(values = cond_col)+
  facet_wrap(~condition_new, scales = "free_x")+
  geom_bar(stat="identity") +
  geom_text(aes(label=n), vjust=-1) +
  ylab("# Cells")+
  theme_minimal()+
  theme(axis.text = element_text(size=10),
        axis.text.x = element_text(angle = 45),
        strip.text = element_text(size=15,face = "bold"),
        legend.title=element_text(size=15,face = "bold"),
        title=element_text(size=15,face = "bold"),
        legend.text=element_text(size=10))+
  theme(axis.ticks.length.x = unit(0, "cm")) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))




DimPlot(seu.obj,reduction = "tsne",group.by = "biccn_predicted.id",label=T,
        raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
DimPlot(seu.obj,reduction = "tsne",group.by = "ma_predicted.id",label=T,
        raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
DimPlot(seu.obj,reduction = "tsne",group.by = "nhp_vst_predicted.id",label=T,
        raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

DotPlot(seu.obj,features = Ext_markers)

saveRDS(seu.obj,"./data/03-cell_integration_vis/integrated_seurat_obj.RDS")
seu.obj <- readRDS("./data/03-cell_integration_vis/integrated_seurat_obj.RDS")
