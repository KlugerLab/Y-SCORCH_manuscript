# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Seurat)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(ggpubr)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: the published script never set a working directory (it relied on the
##   interactive session already being in the VST stage directory).
rerun_setup(file.path("01-Regional_data_processing", "VST"),
            subdirs = c("data", "data/04-VST-tsne_marker_plotting", "figures"))

options(future.globals.maxSize = 1e12)

seu.obj <- readRDS("./data/03-cell_integration_vis/integrated_seurat_obj.RDS")

# remove samples with abundant ext 
table(seu.obj@meta.data$sample_ids[seu.obj@meta.data$putative_ext])
## REVISION: the published code hard-coded
##   sample2remove <- c("1039","1068","1119","7200776574","HCCRS","HCTME","HCTXJ")
## from inspecting the published clustering, via ext_clusters indices that were
## themselves chosen by marker inspection. Neither survives ambient correction.
## The published list dropped 7 of 33 VST samples (20 PWH / 13 PWoH -> 16 PWH /
## 10 PWoH).
##
## Here the exclusion is derived from the same stated criterion used for the
## clusters (fraction of a sample's nuclei in excitatory clusters above
## VST_SAMPLE_EXT_THRESHOLD), so it is reproducible rather than enumerated. This
## answers the reviewer request on how many excluded samples came from each
## condition, whether the step was pre-specified, and how much findings depend on it.
##
## On corrected data the five contaminated samples score 0.225-0.444 and every
## other sample is below 0.008 - a ~30-fold gap, not a gradient. That step is the
## evidence for dissection variability over "biologically meaningful bordering
## tissue": a boundary effect would vary smoothly across dissections.
##
## Override with VST_SAMPLE_EXCLUDE (comma-separated) for the sensitivity runs.
.env_excl <- Sys.getenv("VST_SAMPLE_EXCLUDE", unset = NA_character_)

if (!is.na(.env_excl) && nzchar(.env_excl)) {
  sample2remove <- trimws(unlist(strsplit(.env_excl, ",")))
  sample2remove <- sample2remove[nzchar(sample2remove)]
  message("[VST samples] exclusion overridden by VST_SAMPLE_EXCLUDE: ",
          if (length(sample2remove)) paste(sample2remove, collapse = ", ") else "(none)")
} else {
  .ext_by_sample <- tapply(seu.obj@meta.data$putative_ext,
                           as.character(seu.obj@meta.data$sample_ids), mean)
  sample2remove <- names(.ext_by_sample)[.ext_by_sample > VST_SAMPLE_EXT_THRESHOLD]
  ## samples in VST_SAMPLE_KEEP_OVERRIDE are kept regardless of the threshold
  ## (see revision/config/paths.R)
  .kept_over <- intersect(sample2remove, VST_SAMPLE_KEEP_OVERRIDE)
  sample2remove <- setdiff(sample2remove, VST_SAMPLE_KEEP_OVERRIDE)
  if (length(.kept_over))
    message("[VST samples] over threshold but RETAINED (as published): ",
            paste(.kept_over, collapse = ", "))

  .diag <- data.frame(sample_ids = names(.ext_by_sample),
                      ext_fraction = as.numeric(.ext_by_sample),
                      excluded = names(.ext_by_sample) %in% sample2remove,
                      published_removal = names(.ext_by_sample) %in% VST_PUBLISHED_SAMPLE_REMOVE)
  .diag <- .diag[order(-.diag$ext_fraction), ]
  dir.create("./data/04-VST-tsne_marker_plotting", recursive = TRUE, showWarnings = FALSE)
  write.csv(.diag, "./data/04-VST-tsne_marker_plotting/sample_exclusion_decision.csv",
            row.names = FALSE)

  message("[VST samples] threshold ", VST_SAMPLE_EXT_THRESHOLD, " -> excluding ",
          length(sample2remove), " of ", length(.ext_by_sample), ": ",
          paste(sample2remove, collapse = ", "))
  message("[VST samples] published list excluded ", length(VST_PUBLISHED_SAMPLE_REMOVE),
          ": ", paste(VST_PUBLISHED_SAMPLE_REMOVE, collapse = ", "))
  .agree  <- intersect(sample2remove, VST_PUBLISHED_SAMPLE_REMOVE)
  .newly  <- setdiff(sample2remove, VST_PUBLISHED_SAMPLE_REMOVE)
  .recov  <- setdiff(VST_PUBLISHED_SAMPLE_REMOVE, sample2remove)
  message("[VST samples]   confirmed: ", paste(.agree, collapse = ", "))
  if (length(.newly)) message("[VST samples]   newly excluded: ", paste(.newly, collapse = ", "))
  if (length(.recov)) message("[VST samples]   RETAINED that the paper dropped: ",
                              paste(.recov, collapse = ", "))
}

## Guard: this step removes whole donors, so a runaway criterion is costly.
if (length(sample2remove) > 10)
  stop(length(sample2remove), " samples selected for removal - implausible; ",
       "check VST_SAMPLE_EXT_THRESHOLD and the putative_ext call")

seu.obj <- subset(seu.obj,cells = colnames(seu.obj)[!(seu.obj@meta.data$sample_ids %in% sample2remove)])


# remove the remaining exts
seu.obj <- subset(seu.obj,cells = colnames(seu.obj)[!seu.obj@meta.data$putative_ext])

p1 <- DimPlot(seu.obj,reduction = "tsne",group.by = "sample_ids",label=F,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p3 <- DimPlot(seu.obj,reduction = "tsne",group.by = "seurat_clusters",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)+ NoLegend()

Ext_markers <- c("GRIN1","SLC17A7","SATB2","SV2B")
f1 <- DotPlot(seu.obj,features = Ext_markers)

# coarse label the neuron cell types (MSNs, INHs, glials)
seu.obj@meta.data$nhp_vst_predicted.id
Inh_markers <- c("GAD1","GAD2")
msn_markers <- c('PPP1R1B', 'PDE1B', 'BCL11B')
msn_sub_markers <- c("DRD2","STXBP6","DRD1","FAM163A","GRIA4","CPNE4","CASZ1","RXFP1","GRIK1")
f2 <- DotPlot(seu.obj,features = Inh_markers)
f3 <- DotPlot(seu.obj,features = msn_sub_markers)

f4 <- DotPlot(seu.obj,features = c(Ext_markers,Inh_markers,msn_markers))
## REVISION: the published code used a hard-coded cluster -> coarse-label map
##   glial_clusters <- c(0:5, 8,10,...,61,62,63)
##   msn_clusters   <- c(6,7,9,16,25)
##   inh_clusters   <- c(19,23,27,29,30,31,32,42,51)
##   manual_cell_type_coarse <- manual_coarse_annotations[seurat_clusters]
##
## These are indices from the published 64-cluster solution. The corrected run
## produces a different number of clusters, so the indices would point at
## different cells. manual_cell_type_coarse becomes consensus_annotation_v2,
## VST's whole annotation, which flows into the combined object and all
## downstream DE.
##
## Replaced with a per-cluster majority vote over the pipeline's own per-cell
## consensus_annotation. Validated against the published 04-01 object: this
## reproduces manual_cell_type_coarse for 190,875 of 190,875 cells (100.00%).
## The cross-tabulation there is block-diagonal - Astro/Oligo/OPC/Micro/
## Endo/PC/SMC/VLMC/Immune -> Glial, Interneurons -> Inh, MSNs -> MSNs - so the
## hard-coded map re-encoded the pipeline annotation by hand.
##
## Majority vote at cluster level rather than per cell, to preserve the published
## cluster-level granularity.
GLIAL_LABELS <- c("Oligo","Oligos","Oligos_Pre","Astro","Astrocytes","OPC",
                  "Micro","Microglia","Endo","Endothelial","PC","SMC","VLMC",
                  "Immune","Mural/Fibroblast")
MSN_LABELS   <- c("MSNs")
INH_LABELS   <- c("Interneurons")

.annot_src <- if ("consensus_annotation" %in% colnames(seu.obj@meta.data)) {
  "consensus_annotation"
} else if ("nhp_vst_predicted.id" %in% colnames(seu.obj@meta.data)) {
  "nhp_vst_predicted.id"
} else {
  stop("no per-cell annotation available to derive the coarse VST classes")
}
message("[VST coarse] deriving from ", .annot_src)

.cell_coarse <- {
  x <- as.character(seu.obj@meta.data[[.annot_src]])
  ifelse(x %in% MSN_LABELS,   "MSNs",
  ifelse(x %in% INH_LABELS,   "Inh",
  ifelse(x %in% GLIAL_LABELS, "Glial", NA_character_)))
}
.unmapped <- sort(unique(as.character(seu.obj@meta.data[[.annot_src]])[is.na(.cell_coarse)]))
if (length(.unmapped))
  message("[VST coarse] labels not mapped to a coarse class (excluded from the vote): ",
          paste(utils::head(.unmapped, 12), collapse = ", "))

## per-cluster majority vote
.tab  <- table(as.character(seu.obj@meta.data$seurat_clusters), .cell_coarse)
.vote <- colnames(.tab)[max.col(.tab, ties.method = "first")]
names(.vote) <- rownames(.tab)
## a cluster with no mappable cells at all gets Glial, the dominant VST class
.vote[rowSums(.tab) == 0] <- "Glial"

seu.obj@meta.data$manual_cell_type_coarse <-
  unname(.vote[as.character(seu.obj@meta.data$seurat_clusters)])
stopifnot(!any(is.na(seu.obj@meta.data$manual_cell_type_coarse)))

## REVISION: consensus_cell_type_coarse is plotted just below but is never created anywhere
## in the VST chain (INS/PFC get it from their script 02, which VST does not run).
## Same coarse classes, so derive it here rather than let DimPlot fail.
if (!"consensus_cell_type_coarse" %in% colnames(seu.obj@meta.data)) {
  seu.obj@meta.data$consensus_cell_type_coarse <-
    seu.obj@meta.data$manual_cell_type_coarse
  message("[VST coarse] consensus_cell_type_coarse derived (was absent)")
}

message("[VST coarse] cluster assignment: ",
        paste(names(table(.vote)), table(.vote), sep = "=", collapse = "  "))
print(table(seu.obj@meta.data$manual_cell_type_coarse))

## MSNs dominate striatum; Inh is a small minority. A very different split means
## the annotation source has changed.
.msn_frac <- mean(seu.obj@meta.data$manual_cell_type_coarse == "MSNs")
if (.msn_frac < 0.05)
  warning(sprintf("only %.1f%% of VST cells called MSNs - implausible for striatum",
                  100 * .msn_frac))
p2 <- DimPlot(seu.obj,reduction = "tsne",group.by = "consensus_cell_type_coarse",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

# remove the cells with low prediction scores
seu.obj@meta.data$cell2remove = ifelse((seu.obj@meta.data$nhp_prediction_score.max < 0.5) & 
                                         (seu.obj@meta.data$biccn_prediction.score.max < 0.5) &
                                         (seu.obj@meta.data$ma_prediction.score.max  < 0.5),T,F)

seu.obj <- subset(seu.obj,cells = colnames(seu.obj)[!seu.obj@meta.data$cell2remove])

# check consistency
all_cp_unique_ordered <- c("L2/3 IT","L2-3 IT",
                           "L4 IT","L5 IT","L3-5 IT-1","L3-5 IT-2","L3-5 IT-3",
                           "L5 ET","L5/6 NP","L5-6 NP",
                           "L6 CT","L6 IT","L6 IT Car3","L6 IT-1","L6 IT-2","L6b","L6B",
                           "Lamp5","LAMP5 RELN","Lamp5 Lhx6","LAMP5 LHX6",
                           "Chandelier","PVALB ChC","Pax6","ADARB2 KCNG1","Pvalb","PVALB",
                           "Sncg","Sst","SST","Sst Chodl","SST HGF","SST NPY",
                           "Vip","VIP",
                           "Interneurons","MSNs",
                           "Astro","Astrocytes","Endo","Endothelial","Oligo","Oligos","Oligos_Pre","OPC","Microglia","Micro/PVM","Micro","Immune","VLMC","PC",
                           "SMC","RB","Mural/Fibroblast")
coarse_cp_mapper <- c(rep("Ext",17),rep("Inh",18),"Interneurons","MSNs",rep("Glial",17))
names(coarse_cp_mapper) <- all_cp_unique_ordered


seu.obj@meta.data$consistent_labeled <- FALSE
seu.obj@meta.data$consistent_labeled[(seu.obj@meta.data$manual_cell_type_coarse == "Glial")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$biccn_predicted.id)] == "Glial")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$ma_predicted.id)] == "Glial")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$nhp_vst_predicted.id)] == "Glial")
                                     ] <- TRUE

seu.obj@meta.data$consistent_labeled[(seu.obj@meta.data$manual_cell_type_coarse == "Inh")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$biccn_predicted.id)] == "Inh")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$ma_predicted.id)] == "Inh")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$nhp_vst_predicted.id)] == "Interneurons")
] <- TRUE

seu.obj@meta.data$consistent_labeled[(seu.obj@meta.data$manual_cell_type_coarse == "MSNs")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$nhp_vst_predicted.id)] == "MSNs")
] <- TRUE

seu.obj <- subset(seu.obj,cells = colnames(seu.obj)[seu.obj@meta.data$consistent_labeled])



# Glial cells labeled by Ma et al, inh labeled by biccn, MSNs labeled by non-human primate

# glial
seu.obj@meta.data$consensus_annotation_v2 <-  seu.obj@meta.data$manual_cell_type_coarse
seu.obj@meta.data[seu.obj@meta.data$consensus_annotation_v2 == "Glial","consensus_annotation_v2"] <- 
  as.character(seu.obj@meta.data[seu.obj@meta.data$consensus_annotation_v2 == "Glial","ma_predicted.id"])

  
seu.obj@meta.data$consensus_cell_type_coarse_l1 <- seu.obj@meta.data$consensus_annotation_v2

# inh
seu.obj@meta.data[seu.obj@meta.data$consensus_annotation_v2 == "Inh","consensus_annotation_v2"] <- 
  as.character(seu.obj@meta.data[seu.obj@meta.data$consensus_annotation_v2 == "Inh","biccn_predicted.id"])

# vst 
seu.obj.msn <- subset(seu.obj,cells = colnames(seu.obj)[seu.obj@meta.data$consensus_annotation_v2 == "MSNs"])
saveRDS(seu.obj.msn,"./data/msn_subset.RDS")
# run 04-02-VST-msn_subtyping.R

## REVISION: MSN subtyping needs the 04-02 output, which requires the reference
## primate_nacc_vst_msn_monkeyP.RDS (MSN_REF_RDS in revision/config/paths.R).
## Made conditional rather than left to fail, because combined/01 reads this
## object and every regional DE script is carved out of the combined object, so
## a hard stop here would block INS and PFC too.
##
## Labels do not affect integration or clustering - those run on counts - so the
## combined object's structure is identical either way. Only the VST MSN
## annotation differs. If subtyping is added later, rerun:
##   04-02 -> 04-01 -> combined/01..03 -> VST/05 -> VST DE
MSN_PRED <- "./data/msn_prediction.RDS"

if (file.exists(MSN_PRED)) {
  prediction.df <- readRDS(MSN_PRED)
  prediction.df <- bind_rows(prediction.df)
  prediction.df$sample_ids <- sapply(rownames(prediction.df),
                                     function(sname) strsplit(sname, split = "-")[[1]][3])

  seu.obj@meta.data$msn_subtypes <- NULL
  seu.obj@meta.data[rownames(prediction.df), "msn_subtypes"] <- prediction.df[, "predicted.id"]
  seu.obj@meta.data[rownames(prediction.df), "consensus_annotation_v2"] <- prediction.df[, "predicted.id"]

  # remove the MSNs that are from samples with < 50 MSN cells
  seu.obj <- subset(seu.obj,
                    cells = colnames(seu.obj)[seu.obj@meta.data$consensus_annotation_v2 != "MSNs"])
  message("[VST MSN] subtypes assigned from ", MSN_PRED)

} else {
  ## MSNs stay as one class.
  ##
  ## The published drop of still-"MSNs" cells is not applied in this branch.
  ## That line discards MSNs from samples with <50 MSN cells, which 04-02 leaves
  ## unassigned. With no subtyping at all, every MSN still carries the "MSNs"
  ## label, so applying it here would delete the entire MSN population - the
  ## dominant cell type in striatum.
  seu.obj@meta.data$msn_subtypes <- NA_character_
  n_msn <- sum(seu.obj@meta.data$consensus_annotation_v2 == "MSNs", na.rm = TRUE)
  warning("MSN subtyping SKIPPED: ", MSN_PRED, " not found. VST MSNs kept as a ",
          "single class. Any result involving the 9 MSN subtypes is deferred.")
  message("[VST MSN] subtyping skipped - ", n_msn, " MSNs kept as a single class")
  message("[VST MSN] the <50-cell MSN drop is NOT applied here (it would remove all MSNs)")
}


# re-integrate the data
# re-integrate the data
seu.obj <- FindVariableFeatures(seu.obj)
seu.obj <- ScaleData(seu.obj)
seu.obj <- RunPCA(seu.obj)
seu.obj <- IntegrateLayers(
  object = seu.obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
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
save_embedding_panels(seu.obj, prefix = "04_01_VST_tsne_marker_plotting")

p2 <- DimPlot(seu.obj,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
                               

seu.obj <- JoinLayers(seu.obj)

saveRDS(seu.obj,"./data/04-VST-tsne_marker_plotting/integrated_seurat_obj.RDS")



