# ==============================================================================
# "Other immune" subset: re-embedded tSNE and marker feature plots (supplemental figure).
#
# Two versions, selected with --labels:
#   Immune              2,226 cells - the label the figures show as "Other immune"
#   Immune,Macrophages  4,399 cells - with the myeloid compartment
# Micro (36,677 cells) is excluded from both; it is an order of magnitude larger
# and would dominate any shared embedding.
#
# The marker set spans lymphoid (CD3D/E/G, CD8A, SKAP1, THEMIS) and myeloid
# (CD14, FCGR3A, CD83) genes, so sparse myeloid panels in the Immune-only
# version are a result about that population, not a defect.
#
# Method: the combined embedding was computed over 641,718 cells of every type,
# so 2,226 immune nuclei occupy a small blob in it with no visible internal
# structure. The subset is therefore re-embedded on its own: Normalize ->
# VariableFeatures -> Scale -> PCA on the joined layer, Harmony on sample_name,
# FIt-SNE on the Harmony reduction (same batch variable as the astrocyte
# analysis). The subset is not clustered or annotated.
#
# The Immune label is skewed toward PWH - 1,869 of 2,226 cells (84%) against a
# cohort that is ~60% PWH (INS 878/74, PFC 383/61, VST 608/222) - so structure
# in the embedding partly tracks condition rather than cell identity.
#
# Input:  OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3_no1040_msn.RDS
# Output: OUT_ROOT/07-immune_analysis/data/immune_subset_reembedded_<tag>.RDS
#         OUT_ROOT/07-immune_analysis/figures/immune_feature_grid_<tag>.{png,pdf},
#         immune_tsne_by_<variable>_<tag>.png
#
# Usage:
#   revision/docker/run.sh revision/tools/immune_subset_features.R            # Immune
#   revision/docker/run.sh revision/tools/immune_subset_features.R --labels Immune,Macrophages
#   add --replot to redraw figures from the saved embedding
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(ggplot2); library(patchwork)
  library(dplyr)
})

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_cluster_features.R"))  # feature_scale()

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }

## --labels takes a comma-separated set of consensus_annotation_v3 labels.
LABELS <- strsplit(getarg("--labels", "Immune"), ",")[[1]]
DIMS   <- as.integer(getarg("--dims", "30"))
TAG    <- getarg("--tag", paste(gsub("[^A-Za-z]", "", LABELS), collapse = "_"))

## Marker set. CX3CR1 is used rather than CXCR1, which is detected in 0 of 2,226
## cells (CX3CR1 in 117). Any gene with zero detection is still drawn grey and
## titled "(not detected)" rather than silently dropped.
FEATURES <- c("CD8A","CD3D","CD3E","CD3G","SKAP1","THEMIS",
              "CD14","CX3CR1","FCGR3A","CD83","ELL2","GLUL")

rerun_setup("07-immune_analysis", subdirs = c("data", "figures"))

IN <- file.path(OUT_ROOT, "01-Regional_data_processing", "combined", "data",
                "seurat_combined_v3_no1040_msn.RDS")
if (!file.exists(IN)) {
  alt <- file.path(OUT_ROOT, "01-Regional_data_processing", "combined", "data",
                   "seurat_combined_v3_no1040.RDS")
  if (file.exists(alt)) { IN <- alt
    warning("using ", basename(alt), " - no MSN subtypes, harmless here since ",
            "this subset contains no MSNs", immediate. = TRUE)
  } else stop("no combined object found")
}
message("[immune] input: ", basename(IN))

## --replot redraws the figures from the saved embedding without redoing the
## subset/Harmony/tSNE, which is the expensive part.
REPLOT <- "--replot" %in% args
SAVED  <- sprintf("./data/immune_subset_reembedded_%s.RDS", TAG)
if (REPLOT) {
if (!file.exists(SAVED)) stop("--replot needs ", SAVED, " - run without --replot first")
message("[immune] --replot: loading saved embedding ", SAVED)
imm <- readRDS(SAVED)
} else {
obj <- readRDS(IN)
Idents(obj) <- "consensus_annotation_v3"
present <- unique(as.character(obj@meta.data$consensus_annotation_v3))
miss_lab <- setdiff(LABELS, present)
if (length(miss_lab))
  stop("label(s) not present: ", paste(miss_lab, collapse = ", "),
       "\nAvailable: ", paste(sort(present), collapse = ", "))

imm <- subset(obj, idents = LABELS)
rm(obj); gc(verbose = FALSE)

md <- imm@meta.data
message("[immune] ", ncol(imm), " cells across label(s): ", paste(LABELS, collapse = ", "))
if (length(LABELS) > 1) { cat("\nby label:\n"); print(table(as.character(md$consensus_annotation_v3))) }
cat("\nby region:\n");    print(table(as.character(md$region)))
cat("\nby condition:\n"); print(table(as.character(md$region), as.character(md$condition_new)))
cat("\ndonors:", length(unique(as.character(md$sample_ids))), "\n")

## ---- re-embed ----------------------------------------------------------------
## Deviation from the astrocyte workflow: the astrocyte notebook splits the RNA
## assay into one layer per sample_name and lets IntegrateLayers() treat layers
## as batches. With 2,226 cells over 86 region-sample layers, 38 layers have
## under 10 cells and several have one, and FindVariableFeatures' per-layer
## loess fails ("simpleLoess ... invalid 'x'"). So normalisation, variable
## features, scaling and PCA run on the joined layer, and Harmony corrects the
## PCA directly with sample_name as a metadata column.
imm <- JoinLayers(imm)
imm <- NormalizeData(imm)
imm <- FindVariableFeatures(imm)
imm <- ScaleData(imm)
imm <- RunPCA(imm)

npc <- min(DIMS, ncol(Embeddings(imm, "pca")))
if (npc < DIMS) message("[immune] only ", npc, " PCs available; using those")

BATCH <- getarg("--batch", "sample_name")
message("[immune] Harmony batch variable: ", BATCH, " (",
        length(unique(as.character(imm@meta.data[[BATCH]]))), " levels)")
imm <- harmony::RunHarmony(imm, group.by.vars = BATCH, reduction.use = "pca",
                           dims.use = 1:npc, reduction.save = "harmony_immune")

imm <- RunTSNE(imm, reduction = "harmony_immune", dims = 1:npc,
               tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE,
               seed.use = 3, fast_tsne_path = FAST_TSNE_PATH,
               ann_not_vptree = FALSE, nthreads = 16,
               reduction.name = "tsne_immune")

imm <- JoinLayers(imm)
saveRDS(imm, sprintf("./data/immune_subset_reembedded_%s.RDS", TAG))
message("[immune] wrote ./data/immune_subset_reembedded_", TAG, ".RDS")
}  ## end of the subset + re-embed block skipped by --replot

## ---- feature plots -----------------------------------------------------------
have    <- FEATURES[FEATURES %in% rownames(imm)]
missing <- setdiff(FEATURES, have)
if (length(missing))
  message("[immune] NOT in the matrix: ", paste(missing, collapse = ", "))
## detected at all? a gene present but all-zero plots as an empty panel
det <- vapply(have, function(g) sum(GetAssayData(imm, layer = "data")[g, ] > 0), integer(1))
message("[immune] cells expressing each feature:")
print(sort(det, decreasing = TRUE))
if (any(det == 0))
  message("[immune] zero-detection (panels will be blank): ",
          paste(names(det)[det == 0], collapse = ", "))

fp <- FeaturePlot(imm, features = have, reduction = "tsne_immune",
                  order = TRUE, combine = FALSE)
names(fp) <- have
fp <- lapply(have, function(g) {
  ## Each panel keeps its OWN expression legend: genes differ in range, so a
  ## shared or absent key would make a dim panel look like low expression when it
  ## may just be a different scale. Gene names large and bold for the figure.
  p <- fp[[g]] + theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text  = element_text(size = 10),
          legend.key.height = unit(0.9, "lines"),
          legend.key.width  = unit(0.45, "lines"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 22,
                                    colour = "black", margin = margin(b = 6)))
  if (det[[g]] == 0) {
    ## A gene detected in NO cell has a single-valued colour scale, which ggplot
    ## maps to the middle of the ramp, so the panel would look uniformly
    ## expressed. Force the floor colour and say so in the title.
    suppressMessages(p + scale_colour_gradient(low = "grey85", high = "grey85",
                                               limits = c(0, 1)) +
                       ggtitle(paste0(g, " (not detected)")) +
                       ## no legend: there is no expression range to key
                       theme(legend.position = "none"))
  } else {
    suppressMessages(p + feature_scale("grey_red", name = "Expression"))
  }
})

## 3 columns x 4 rows, portrait
NCOL <- 3
grid <- wrap_plots(fp, ncol = NCOL)
ggsave(sprintf("./figures/immune_feature_grid_%s.png", TAG), grid, width = 15,
       height = 4 * ceiling(length(fp) / NCOL), dpi = 200)
ggsave(sprintf("./figures/immune_feature_grid_%s.pdf", TAG), grid, width = 15,
       height = 4 * ceiling(length(fp) / NCOL))
message("[immune] wrote ./figures/immune_feature_grid_", TAG, ".{png,pdf}")

## the embedding itself, coloured by region and by condition - not a labelled
## annotation, just orientation for the reader
for (v in c(if (length(LABELS) > 1) "consensus_annotation_v3", "region", "condition_new")) {
  p <- DimPlot(imm, reduction = "tsne_immune", group.by = v, shuffle = TRUE) +
    theme_void() + ggtitle(v)
  ggsave(sprintf("./figures/immune_tsne_by_%s_%s.png", v, TAG), p, width = 7, height = 6, dpi = 200)
}
message("[immune] wrote the two orientation tSNEs")
message("\n[immune] done - ", file.path(OUT_ROOT, "07-immune_analysis"))
