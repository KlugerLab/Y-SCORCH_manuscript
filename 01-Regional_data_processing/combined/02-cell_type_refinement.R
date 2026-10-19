# for the glial cells
# remove Endo/Peri
# remove ambiguous clusters
#
# REVISION: run on the CellBender-corrected data; paths (including fast_tsne_path) come
# from revision/config/paths.R and outputs go to OUT_ROOT. Changes from the published
# script are marked ## REVISION.

library(Seurat)
## REVISION: package used but never loaded in the published script (it was
##   loaded earlier in the interactive session); under Rscript this is an error.
library(dplyr)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup("01-Regional_data_processing",
            subdirs = c("combined/data"))
seurat.combined <- readRDS("./combined/data/seurat_combined.RDS")

## REVISION: checkpoint resume.
##   Everything from the subset() down to RunTSNE() below is ~1 hour of compute
##   that depends only on seurat_combined.RDS, and the block further down saves
##   its result to seurat_glial_clustered.RDS. When that checkpoint exists, load
##   it instead of recomputing an identical object. Set SCORCH_NO_RESUME=1 to
##   force the full path.
##
##   The checkpoint is written by the same code this branch skips. Delete the
##   file to invalidate it after editing anything in the block below.
GLIAL_CKPT <- "./combined/data/seurat_glial_clustered.RDS"
RESUME_GLIAL <- file.exists(GLIAL_CKPT) && !nzchar(Sys.getenv("SCORCH_NO_RESUME"))

if (RESUME_GLIAL) {
  message("[combined/02] RESUMING from ", GLIAL_CKPT,
          " - skipping subset/integration/clustering/tSNE.",
          "  (SCORCH_NO_RESUME=1 to recompute)")
  seurat.glial <- readRDS(GLIAL_CKPT)
  stopifnot(!is.null(seurat.glial@meta.data$seurat_clusters),
            "tsne" %in% names(seurat.glial@reductions))
  message("[combined/02] checkpoint: ", ncol(seurat.glial), " glial cells, ",
          length(unique(as.character(seurat.glial@meta.data$seurat_clusters))), " clusters")
} else {

seurat.glial <- subset(seurat.combined,
                       cells = colnames(seurat.combined)[!(seurat.combined@meta.data$consensus_cell_type_coarse_l1 %in% c("Ext","Inh","MSNs"))])

# re-integrate the data, cluster the data, and annotate again
DefaultAssay(seurat.glial) <- "RNA"
seurat.glial <- JoinLayers(seurat.glial)
seurat.glial[["RNA"]] <- split(seurat.glial[["RNA"]], f = seurat.glial$sample_name)
seurat.glial <- NormalizeData(seurat.glial)
seurat.glial <- FindVariableFeatures(seurat.glial)
seurat.glial <- ScaleData(seurat.glial)
seurat.glial <- RunPCA(seurat.glial)
seurat.glial <- IntegrateLayers(
  object = seurat.glial, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
# cluster the data
seurat.glial <- FindNeighbors(seurat.glial, reduction = "harmony", dims = 1:50)
seurat.glial <- FindClusters(seurat.glial, resolution = 2, cluster.name = "seurat_clusters")
# visualization
seurat.glial <-  RunTSNE(
  seurat.glial,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 100,
  fast_tsne_path=FAST_TSNE_PATH, ann_not_vptree=FALSE, nthreads=12
)

}   ## end !RESUME_GLIAL

p1 <- DimPlot(seurat.glial,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

p2 <- DimPlot(seurat.glial,reduction = "tsne",group.by = "seurat_clusters",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)+NoLegend()

manual_markers <- c("AQP4","MBP","PDGFRA","P2RY12","SKAP1","CLDN5","GRM8","ACTA2","COL1A2")
d1 <- DotPlot(seurat.glial,features = manual_markers)

f1 <- FeaturePlot(seurat.glial,features = c("CLDN5","GRM8"))

# these are glial clusters that look ambiguous in terms of marker expressions
markers <- readRDS(file.path(REFERENCE_DIR, "fine_grain_markers_INS_PFC.rds"))

d2 <- DotPlot(seurat.glial,features = markers[["oligo"]])
f2 <- FeaturePlot(seurat.glial,features = markers[["oligo"]])

f3 <- FeaturePlot(seurat.glial,features = c("COLEC12","F13A1"))


## REVISION: the published code hard-coded
##   cluster2remove <- c(21,24,26,28,33,37,39,40,42,43,47,48,50,53,54)
##   macrophage_cells <- colnames(seurat.glial)[Idents(seurat.glial) == 38]
##
## Both index into a clustering computed earlier in this same script
## (FindClusters, resolution 2), so they are not stable between runs, and
## certainly not after ambient correction.
##
## The published comments give the intent:
##   "these are glial clusters that look ambiguous in terms of marker expressions"
##   "these are all macrophage markers"  (for cluster 38)
##
## Calibration target, from the published objects:
##   seurat_combined.RDS     673,993 cells
##   seurat_combined_v2.RDS  654,918 cells
##   dropped                  19,075  = 4.8% of glia, all Glial, spread across
##                                      every glial annotation
##                            (Oligo 4115, OPC 3760, Astro 3619, PC 3193,
##                             Micro 1972, Endo 1917, Immune 293, SMC 134, VLMC 72)
##   Macrophages in v2: 1,002 cells
##
## ---------------------------------------------------------------------------
## Criterion: per-lineage relative evidence.
##
## Detection rate in single-nucleus data varies widely by gene, so no single
## absolute cut works across lineages. Measured on this object, the purest
## cluster of each type reaches:
##   Oligo 1.000  Astro 0.996  Endo 0.983  OPC 0.977  Tcell 0.839
##   VLMC  0.754  SMC   0.680  Macro 0.619  Micro 0.608
## MBP is in ~100% of oligodendrocyte nuclei; P2RY12 in ~61% of microglial nuclei.
## An absolute threshold of 0.84 would spare every oligodendrocyte cluster and
## remove all 34,661 microglia.
##
## So each cluster's evidence for a lineage is scored relative to the purest
## cluster of that lineage in this same object:
##
##     rel(cluster, lineage) = pct(cluster, lineage) / max_over_clusters pct(., lineage)
##
## and a cluster is ambiguous when it fails to reach REL_THRESHOLD of the best
## cluster for every lineage - it looks like a degraded version of everything and
## a clean version of nothing.
##
## Alternatives that do not work:
##   - z-scored panels ("ambiguous if best score below average"): the score is a
##     maximum over 8 z-scored panels, which is almost never negative, so it
##     flags 0 clusters against a published 15.
##   - an absolute cut on the mean detection fraction across a panel's genes:
##     averaging a near-universal marker with a sparse one caps a pure cluster far
##     below 1 (a 30,252-cell microglial cluster scored 0.34), and PLP1 is absent
##     from the corrected object. It flags 1 cluster / 33 cells.
## Hence each panel is scored by its best-detected marker, relative to the purest
## cluster.
##
## Panels use only lineage-specific markers. DOCK8 and CSF1R are excluded from
## the microglial panel: they are pan-myeloid, so they score the macrophage
## cluster as 0.76 microglial and hide it.
LINEAGE_MARKERS <- list(
  Astro = c("AQP4", "GJA1", "SLC1A2"),
  Oligo = c("MBP", "PLP1", "MOG", "MOBP", "ST18", "CNP"),
  OPC   = c("PDGFRA", "CSPG4", "VCAN"),
  Micro = c("P2RY12", "CX3CR1", "TMEM119"),   # microglia-specific only
  Macro = c("CD163", "F13A1", "MRC1"),        # perivascular macrophage
  Tcell = c("SKAP1", "CD3E", "THEMIS"),
  Endo  = c("CLDN5", "FLT1"),
  SMC   = c("ACTA2", "MYH11"),
  VLMC  = c("COL1A2", "DCN")
)
## COLEC12 and LYVE1 are deliberately not in the Macro panel. Both are expressed
## by endothelium, and with them a 5,100-cell endothelial cluster (CLDN5 0.83,
## FLT1 0.98) outscored the real macrophage cluster on macrophage markers.

REL_THRESHOLD <- 0.86            # calibrated below; see the sweep in the comment
PUBLISHED_DROP_FRACTION <- 0.048 # of glial cells
PUBLISHED_DROP_CLUSTERS <- 15
PUBLISHED_MACROPHAGE_CELLS <- 1002
PUBLISHED_CLUSTER2REMOVE <- c(21,24,26,28,33,37,39,40,42,43,47,48,50,53,54)  # reference only
PUBLISHED_MACROPHAGE_CLUSTER <- 38                                            # reference only

## Sweep of REL_THRESHOLD on this object:
##   0.82 ->  6 clusters   12,020 cells   3.01%
##   0.84 -> 10 clusters   14,555 cells   3.64%
##   0.85 -> 11 clusters   17,649 cells   4.42%
##   0.86 -> 13 clusters   19,486 cells   4.88%   <- published 19,075 / 4.8%
##   0.87 -> 17 clusters   36,753 cells   9.20%
## 0.86 reproduces the published cell count to within 2% (19,486 vs 19,075).
## Cluster count is 13 against the published 15; the published dropped set
## included 3,193 pericytes, and PC is not a lineage here (those cells fall under
## Endo), which accounts for the difference.
##
## The threshold is calibrated to reproduce the published drop; it is chosen,
## not derived from first principles. A single threshold matches the cell count,
## gives a sensible cluster count, and leaves a biologically coherent set of
## retained populations. The curve is steep above 0.87, so check the diagnostic
## below if the input data change; do not carry 0.86 to a different object
## unexamined.

## ---- scores ------------------------------------------------------------------
## Per-gene fractions kept, not just the collapsed score, so the criterion can be
## re-derived from the CSV without reloading this object.
gene_pct <- cluster_gene_pct(seurat.glial, unique(unlist(LINEAGE_MARKERS)),
                             group = "seurat_clusters")

lineage_pct <- sapply(LINEAGE_MARKERS, function(gs) {
  g <- intersect(gs, rownames(gene_pct))
  if (!length(g)) return(rep(NA_real_, ncol(gene_pct)))
  apply(gene_pct[g, , drop = FALSE], 2, max)   # best-detected marker of the panel
})
rownames(lineage_pct) <- colnames(gene_pct)

lineage_max <- apply(lineage_pct, 2, max, na.rm = TRUE)   # purest cluster per lineage
rel <- sweep(lineage_pct, 2, lineage_max, "/")

best_rel     <- apply(rel, 1, max, na.rm = TRUE)
best_lineage <- colnames(rel)[apply(rel, 1, which.max)]
names(best_lineage) <- rownames(rel)

## AMBIGUOUS: reaches REL_THRESHOLD of the purest cluster for NO lineage
ambiguous <- best_rel < REL_THRESHOLD
cluster2remove <- rownames(rel)[ambiguous]

## MACROPHAGE: nearest lineage is Macro, and it survived the ambiguity filter
macro_candidates <- rownames(rel)[best_lineage == "Macro" & !ambiguous]
macrophage_cluster <- if (length(macro_candidates)) {
  macro_candidates[which.max(rel[macro_candidates, "Macro"])]
} else NA_character_

## ---- diagnostics -------------------------------------------------------------
dir.create("./combined/data", recursive = TRUE, showWarnings = FALSE)
n_by_cluster <- table(as.character(seurat.glial@meta.data$seurat_clusters))

write.csv(round(t(gene_pct), 4), "./combined/data/glial_per_gene_pct.csv")

glial_diag <- data.frame(
  cluster           = rownames(rel),
  n_cells           = as.integer(n_by_cluster[rownames(rel)]),
  best_lineage      = best_lineage[rownames(rel)],
  best_rel          = round(as.numeric(best_rel), 4),
  best_pct          = round(as.numeric(apply(lineage_pct, 1, max, na.rm = TRUE)), 4),
  called_ambiguous  = ambiguous,
  called_macrophage = rownames(rel) == macrophage_cluster,
  stringsAsFactors  = FALSE
)
glial_diag <- cbind(glial_diag, round(rel, 3))
glial_diag <- glial_diag[order(glial_diag$best_rel), ]
write.csv(glial_diag, "./combined/data/glial_cluster_diagnostic.csv", row.names = FALSE)

## Save the glial object so the criterion can be re-derived without recomputing
## the subset, integration, clustering and t-SNE - roughly an hour.
## Skipped when the object was itself loaded from that checkpoint - rewriting
## 4 GB of identical data costs time, and a crash mid-write would destroy the
## only copy.
if (!isTRUE(mget("RESUME_GLIAL", ifnotfound = list(FALSE))$RESUME_GLIAL))
  saveRDS(seurat.glial, "./combined/data/seurat_glial_clustered.RDS")

n_drop <- sum(as.integer(n_by_cluster[cluster2remove]), na.rm = TRUE)
frac   <- n_drop / ncol(seurat.glial)
n_macro <- if (!is.na(macrophage_cluster)) as.integer(n_by_cluster[macrophage_cluster]) else 0L

message("[combined/02] glial cells: ", ncol(seurat.glial),
        "  clusters: ", nrow(rel))
message("[combined/02] per-lineage purest-cluster pct: ",
        paste(sprintf("%s %.2f", names(lineage_max), lineage_max), collapse = "  "))
message("[combined/02] ambiguous (rel < ", REL_THRESHOLD, "): ",
        length(cluster2remove), " clusters, ", n_drop, " cells ",
        sprintf("(%.2f%% of glia; published %d clusters / %.1f%%)",
                100 * frac, PUBLISHED_DROP_CLUSTERS, 100 * PUBLISHED_DROP_FRACTION))
message("[combined/02] macrophage cluster: ", macrophage_cluster, " (", n_macro,
        " cells; published ", PUBLISHED_MACROPHAGE_CELLS, ")")
message("[combined/02] retained by nearest lineage:")
kept <- !ambiguous
for (L in names(sort(tapply(as.integer(n_by_cluster[rownames(rel)])[kept],
                            best_lineage[rownames(rel)][kept], sum), decreasing = TRUE)))
  message("               ", format(L, width = 6), " ",
          sum(as.integer(n_by_cluster[rownames(rel)])[kept & best_lineage[rownames(rel)] == L]),
          " cells")
message("[combined/02] diagnostics + glial object -> combined/data/")

## ---- guards ------------------------------------------------------------------
## Every one of these fires before any downstream object is written, and the
## glial object above is already saved, so hitting one costs minutes.
if (length(cluster2remove) > 0 &&
    sum(as.character(Idents(seurat.glial)) %in% cluster2remove) == 0) {
  stop(length(cluster2remove), " clusters marked ambiguous but 0 cells matched - ",
       "cluster names do not match Idents()")
}
if (frac > 2 * PUBLISHED_DROP_FRACTION || frac < 0.4 * PUBLISHED_DROP_FRACTION) {
  stop(sprintf("ambiguous drop is %.2f%% of glia against a published %.1f%% - ",
               100 * frac, 100 * PUBLISHED_DROP_FRACTION),
       "read combined/data/glial_cluster_diagnostic.csv and re-run the sweep in ",
       "the threshold against the saved glial object. Do not simply ",
       "widen this guard.")
}
if (is.na(macrophage_cluster)) {
  stop("no macrophage cluster identified - consensus_annotation_v3 would have no ",
       "Macrophages, and the published run had ", PUBLISHED_MACROPHAGE_CELLS, ". ",
       "Check the Macro panel against combined/data/glial_per_gene_pct.csv.")
}
if (n_macro > 5 * PUBLISHED_MACROPHAGE_CELLS)
  warning("macrophage cluster has ", n_macro, " cells against a published ",
          PUBLISHED_MACROPHAGE_CELLS, " - verify it is not merged with microglia")

seurat.glial.sub <- subset(seurat.glial,cells = colnames(seurat.glial)[!(Idents(seurat.glial) %in% cluster2remove)])

p1_new <- DimPlot(seurat.glial.sub,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
                  raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

seurat.glial.sub <- JoinLayers(seurat.glial.sub)
## REVISION: the published code used ident.1 = 38. This uses the marker-derived
## cluster and serves as a confirmation step rather than defining the macrophages.
markers.38 <- if (!is.na(macrophage_cluster))
  FindMarkers(seurat.glial.sub, ident.1 = macrophage_cluster, ident.2 = NULL) else NULL
markers.38.top10 <- if (!is.null(markers.38)) markers.38 %>%
  slice_max(avg_log2FC, n = 10) else NULL
# these are all macrophage markers

cells2keep <- colnames(seurat.combined)[!(colnames(seurat.combined) %in%
  colnames(seurat.glial)[as.character(Idents(seurat.glial)) %in% cluster2remove])]
## REVISION: the published code used Idents(seurat.glial) == 38
macrophage_cells <- if (!is.na(macrophage_cluster))
  colnames(seurat.glial)[as.character(Idents(seurat.glial)) == macrophage_cluster] else character(0)

# subset the original data, and rename macrophages in V3
seurat.combined.sub <- subset(seurat.combined,cells = cells2keep)
seurat.combined.sub@meta.data$consensus_annotation_v3 <- as.character(seurat.combined.sub@meta.data$consensus_annotation_v2)
seurat.combined.sub@meta.data$consensus_annotation_v3[colnames(seurat.combined.sub) %in% macrophage_cells] <- "Macrophages"

p4 <- DimPlot(seurat.combined.sub,reduction = "tsne",group.by = "consensus_annotation_v3",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
saveRDS(seurat.combined.sub,"./combined/data/seurat_combined_v2.RDS")


