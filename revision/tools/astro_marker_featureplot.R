# ==============================================================================
# Feature plots of contamination markers (microglial APBB1IP, endothelial CLDN5 and
# FLT1 by default) on the astrocyte embeddings, with % of nuclei expressing per subtype.
#
# Method
#   The Endo-Astro cluster (Rmd cluster 16) was called with the panel CLDN5, FLT1,
#   PECAM1, VWF, ABCB1, PTPRB (astro_vascular_check in
#   revision/tools/astro_cluster_features.R); CLDN5 and FLT1 are the two primary
#   ones. Micro-Astro (clusters 8 and 18) was called with P2RY12, CX3CR1, TMEM119,
#   CSF1R, C1QB, DOCK8, AIF1; APBB1IP belongs to the same microglial set.
#   All clusters are plotted, so the Oligo and Endo-Astro clusters that downstream
#   analyses drop remain visible. Expression is the CellBender-corrected assay, as
#   for every other marker visualisation. The UMAP is the cached umap_astro
#   embedding (harmony_astro dims 1:30, as Fig1D); the object's stock "umap" is a
#   different reduction.
#
# Inputs   OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled_no1040vst.RDS
#          OUT_ROOT/06-astrocyte_analysis/data/astro_umap_embeddings_no1040vst.rds
# Outputs  OUT_ROOT/06-astrocyte_analysis/figures/
#            markers_contamination_<reduction>[_withlabels].png/.pdf  one panel per gene
#            markers_contamination_pct.csv      % of nuclei expressing, per subtype
#
# Usage: revision/docker/run.sh revision/tools/astro_marker_featureplot.R
#          [--genes APBB1IP,CLDN5,FLT1] [--reductions tsne,umap_astro] [--with-labels]
# ==============================================================================

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork); library(dplyr) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_cluster_features.R"))
## ASTRO_SUBTYPE_COLOURS / ASTRO_EXTRA_COLOURS live here
source(file.path(SCORCH_REPO, "revision/tools/astro_reviewer_figures.R"))

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
GENES <- trimws(strsplit(getarg("--genes", "APBB1IP,CLDN5,FLT1"), ",")[[1]])
REDS  <- trimws(strsplit(getarg("--reductions", "tsne,umap_astro"), ",")[[1]])
WITH_LABELS <- "--with-labels" %in% args
LAB   <- "subcluster_labels"
BASE  <- file.path(OUT_ROOT, "06-astrocyte_analysis")
OUT   <- file.path(BASE, "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(BASE, "data", "astro_all_clusters_labeled_no1040vst.RDS"))
message("[mk] ", ncol(obj), " nuclei; reductions: ", paste(names(obj@reductions), collapse = ", "))
## The object carries a stock "umap" inherited from upstream, not the astrocyte UMAP.
## The astrocyte one is the cached umap_astro embedding (harmony_astro, dims 1:30)
## that Fig1D uses; always take that.
f <- file.path(BASE, "data", "astro_umap_embeddings_no1040vst.rds")
if (file.exists(f)) {
  em <- as.matrix(readRDS(f))[colnames(obj), , drop = FALSE]
  colnames(em) <- paste0("umapastro_", seq_len(ncol(em)))
  obj[["umap_astro"]] <- CreateDimReducObject(embeddings = em, key = "umapastro_", assay = DefaultAssay(obj))
  message("[mk] umap_astro added from ", basename(f))
} else stop("astro UMAP cache missing: ", f)

miss <- setdiff(GENES, rownames(obj))
if (length(miss)) message("[mk] not in the object, skipped: ", paste(miss, collapse = ", "))
GENES <- intersect(GENES, rownames(obj))

## % of nuclei expressing, per subtype - the numbers behind the panels
pct <- cluster_gene_pct(obj, GENES, group = LAB)
utils::write.csv(round(100 * as.data.frame(pct), 2), file.path(OUT, "markers_contamination_pct.csv"))
print(round(100 * t(as.data.frame(pct)), 1))

cols <- c(ASTRO_SUBTYPE_COLOURS, ASTRO_EXTRA_COLOURS)
Idents(obj) <- LAB
for (rd in intersect(REDS, names(obj@reductions))) {
  ## The labelled subtype panel is drawn only with --with-labels.
  p0 <- if (WITH_LABELS) DimPlot(obj, reduction = rd, group.by = LAB, label = TRUE, label.size = 3, repel = TRUE,
                cols = cols[intersect(names(cols), unique(as.character(obj@meta.data[[LAB]])))], raster = FALSE) +
    ggtitle("astrocyte subtypes") + NoAxes() + NoLegend() else NULL
  ps <- lapply(GENES, function(g)
    FeaturePlot(obj, features = g, reduction = rd, order = TRUE, raster = FALSE, pt.size = 0.2) +
      scale_colour_gradientn(colours = c("grey88", "#FDBF6F", "#B2182B")) +
      ggtitle(g) + NoAxes() +
      theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 13)))
  panels <- if (is.null(p0)) ps else c(list(p0), ps)
  ncol <- if (is.null(p0)) length(panels) else 2
  fig <- wrap_plots(panels, ncol = ncol)
  w <- if (is.null(p0)) 3.1 * ncol else 9
  h <- if (is.null(p0)) 3.4 else 8
  f <- file.path(OUT, sprintf("markers_contamination_%s%s.png", rd, if (WITH_LABELS) "_withlabels" else ""))
  ggsave(f, fig, width = w, height = h, dpi = 300, bg = "white")
  ggsave(sub("\\.png$", ".pdf", f), fig, width = w, height = h, device = grDevices::cairo_pdf, bg = "white")
  message("[mk] wrote ", f)
}
message("[mk] done")
