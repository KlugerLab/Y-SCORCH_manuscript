# ==============================================================================
# Module-defined microglial subgroups (supplemental module_subgroup panel) on UMAP.
#
# Method: the subgroup labels come from Ruiqi's module analysis and are stored
#   as `module_subgroup` on the PWH microglia object: AUCell AUC per module,
#   thresholded by AUCell_exploreThresholds, composite label M1 / M2 / M3 /
#   M1+M2 ... / Low, with "M2+M3" folded into "M3". Nothing is recomputed here;
#   only the embedding changes.
#   Embedding: the cached UMAP from revision/tools/micro_module_umap.R (RunUMAP
#   on the object's own harmony reduction, dims 1:30, seed 3). If the cache is
#   missing the script stops rather than making a second, differently-seeded
#   embedding.
#   Style as the published tSNE panel: M1 #2166AC, M2 #B2182B, M3 #D4A017, Low
#   grey80; plotting order Low, M1, M2, M3 so labelled cells sit on top; no
#   axes; 7 x 6 in. The same plot is drawn on the object's tSNE as a check.
#
# Inputs: OUT_ROOT/microglia_module_analysis/umap/umap_coords_PWH.csv
#         OUT_ROOT/microglia_module_analysis/input/scorch_microglia_HIV.rds
# Output: OUT_ROOT/microglia_module_analysis/umap/
#           fig5_subgroup_umap.png/.pdf, fig5_subgroup_tsne_check.png,
#           fig5_subgroup_counts.csv, module_subgroup_labels_PWH.csv (label cache)
#
# Usage: revision/docker/run.sh revision/tools/micro_subgroup_umap.R
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
BASE <- file.path(OUT_ROOT, "microglia_module_analysis")
OUT  <- file.path(BASE, "umap")
COLS <- c(M1 = "#2166AC", M2 = "#B2182B", M3 = "#D4A017", Low = "grey80")

f_umap <- file.path(OUT, "umap_coords_PWH.csv")
if (!file.exists(f_umap)) stop("no cached UMAP: run revision/tools/micro_module_umap.R first")
u <- utils::read.csv(f_umap, stringsAsFactors = FALSE)

## the labels are cached on the first run, so a restyle does not reload the object
f_lab <- file.path(OUT, "module_subgroup_labels_PWH.csv")
o <- NULL
if (file.exists(f_lab)) {
  md <- utils::read.csv(f_lab, stringsAsFactors = FALSE); rownames(md) <- md$cell
  message("[subumap] labels read from ", basename(f_lab))
} else {
  o <- readRDS(file.path(BASE, "input", "scorch_microglia_HIV.rds"))
  md <- o@meta.data
}
if (!"module_subgroup" %in% colnames(md))
  stop("no module_subgroup column on the PWH object; it is written by the microglia module pipeline (Fig 5d-e)")
md$cell <- rownames(md)
if (!file.exists(f_lab)) utils::write.csv(md[, c("cell", "module_subgroup")], f_lab, row.names = FALSE)
d <- dplyr::inner_join(u, md[, c("cell", "module_subgroup")], by = "cell")
message(sprintf("[subumap] %d of %d PWH nuclei matched the cached UMAP", nrow(d), nrow(md)))
lev <- c("Low", setdiff(sort(unique(d$module_subgroup)), "Low"))
d$module_subgroup <- factor(d$module_subgroup, levels = lev)
print(table(d$module_subgroup))
extra <- setdiff(levels(d$module_subgroup), names(COLS))
if (length(extra)) {                       # composite labels the published panel did not keep
  message("[subumap] extra labels kept with hue colours: ", paste(extra, collapse = ", "))
  COLS <- c(COLS, setNames(scales::hue_pal()(length(extra)), extra))
}
cnt <- as.data.frame(table(subgroup = d$module_subgroup))
cnt$pct <- round(100 * cnt$Freq / sum(cnt$Freq), 1)
utils::write.csv(cnt, file.path(OUT, "fig5_subgroup_counts.csv"), row.names = FALSE)

draw <- function(df, x, y, title) {
  df <- df[order(df$module_subgroup), ]          # Low first, labelled subgroups on top
  ggplot(df, aes(.data[[x]], .data[[y]], colour = module_subgroup)) +
    geom_point(size = 0.5, stroke = 0) +
    scale_colour_manual(values = COLS, name = NULL, drop = FALSE) +
    ggtitle(title) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.text = element_text(size = 10)) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
}
p <- draw(d, "umap_1", "umap_2", "module_subgroup")   # the published panel title
ggsave(file.path(OUT, "fig5_subgroup_umap.png"), p, width = 7, height = 6, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig5_subgroup_umap.pdf"), p, width = 7, height = 6, device = grDevices::cairo_pdf, bg = "white")

if (!is.null(o) && "tsne" %in% names(o@reductions)) {
  ts <- as.data.frame(Embeddings(o, "tsne")); names(ts) <- c("tsne_1", "tsne_2"); ts$cell <- rownames(ts)
  t2 <- dplyr::inner_join(ts, md[, c("cell", "module_subgroup")], by = "cell")
  t2$module_subgroup <- factor(t2$module_subgroup, levels = lev)
  ggsave(file.path(OUT, "fig5_subgroup_tsne_check.png"),
         draw(t2, "tsne_1", "tsne_2", "module_subgroup (tSNE check)"),
         width = 7, height = 6, dpi = 300, bg = "white")
}
message("[subumap] done -> ", OUT)
