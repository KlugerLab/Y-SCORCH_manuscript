# ==============================================================================
# Figure 2E: RNA vs gene-activity log2 FC scatter, split by cell class.
#
# Same genes and values as fig2E_rna_vs_geneactivity_scatter_combined
# (revision/tools/fig2_atac_summary.R), which puts ~100 points in 15 similar
# cell-type colours on one plot. Here there is one panel per cell class
# (excitatory neurons / inhibitory neurons & MSNs / glia & vascular), cell types
# in the published palette (cp_col), region as shape, free axes per panel.
# Every point is labelled and no gene is forced.
#
# INPUT   OUT_ROOT/04-ATAC_DE/figures/fig2/fig2E_rna_vs_geneactivity_scatter.csv
# OUTPUT  same folder, fig2E_rna_vs_geneactivity_by_class.png/.pdf/_caption.txt (2 rows)
#         and fig2E_rna_vs_geneactivity_by_class_row.* (one row)
#
# Usage: revision/docker/run.sh revision/tools/fig2E_scatter_by_class.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

DIR <- file.path(OUT_ROOT, "04-ATAC_DE", "figures", "fig2")
sc  <- utils::read.csv(file.path(DIR, "fig2E_rna_vs_geneactivity_scatter.csv"), stringsAsFactors = FALSE)
REGIONS <- c("PFC", "INS", "VST")
REGION_SHAPES <- c(PFC = 16, INS = 17, VST = 15)
EXC <- c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L6b", "L5 ET", "L5/6 NP", "L6 CT")
GLIA <- c("Astro", "Oligo", "OPC", "Micro", "Macrophages", "Other immune", "PC", "Endo", "SMC", "VLMC")
CLASSES <- c("Excitatory neurons", "Inhibitory neurons & MSNs", "Glia & vascular")
sc$class <- ifelse(sc$cell_type %in% EXC, CLASSES[1], ifelse(sc$cell_type %in% GLIA, CLASSES[3], CLASSES[2]))
## Published cell-type palette (cp_col, 02-ATAC_processing/plotting/01-TSNE_markers.R).
## Points are filled shapes with a thin dark outline so the pale entries (khaki2 OPC,
## skyblue1 Sst) stay visible.
cp_col <- setNames(c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1",
  "purple1", "lightslateblue","orchid3","pink3","plum1",  "tomato1", "tomato3","sienna2", "burlywood3",
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"),
  c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
    "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
    "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
    "Astro","Oligo","OPC","Micro","Macrophages","Other immune","PC","Endo","SMC","VLMC"))
stopifnot(all(sc$cell_type %in% names(cp_col)))
FILLED_SHAPES <- c(PFC = 21, INS = 24, VST = 22)

n_reg <- table(factor(sc$region, levels = REGIONS))
## Layout: neurons side by side on top, glia & vascular (the most points) full
## width below. Each panel carries its own cell-type key on the right; the region
## (shape) key appears once, with the excitatory neurons.
panel <- function(cl, region_key = FALSE, row = FALSE, lab_size = 2.1) {
  d <- sc[sc$class == cl, ]
  lv <- intersect(names(cp_col), unique(d$cell_type))                # published order
  d$cell_type <- factor(d$cell_type, levels = lv)
  d$region <- factor(d$region, levels = REGIONS)
  ggplot(d, aes(x = logFC_rna, y = logFC_ga, fill = cell_type, shape = region)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey40") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey40") +
    geom_point(size = 1.9, stroke = 0.25, colour = "grey15") +
    ## a key is only drawn for a region present in the panel, and the excitatory
    ## panel has no VST points - one invisible point per region fixes the key
    geom_point(data = transform(d[rep(1, length(REGIONS)), ], region = factor(REGIONS, levels = REGIONS)),
               fill = "grey50", colour = "grey15", alpha = 0, size = 1.9, stroke = 0.25,
               show.legend = region_key, inherit.aes = TRUE) +
    ggrepel::geom_text_repel(aes(label = gene), size = lab_size, fontface = "italic", colour = "black",
                             max.overlaps = Inf, box.padding = 0.22, point.padding = 0.12, force = 4,
                             force_pull = 0.5, max.time = 5, max.iter = 1e6,
                             min.segment.length = 0.15, segment.size = 0.2, segment.colour = "grey55",
                             seed = 1, show.legend = FALSE) +
    scale_fill_manual(values = cp_col[lv], name = NULL) +
    scale_shape_manual(values = FILLED_SHAPES, breaks = REGIONS, drop = FALSE,
                       labels = REGIONS, name = "Region") +          # counts are in the caption
    labs(title = cl,                                   # no counts in titles
         x = if (!row || cl == CLASSES[2]) "log2 FC, gene expression (RNA)" else NULL,
         y = if (!row || cl == CLASSES[1]) "log2 FC, gene activity (ATAC)" else NULL) +
    theme_astro_pub(base_size = 7) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
          legend.position = if (row) "bottom" else "right", legend.box = "horizontal",
          legend.justification = if (row) "left" else "center", legend.text = element_text(size = 6.5),
          legend.title = element_text(size = 7, face = "bold"), legend.key.size = unit(0.28, "cm")) +
    guides(shape = if (region_key) guide_legend(order = 1, ncol = 1,
                                               override.aes = list(size = 2, shape = unname(FILLED_SHAPES), fill = "grey50", alpha = 1)) else "none",
           fill = guide_legend(order = 2, ncol = if (row) 2 else 1, override.aes = list(size = 2.4, shape = 21)))
}
p <- wrap_plots(A = panel(CLASSES[1], region_key = TRUE), B = panel(CLASSES[2]), C = panel(CLASSES[3]),
                design = "AB\nCC", heights = c(1, 1.15))
f <- file.path(DIR, "fig2E_rna_vs_geneactivity_by_class.png")
save_pub(p, f, width = 7.2, height = 6.4)
cap <- paste0(
  "Effect sizes of genes changed in PWH in both gene expression (x, RNA pseudobulk edgeR) and gene activity (y, ",
  "ATAC reads over gene body and promoter, pseudobulk edgeR) in the same cell type and region (each FDR < 0.05, ",
  "|log2FC| > 0.25, detected in > 5% of nuclei), one panel per cell class; colour, cell type; shape, region. ", round(100 * mean(sc$concordant)), "% of points change in the ",
  "same direction in both modalities. Every gene is labelled. Genes per region: ",
  paste(sprintf("%s %d", REGIONS, as.integer(n_reg)), collapse = ", "), ".")
astro_caption(f, cap)

## ONE-ROW version: the three classes side by side at full page width, glia
## wider (most points); each cell-type key under its panel, region key with glia;
## y title on the first panel, x title on the middle one; labels slightly smaller.
pr <- wrap_plots(panel(CLASSES[1], region_key = TRUE, row = TRUE, lab_size = 1.9),
                 panel(CLASSES[2], row = TRUE, lab_size = 1.9),
                 panel(CLASSES[3], row = TRUE, lab_size = 1.9),
                 nrow = 1, widths = c(1, 1, 1.8))   # glia wider
fr <- file.path(DIR, "fig2E_rna_vs_geneactivity_by_class_row.png")
## neuron panels 2.15 in each (7.2 / 3.35), glia panel 1.8x that width
save_pub(pr, fr, width = round(7.2 / 3.35 * 3.8, 1), height = 3.6)
astro_caption(fr, cap)
message("[fig2E-class] ", paste(sprintf("%s %d", CLASSES, as.integer(table(factor(sc$class, levels = CLASSES)))), collapse = "; "))
message("[fig2E-class] done")
