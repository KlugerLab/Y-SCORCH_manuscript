# ==============================================================================
# Figure 4: region-pair logFC scatters, rebuilt for readability (reviewer request).
#
# The published panels (03-DE_analysis/06-compare_logFC_cross_region.R,
# corr_<A>_<B>_logfc_labeled.png) are hard to read for four reasons, all fixed here:
#   1 every gene passing the label rule was passed to geom_label_repel - 1,095 of
#     them in INS-PFC - and max.overlaps silently dropped all but ~50, so WHICH
#     genes got labelled was arbitrary. Here a fixed rule picks N_LABEL genes.
#   2 the same gene was labelled once per cell type (AC105402.3 six times). Here
#     labels are one per gene, at its most extreme cell type.
#   3 filled label boxes in 36 cell-type colours dominated the panel. Here: plain
#     italic text with a thin leader line.
#   4 the figure was printed with NoLegend(), so the colours could not be read.
#     Here there is a legend, and the per-panel correlation is printed in the corner.
# Non-significant genes are grey, so the coloured points are the ones that matter.
#
# Two colour variants are written:
#   celltype  DE-in-both points in the published cell-type palette (cp_col)
#   class     the same points in the published coarse palette (cp_col2:
#             Ext / MSNs / Inh / Glial) - four colours instead of sixteen
# A published-style version (same look as the published panels, fewer labels)
# and a separate cell-type key are also written.
#
# INPUT   OUT_ROOT/03-DE_analysis/<R>/results/<R>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#         (regional DE with donor 1040's VST cells excluded)
# OUTPUT  OUT_ROOT/03-DE_analysis/figures/fig4/
#           fig4_region_pair_scatter_{celltype,class}.png/.pdf/_caption.txt
#           fig4_region_pair_scatter_labels.csv   every labelled gene
#           fig4_region_pair_stats.csv            n and r per pair
#           fig4_<A>_<B>_published_style[sfx].png/.pdf
#           fig4_published_style_legend[_narrow|_column].png/.pdf
#
# Usage: revision/docker/run.sh revision/tools/fig4_cross_region_scatter.R [--n-label 20]
#          [--n-pub 20[,40]] [--label-size 5] [--force Oligo,OPC] [--per-type 1]
#          [--genes CD163:Micro] [--updown]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

args    <- commandArgs(trailingOnly = TRUE)
getarg  <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
N_LABEL <- as.integer(getarg("--n-label", "20"))
OUT <- file.path(OUT_ROOT, "03-DE_analysis", "figures", "fig4")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(SCORCH_REPO, "revision/tools/fig4_common.R"))   # palette, get(), prep(), pub_lab()

REG <- lapply(setNames(c("INS", "PFC", "VST"), c("INS", "PFC", "VST")), get)

D <- lapply(PAIRS, prep)
names(D) <- vapply(PAIRS, function(p) paste(p, collapse = "-"), character(1))

## labels: the N_LABEL strongest genes that are DE in BOTH regions, one label per
## gene (its most extreme cell type), so no gene is repeated - no forced labels
lab_of <- function(d) {
  d %>% dplyr::filter(state == "both") %>% arrange(desc(size_metric)) %>%
    distinct(gene, .keep_all = TRUE) %>% head(N_LABEL)
}
LAB <- lapply(D, lab_of)
utils::write.csv(bind_rows(LAB) %>%
                   dplyr::select(pair, gene, cell_type, logFC.x, logFC.y, p_adj.glb.x, p_adj.glb.y, concordant),
                 file.path(OUT, "fig4_region_pair_scatter_labels.csv"), row.names = FALSE)
st <- bind_rows(lapply(names(D), function(k) {
  d <- D[[k]]
  data.frame(pair = k, n = nrow(d), n_both = sum(d$state == "both"),
             r = cor(d$logFC.x, d$logFC.y),
             pct_concordant_both = round(100 * mean(d$concordant[d$state == "both"])))
}))
utils::write.csv(st, file.path(OUT, "fig4_region_pair_stats.csv"), row.names = FALSE)
print(st, row.names = FALSE, digits = 2)

panel <- function(k, by = c("celltype", "class"), legend = FALSE) {
  by <- match.arg(by)
  d  <- D[[k]]; l <- LAB[[k]]
  bg <- d[d$state != "both", ]; fg <- d[d$state == "both", ]
  pal <- if (by == "celltype") cp_col[intersect(LEV, unique(unlist(lapply(D, function(z) as.character(z$cell_type[z$state == "both"])))))] else cp_col2
  fg$col <- if (by == "celltype") fg$cell_type else fg$coarse
  l$col  <- if (by == "celltype") l$cell_type  else l$coarse
  rr <- st$r[st$pair == k]
  ggplot(mapping = aes(x = logFC.x, y = logFC.y)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_point(data = bg, colour = "grey85", size = 0.5, stroke = 0, alpha = 0.6) +
    geom_point(data = fg, aes(colour = col), size = 1.1, stroke = 0, alpha = 0.95) +
    ggrepel::geom_text_repel(data = l, aes(label = gene), size = 2.1, fontface = "italic",
                             colour = "grey10", max.overlaps = Inf, box.padding = 0.3,
                             point.padding = 0.15, force = 4, max.time = 5, max.iter = 1e6,
                             min.segment.length = 0.15, segment.size = 0.2, segment.colour = "grey55",
                             seed = 1) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.6, size = 2.4,
             label = sprintf("r = %.2f", rr)) +
    ## limits, not just values: a key is otherwise drawn only for the levels present
    ## in that panel, the three legends then differ and patchwork cannot collect them
    scale_colour_manual(values = pal, limits = names(pal), name = NULL, drop = FALSE) +
    labs(title = sub("-", " vs ", k),
         x = sprintf("log2 FC, %s", D[[k]]$x_reg[1]), y = sprintf("log2 FC, %s", D[[k]]$y_reg[1])) +
    theme_astro_pub(base_size = 7) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
          legend.position = if (legend) "bottom" else "none",
          legend.text = element_text(size = 6.5), legend.key.size = unit(0.25, "cm")) +
    guides(colour = guide_legend(nrow = if (by == "celltype") 4 else 1,
                                 override.aes = list(size = 2, alpha = 1)))
}

for (by in c("celltype", "class")) {
  ## every panel carries the same scale, so patchwork collects one shared legend
  p <- wrap_plots(lapply(names(D), panel, by = by, legend = TRUE), nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.justification = "center")
  f <- file.path(OUT, paste0("fig4_region_pair_scatter_", by, ".png"))
  save_pub(p, f, width = 7.2, height = if (by == "celltype") 3.9 else 3.4)
  astro_caption(f, paste0(
    "Fold changes in PWH compared between regions, for every gene tested in the same cell type in both regions ",
    "(detected in > 5% of nuclei). Coloured points are significant in BOTH regions (FDR < 0.05, |log2FC| > 0.25) ",
    "and are coloured by ", if (by == "celltype") "cell type" else "cell class", "; grey points are not. ",
    "The ", N_LABEL, " strongest such genes per panel are labelled, once each, at their most extreme cell type. ",
    "r, Pearson correlation over all points: ",
    paste(sprintf("%s %.2f", st$pair, st$r), collapse = ", "), "."))
  message("[fig4] ", by, " written")
}

## ---- PUBLISHED LOOK, fewer labels ---------------------------------------------------
## Exactly the published panel - every gene coloured by cell type, non-DE genes faded
## by the alpha scale, filled label boxes in the cell-type colour - with two changes:
##   * the label set is the N_PUB strongest genes (|logFC| in both) among those the
##     published rule selected, one label per gene, so ggrepel no longer drops
##     hundreds of labels at random
##   * drawn from the regional DE without donor 1040's VST cells
N_PUB <- as.integer(getarg("--n-pub", "20"))   # ~20 labels
## --n-pub takes a comma list (one set of panels per value) and --label-size sets
##   the label text (published 5). Non-default runs get a _n<N> suffix so the
##   20-label panels are not overwritten.
N_PUBS <- as.integer(strsplit(getarg("--n-pub", "20"), ",")[[1]])
LSIZE  <- as.numeric(getarg("--label-size", "5"))
## --force "Oligo,OPC" [--per-type 1]: add labels for these cell types on top of the
##   N_PUB (see pub_lab). Forced runs get a _force suffix; the default panels stay.
FORCE    <- trimws(strsplit(getarg("--force", ""), ",")[[1]])
PER_TYPE <- as.integer(getarg("--per-type", "1"))
## --genes "CD163:Micro,STC1:Endo": label these gene x cell type pairs as well
FGENES   <- trimws(strsplit(getarg("--genes", ""), ",")[[1]])
## --updown: --n-pub labels up AND --n-pub down (--n-pub 20 --updown = 20 + 20)
UPDOWN   <- "--updown" %in% args
pub_panel <- function(k) {
  d <- D[[k]]; l <- pub_lab(d, force = FORCE, per_type = PER_TYPE, genes = FGENES, updown = UPDOWN)
  d$isDE <- d$state != "neither"
  pad  <- function(v) { r <- range(v); r + c(-1, 1) * 0.04 * diff(r) }
  lims <- list(x = pad(d$logFC.x[d$isDE]), y = pad(d$logFC.y[d$isDE]))
  ggplot(d, aes(x = logFC.x, y = logFC.y, colour = cell_type, alpha = isDE)) +
    geom_point(size = 2) +
    scale_alpha_manual(values = c(`FALSE` = 0.01, `TRUE` = 0.9), guide = "none") +
    scale_colour_manual(values = cp_col, limits = names(cp_col), name = NULL, drop = FALSE) +
    ggrepel::geom_label_repel(data = l, aes(label = gene, fill = cell_type), colour = "black",
                              fontface = "bold", size = LSIZE, alpha = 1, max.overlaps = Inf,
                              box.padding = 0.5, force = 5, max.time = 5, max.iter = 1e6,
                              seed = 1, show.legend = FALSE, label.size = 0.2) +
    scale_fill_manual(values = cp_col, limits = names(cp_col), guide = "none", drop = FALSE) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.5) +
    xlab(paste0("logFC - ", d$x_reg[1])) + ylab(paste0("logFC - ", d$y_reg[1])) +
    ## Axes from the VISIBLE points only. Genes significant in neither region are
    ## drawn at alpha 0.01 - invisible, but they would still set the scale: PFC-VST
    ## is stretched to +13.7 by CFAP299 in Vip (FDR 0.27), squashing everything
    ## real into a corner. The published script handled the same problem with a
    ## hard-coded ylim(-7.5, 7.5) on the VST panels.
    coord_cartesian(xlim = lims$x, ylim = lims$y) +
    theme_minimal() +
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20),
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 6, alpha = 1)))
}
for (N_PUB in N_PUBS) for (k in names(D)) {
  g <- pub_panel(k)
  sfx <- if (N_PUB == 20 && LSIZE == 5) "" else sprintf("_n%d", N_PUB)
  if (length(FORCE)) sfx <- paste0(sfx, "_force")
  if (length(FGENES)) sfx <- paste0(sfx, "_textgenes")
  if (UPDOWN) sfx <- sprintf("_up%d_down%d", N_PUB, N_PUB)   # never overwrites the top-N panels
  f <- file.path(OUT, paste0("fig4_", sub("-", "_", k), "_published_style", sfx, ".png"))
  grDevices::png(f, width = 3500, height = 3500, res = 300)
  print(g + theme(legend.position = "none")); grDevices::dev.off()
  grDevices::cairo_pdf(sub("\\.png$", ".pdf", f), width = 11.67, height = 11.67)
  print(g + theme(legend.position = "none")); grDevices::dev.off()
  message("[fig4] published style: ", basename(f), " (",
          nrow(pub_lab(D[[k]], force = FORCE, per_type = PER_TYPE, genes = FGENES, updown = UPDOWN)), " labels)")
}
## ---- the cell-type key the published panels never showed ---------------------------
## Drawn separately rather than as a ggplot legend: a legend lists all 37 palette
## entries, including types with no point in any panel - all D1/D2 MSNs exist
## only in VST and so never have a point in a two-region comparison. This key
## shows exactly the cell types with a VISIBLE point (DE in either region) in ANY
## of the three panels, each with its dot, grouped by class, drawn at the panels'
## scale (20 pt text, as in pub_panel) and the panels' width (11.67 in).
shown <- unique(unlist(lapply(D, function(d) as.character(d$cell_type[d$state != "neither"]))))
CLASS_LAB <- c(Ext = "Excitatory", MSNs = "MSNs", Inh = "Inhibitory", Glial = "Non-neuronal")
key <- data.frame(cell_type = LEV[LEV %in% shown], stringsAsFactors = FALSE)
key$class <- factor(CLASS_LAB[COARSE[key$cell_type]], levels = CLASS_LAB)
key <- droplevels(key)
KEY_ROWS <- max(table(key$class))            # one column per class: each class reads top to bottom
key <- key %>% group_by(class) %>% mutate(row = dplyr::row_number() - 1) %>% ungroup() %>%
  mutate(col = as.integer(class) - 1)
heads <- key %>% dplyr::distinct(class, col)
kp <- ggplot(key) +
  geom_point(aes(x = col, y = -row, colour = cell_type), size = 7) +
  geom_text(aes(x = col + 0.07, y = -row, label = cell_type), hjust = 0, size = 20 / .pt) +
  geom_text(data = heads, aes(x = col - 0.04, y = 1.2, label = class), hjust = 0,
            fontface = "bold", size = 22 / .pt) +
  scale_colour_manual(values = cp_col, guide = "none") +
  scale_x_continuous(limits = c(-0.08, nlevels(key$class)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-(KEY_ROWS - 0.4), 1.8), expand = c(0, 0)) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
KEY_H <- 0.95 + 0.42 * KEY_ROWS
ggsave(file.path(OUT, "fig4_published_style_legend.png"), kp, width = 11.67, height = KEY_H, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig4_published_style_legend.pdf"), kp, width = 11.67, height = KEY_H,
       device = grDevices::cairo_pdf, bg = "white")
## ---- narrower keys ------------------------------------------------------------------
## Positions in INCHES (x and y scales = the page), so columns are exactly as wide
## as their longest label rather than a third of the page each. 20 pt text is
## ~0.14 in per character; dot + gap 0.45 in; row pitch 0.42 in.
key_inch <- function(k, cols) {
  CH <- 0.14; DOT <- 0.45; GAP <- 0.35; ROW <- 0.42; HEAD <- 0.55
  k$col <- cols[as.character(k$class)]
  cw <- tapply(seq_len(nrow(k)), k$col, function(i)
          max(DOT + CH * max(nchar(k$cell_type[i])), CH * 1.1 * nchar(as.character(k$class[i[1]]))))
  x0 <- c(0, cumsum(cw + GAP))[seq_along(cw)] + 0.1
  ## within a column, classes stack: header, then rows
  k <- k %>% group_by(col) %>% arrange(col, class) %>%
    mutate(blk = cumsum(!duplicated(class)),
           slot = dplyr::row_number() - 1 + (blk - 1) * 1.6 + blk * 0) %>% ungroup()
  k$x <- x0[k$col + 1]; k$y <- -(HEAD + ROW * k$slot + ROW * 0.5 * (k$blk - 1))
  hd <- k %>% group_by(class) %>% summarise(x = first(x), y = max(y) + HEAD * 0.9, .groups = "drop")
  ## +0.45 in on top and +0.3 on the right, so the bold headers and
  ## "Other immune" are not clipped at the edges
  W <- sum(cw) + GAP * (length(cw) - 1) + 0.55; H <- -min(k$y) + ROW * 0.8 + 0.45
  p <- ggplot(k) +
    geom_point(aes(x = x + 0.15, y = y, colour = cell_type), size = 7) +
    geom_text(aes(x = x + DOT, y = y, label = cell_type), hjust = 0, size = 20 / .pt) +
    geom_text(data = hd, aes(x = x, y = y, label = class), hjust = 0, fontface = "bold", size = 22 / .pt) +
    scale_colour_manual(values = cp_col, guide = "none") +
    scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-(H - 0.45), 0.45), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0),
                         plot.background = element_rect(fill = "white", colour = NA))
  list(plot = p, W = W, H = H)
}
k0 <- key %>% select(cell_type, class)
V <- list(
  ## the same three columns, packed to their labels
  narrow = setNames(0:2, levels(key$class)),
  ## one column: all classes stacked, to run down the side of a panel
  column = setNames(rep(0L, nlevels(key$class)), levels(key$class)))
for (nm in names(V)) {
  kk <- key_inch(k0, V[[nm]])
  f <- file.path(OUT, sprintf("fig4_published_style_legend_%s.png", nm))
  ggsave(f, kk$plot, width = kk$W, height = kk$H, dpi = 300, bg = "white")
  ggsave(sub("\\.png$", ".pdf", f), kk$plot, width = kk$W, height = kk$H,
         device = grDevices::cairo_pdf, bg = "white")
  message(sprintf("[fig4] legend (%s): %.1f x %.1f in -> %s", nm, kk$W, kk$H, basename(f)))
}
message(sprintf("[fig4] legend: %d cell types with visible points (%s), 11.67 x %.1f in",
                nrow(key), paste(sprintf("%s %d", levels(key$class), as.integer(table(key$class))),
                                 collapse = ", "), KEY_H))
message("[fig4] done")
