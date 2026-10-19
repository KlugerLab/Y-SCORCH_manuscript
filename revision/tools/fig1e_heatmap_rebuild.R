# ==============================================================================
# Fig1E terms x genes tile heatmap, rebuilt on the donor-aware astrocyte DE.
#
# The published Fig1E was a terms x genes tile heatmap, tiles coloured by fold
# change, with a -log10 adjusted p bar beside it. That layout is kept; only the
# input analyses change. ORA: enrichGO (BP) on genes with FDR < 0.05 and
# logFC > 0.5, background = all tested genes, BH-adjusted p < 0.01. Terms sharing
# >= half their genes with a more significant term are dropped; the top N terms
# by adjusted p are drawn with the N genes of largest fold change in them, plus
# SERPINA3 / OSMR / CHI3L1 where they belong to a drawn term.
#
# FOUR VERSIONS, so the layout can be compared against each source
#   ora_vs_all        enrichGO on genes up vs ALL other astrocytes (published contrast)
#   ora_vs_reactive   enrichGO on genes up vs REACTIVE astrocytes, 24 donors
#   ora_vs_reactive_pwh  the same, 20 PWH donors only
#   gsea_vs_reactive  the GSEA leading edge instead of a thresholded gene list
#
# WHY A GSEA VERSION CAN STILL BE A HEATMAP
#   GSEA reports, for each set, the genes before the running-score peak - the
#   core_enrichment / leading edge. Those are the genes carrying the signal, so
#   they take the place of "the genes in the term that passed the cutoff". The
#   figure then shows a threshold-free analysis in the published layout. Only
#   terms with NES > 0 are drawn, since the colour scale reads as up-regulation.
#
# COLOURS come from the contrast being drawn, never from a different table: a
#   tile in the vs-Reactive panels is that gene's logFC against Reactive
#   astrocytes, which is smaller than its logFC against all astrocytes (SERPINA3
#   1.93 vs all, 1.20 vs Reactive) because reactive cells already express some of
#   these genes. Mixing the two would overstate the vs-Reactive panels.
#
# INPUT   OUT_ROOT/06-astrocyte_analysis/astro_Reactive_PWH_donorDE_no1040vst.csv,
#         astro_Reactive_PWH_vs_Reactive_donorDE_no1040vst.csv (and pwh_only/),
#         figures/Fig1E_gsea/gsea_vs_reactive_24_primary.csv
# OUTPUT  OUT_ROOT/06-astrocyte_analysis/figures/Fig1E_heatmaps/
#   Fig1E_<source>_heatmap_no1040vst.png/.pdf + _caption.txt
#   Fig1E_<source>_terms.csv      the drawn terms and their genes
#
# Usage: revision/docker/run.sh revision/tools/fig1e_heatmap_rebuild.R
#          [--n-terms 10] [--n-genes 40] [--width 7.2] [--title "..."]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
TOOLS <- file.path(SCORCH_REPO, "revision/tools")
for (f in c("astro_pub_theme", "astro_cluster_features")) source(file.path(TOOLS, paste0(f, ".R")))
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(patchwork)
})
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }

D    <- file.path(OUT_ROOT, "06-astrocyte_analysis")
OUT  <- file.path(D, "figures", "Fig1E_heatmaps")
GSEA <- file.path(D, "figures", "Fig1E_gsea")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SFX     <- "_no1040vst"
PADJ    <- 0.01
N_TERMS <- as.integer(getarg("--n-terms", 10))
N_SHOW  <- as.integer(getarg("--n-genes", 40))
## figure width in inches
FIG_W  <- as.numeric(getarg("--width", 7.2))
TITLE_OVERRIDE <- getarg("--title", "")
MARKERS <- c("SERPINA3", "OSMR", "CHI3L1")

DE_FILES <- list(
  vs_all          = file.path(D, "astro_Reactive_PWH_donorDE_no1040vst.csv"),
  vs_reactive     = file.path(D, "astro_Reactive_PWH_vs_Reactive_donorDE_no1040vst.csv"),
  vs_reactive_pwh = file.path(D, "pwh_only", "astro_Reactive_PWH_vs_Reactive_donorDE_no1040vst.csv"))
DE <- lapply(DE_FILES, function(f) if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL)

## the Fig1E gene-overlap rule: drop a term sharing >= half its genes with a
## more significant term already kept
overlap_prune <- function(r) {
  r  <- r[order(r$p.adjust), ]
  gs <- strsplit(r$geneID, "/", fixed = TRUE)
  keep <- logical(nrow(r))
  for (i in seq_len(nrow(r)))
    keep[i] <- !any(vapply(which(keep), function(j)
      length(intersect(gs[[i]], gs[[j]])) / length(union(gs[[i]], gs[[j]])), numeric(1)) >= 0.5)
  r[keep, ]
}

# ------------------------------------------------------------------------------
#' The published Fig1E layout. `r` needs Description / geneID / p.adjust; `de`
#' supplies the colours and the gene ranking.
draw <- function(r, de, tag, fill_lab, caption, width = FIG_W, title = NULL,
                 term_pt = 7.3, gene_pt = 6.3, wrap_n = 34, leg_pt = 6.5) {
  r <- overlap_prune(r)
  r$label <- r$Description
  tr <- utils::head(r[order(r$p.adjust), ], N_TERMS)
  if (!nrow(tr)) { message("[hm] ", tag, ": no terms to draw"); return(invisible(NULL)) }

  in_terms <- unique(unlist(strsplit(tr$geneID, "/", fixed = TRUE)))
  d_in <- de[de$gene %in% in_terms, ]
  show <- utils::head(d_in$gene[order(-d_in$logFC)], N_SHOW)
  forced <- setdiff(intersect(MARKERS, in_terms), show); show <- c(show, forced)
  message("[hm] ", tag, ": markers in drawn terms: ",
          paste(intersect(MARKERS, in_terms), collapse = ", "),
          if (length(forced)) paste0(" (added: ", paste(forced, collapse = ", "), ")") else "",
          "; absent from every drawn term: ",
          paste(setdiff(MARKERS, in_terms), collapse = ", "))

  long <- do.call(rbind, lapply(seq_len(nrow(tr)), function(i) {
    g <- intersect(strsplit(tr$geneID[i], "/", fixed = TRUE)[[1]], show)
    if (length(g)) data.frame(term = tr$label[i], gene = g, stringsAsFactors = FALSE) else NULL }))
  long$logFC <- de$logFC[match(long$gene, de$gene)]
  first <- tapply(match(long$term, tr$label), long$gene, min)
  lfc   <- tapply(long$logFC, long$gene, max)
  gene_order <- names(first)[order(first, -lfc[names(first)])]
  wrap <- function(x) vapply(x, function(s) paste(strwrap(s, wrap_n), collapse = "\n"), character(1))
  term_lv <- rev(wrap(tr$label))
  long$term <- factor(wrap(long$term), levels = term_lv)
  long$gene <- factor(long$gene, levels = gene_order)
  bars <- data.frame(term = factor(wrap(tr$label), levels = term_lv), mlp = -log10(tr$p.adjust))

  p_hm <- ggplot(long, aes(x = gene, y = term, fill = logFC)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_gradient(low = "#FCDCCF", high = "#B2182B", name = fill_lab) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) + theme_astro_pub(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = gene_pt, face = "italic"),
          axis.text.y = element_text(size = term_pt, lineheight = 0.85),
          panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.3),
          legend.position = "right", legend.key.height = unit(0.45, "cm"),
          legend.key.width = unit(0.25, "cm"), legend.title = element_text(size = leg_pt),
          legend.text = element_text(size = leg_pt))
  p_bar <- ggplot(bars, aes(x = mlp, y = term)) + geom_col(fill = "grey55", width = 0.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08)), n.breaks = 3) +
    labs(x = expression(-log[10]~adj.~p), y = NULL) + theme_astro_pub(base_size = 7) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.title.x = element_text(size = 6.5))

  f <- file.path(OUT, paste0("Fig1E_", tag, "_heatmap", SFX, ".png"))
  max_lines <- max(lengths(strsplit(term_lv, "\n")))
  row_in <- max(0.26, 0.05 + max_lines * term_pt / 72 * 0.95)
  p <- p_hm + p_bar + plot_layout(widths = c(1, 0.12), guides = "collect")
  ## bold, centred title; --title overrides it for the vs-Reactive panel
  if (!is.null(title))
    p <- p + plot_annotation(title = title, theme = theme(
      plot.title = element_text(face = "bold", size = term_pt + 2, hjust = 0.5)))
  save_pub(p, f, width = width, height = 1.1 + row_in * nrow(tr) + if (is.null(title)) 0 else 0.3)
  astro_caption(f, caption)
  utils::write.csv(tr[, intersect(c("Description", "p.adjust", "NES", "setSize", "Count", "geneID"),
                                  names(tr))],
                   file.path(OUT, paste0("Fig1E_", tag, "_terms.csv")), row.names = FALSE)
  message("[hm] ", tag, ": drew ", nrow(tr), " terms, ", length(show), " genes -> ", basename(f))
}

ora_of <- function(de) {
  g <- de$gene[de$FDR < 0.05 & de$logFC > 0.5]
  e <- enrichGO(gene = g, universe = unique(de$gene), keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = PADJ, qvalueCutoff = 0.05,
                readable = TRUE)
  list(n = length(g), res = as.data.frame(e))
}

## ---- 1-3. the thresholded versions, in the published layout ------------------------------
SRC <- list(
  list(tag = "ora_vs_all", de = "vs_all",
       title = "GO enrichment: Reactive-PWH vs other astrocytes",
       what = "all other astrocytes", extra = "the contrast used in the published figure"),
  list(tag = "ora_vs_reactive", de = "vs_reactive",
       title = "GO enrichment: Reactive-PWH vs Reactive astrocytes",
       what = "Reactive astrocytes (24 donors)",
       extra = "Reactive astrocytes are the comparator, so generic reactivity is subtracted"),
  list(tag = "ora_vs_reactive_pwh", de = "vs_reactive_pwh",
       title = "GO enrichment: Reactive-PWH vs Reactive astrocytes (PWH donors)",
       what = "Reactive astrocytes (20 PWH donors)",
       extra = "PWH donors only; the 24-donor version is the sensitivity analysis"))

for (s in SRC) {
  de <- DE[[s$de]]
  if (is.null(de)) { message("[hm] missing DE table for ", s$tag); next }
  o <- ora_of(de)
  message("\n[hm] === ", s$tag, ": ", o$n, " up genes, ", nrow(o$res), " terms ===")
  draw(o$res, de, s$tag,
       title = if (nzchar(TITLE_OVERRIDE) && s$tag == "ora_vs_reactive") TITLE_OVERRIDE else s$title,
       fill_lab = paste0("log2 FC\n(Reactive-PWH vs\n", sub(" \\(.*", "", s$what), ")"),
       caption = paste0(
         "GO biological process terms enriched among genes up in Reactive-PWH astrocytes relative to ",
         s$what, " (donor-paired pseudobulk, edgeR quasi-likelihood with donor as a blocking factor; ",
         "FDR < 0.05, log2 fold change > 0.5, n = ", o$n, " genes); background, the ", nrow(de),
         " genes tested in that contrast; BH-adjusted p < ", PADJ, ". ", s$extra,
         ". Terms sharing at least half their genes with a more significant term are removed. Top ",
         N_TERMS, " by adjusted p; tiles show the ", N_SHOW, " genes with the largest fold change in ",
         "those terms, plus SERPINA3, OSMR and CHI3L1 where they belong to a drawn term; ",
         "bars, -log10 adjusted p."))
}

## ---- 4. the same layout, from the GSEA leading edge ---------------------------------------
gf <- file.path(GSEA, "gsea_vs_reactive_24_primary.csv")
if (file.exists(gf)) {
  g <- utils::read.csv(gf, stringsAsFactors = FALSE)
  g <- g[g$NES > 0, ]                       # the colour scale reads as up-regulation
  g$geneID <- g$core_enrichment             # leading edge takes the place of the gene list
  message("\n[hm] === gsea_vs_reactive: ", nrow(g), " up sets ===")
  draw(g, DE$vs_reactive, "gsea_vs_reactive",
       title = "GSEA: Reactive-PWH vs Reactive astrocytes",
       fill_lab = "log2 FC\n(Reactive-PWH vs\nReactive)",
       caption = paste0(
         "GO biological process sets enriched in Reactive-PWH astrocytes relative to Reactive ",
         "astrocytes by gene set enrichment analysis (24 donors; all ", nrow(DE$vs_reactive),
         " tested genes ranked by sign(log2 fold change) x -log10 p, no significance or fold-change ",
         "threshold; BH-adjusted p < 0.05). Tiles show the leading-edge genes - those contributing to ",
         "the enrichment score - for the top ", N_TERMS, " sets with a positive normalised enrichment ",
         "score, limited to the ", N_SHOW, " with the largest fold change, plus SERPINA3, OSMR and ",
         "CHI3L1 where they are in a drawn set; sets sharing at least half their leading-edge genes ",
         "with a more significant set are removed; bars, -log10 adjusted p."))
} else message("[hm] no GSEA table at ", gf, " - the GSEA panel needs a ranked GSEA table (not part of this release)")

message("\n[hm] wrote ", OUT)
