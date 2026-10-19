# ==============================================================================
# Figure 3A: GO pathways from the CROSS-REGION DE, one compact panel for every
# cell type.
#
# WHY IT IS SMALLER  The per-region panels need ~24 in of height for three classes
# and repeat each pathway three times, once
# per region, and most of those dots say the same thing. The cross-region DE has
# no region dimension, so cell type is the only x axis: glia, excitatory and
# inhibitory all fit in one panel.
#
# SAME RULES AS revision/tools/fig3_pathway_panels.R
#   DE source   --de edgeR (default; cross-region edgeR table) or donor_aware
#               (limma-voom + duplicateCorrelation, block = donor)
#   DE genes    p_adj.glb < 0.05, |logFC| > 0.25, detected in > 5% of nuclei
#   enrichment  compareCluster / enrichGO BP, BH, no significance cutoff at the
#               cache stage (selection applies it), readable symbols
#   redundancy  GOSemSim Wang similarity, per cell type over its best 20 terms,
#               then again across the panel's selected terms
#   selection   top N per cell type per direction, then every cell type's
#               measurement of a selected term is drawn
#
# LAYOUT  rows = pathways, up block above down block; columns = cell types
#   grouped into excitatory / inhibitory / glial; colour = mean logFC (diverging,
#   so up reads red and down blue); size = -log10 adjusted p.
#
# INPUT   OUT_ROOT/03-DE_analysis/cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#         or OUT_ROOT/03-DE_analysis/cross_region/donor_aware/cross_region_donor_aware_DE.csv
# OUTPUT  OUT_ROOT/03-DE_analysis/figure_prototypes/fig3_crossregion/
#           fig3_crossregion_pathways_top<N>[_donor_aware].png/.pdf
#           fig3_crossregion_pathways_top<N>[_donor_aware]_plotted.csv   every dot drawn
#           fig3_crossregion_enrichment_<de>.rds                         cached enrichment
#
# Usage: revision/docker/run.sh revision/tools/fig3_crossregion_pathways.R [--top 2]
#          [--de edgeR|donor_aware] [--min-ct 2] [--simplify 0.6] [--refresh]
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
  library(clusterProfiler); library(org.Hs.eg.db); library(Homo.sapiens)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
TOP_N    <- as.integer(getarg("--top", "2"))
DE_SRC   <- match.arg(getarg("--de", "edgeR"), c("edgeR", "donor_aware"))
MIN_CT   <- as.integer(getarg("--min-ct", "2"))
SIMPLIFY <- as.numeric(getarg("--simplify", "0.6"))
SIMPLIFY_POOL <- as.integer(getarg("--simplify-pool", "20"))
MIN_SEL_GENES  <- as.integer(getarg("--min-select-genes", "3"))
MIN_DRAW_GENES <- as.integer(getarg("--min-draw-genes", "1"))
REFRESH  <- "--refresh" %in% args
OUT <- file.path(OUT_ROOT, "03-DE_analysis", "figure_prototypes", "fig3_crossregion")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path(OUT, sprintf("fig3_crossregion_enrichment_%s.rds", DE_SRC))
## figure/table names carry the DE source, so a donor_aware run does not
## overwrite the edgeR panel
SRC_SFX <- if (DE_SRC == "donor_aware") "_donor_aware" else ""

MSN_SUBTYPES <- c("D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP",
                  "D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid")
CLASS <- list(
  Excitatory = c("L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP","L6 IT","L6 IT Car3","L6b","L6 CT"),
  Inhibitory = c("Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                 "MSNs", MSN_SUBTYPES),
  Glial      = c("Astro","Oligo","OPC","Micro","Macrophages","Immune","Other immune","Endo",
                 "PC","SMC","VLMC"))
CT_ORDER <- unlist(CLASS, use.names = FALSE)
class_of <- function(ct) { for (k in names(CLASS)) if (ct %in% CLASS[[k]]) return(k); NA_character_ }

## ---- DE input ----------------------------------------------------------------
load_de <- function() {
  if (DE_SRC == "edgeR") {
    f <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "data",
                   "all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
    readRDS(f) %>% ungroup() %>%
      dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05)) %>%
      dplyr::transmute(gene, cell_type = as.character(cluster_id), logFC, p_adj.glb)
  } else {
    ## donor-aware limma-voom + duplicateCorrelation (revision/tools/cross_region_donor_aware.R);
    ## it carries no detection frequency, the DE table is already filtered upstream
    utils::read.csv(file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "donor_aware",
                              "cross_region_donor_aware_DE.csv"), stringsAsFactors = FALSE) %>%
      dplyr::transmute(gene, cell_type = as.character(cluster_id), logFC, p_adj.glb)
  }
}

compute_enrichment <- function() {
  de <- load_de()
  de$class <- vapply(de$cell_type, class_of, character(1))
  de <- de[!is.na(de$class), ]
  message("[xr-path] DE rows: ", nrow(de), " over ", dplyr::n_distinct(de$cell_type), " cell types (", DE_SRC, ")")

  out <- list()
  for (dirn in c("up", "down")) {
    sig <- de %>% dplyr::filter(p_adj.glb < 0.05, abs(logFC) > 0.25,
                                if (dirn == "up") logFC > 0 else logFC < 0)
    message("[xr-path] ", dirn, ": ", nrow(sig), " significant gene-cell type pairs")
    gene_lists <- list()
    for (cp in unique(sig$cell_type)) {
      sub <- sig[sig$cell_type == cp, ]
      if (nrow(sub) <= 1) next
      conv <- suppressMessages(AnnotationDbi::select(Homo.sapiens, keys = unique(sub$gene),
                                                     keytype = "SYMBOL", columns = "ENTREZID"))
      ent <- conv$ENTREZID[!is.na(conv$ENTREZID)]
      if (length(ent) > 1) gene_lists[[cp]] <- ent
    }
    message("[xr-path] ", dirn, ": ", length(gene_lists), " cell types with a gene list")
    if (!length(gene_lists)) next
    ego <- compareCluster(geneClusters = gene_lists, pAdjustMethod = "BH",
                          pvalueCutoff = 1, qvalueCutoff = 1, keyType = "ENTREZID",
                          OrgDb = org.Hs.eg.db, readable = TRUE, ont = "BP", fun = "enrichGO")
    long <- as.data.frame(ego@compareClusterResult) %>%
      dplyr::rename(cell_type = Cluster) %>%
      separate_longer_delim(geneID, delim = "/") %>%
      mutate(gene = geneID, pathway = Description, cell_type = as.character(cell_type)) %>%
      left_join(sig %>% dplyr::select(cell_type, gene, logFC), by = c("cell_type", "gene")) %>%
      group_by(cell_type, pathway, ID, p.adjust, pvalue) %>%
      summarise(avr_logFC = mean(logFC, na.rm = TRUE), n_genes = dplyr::n(), .groups = "drop") %>%
      mutate(direction = dirn)
    out[[dirn]] <- long
  }
  long <- bind_rows(out)
  saveRDS(long, CACHE)
  message("[xr-path] cached -> ", CACHE)
  long
}

long <- if (file.exists(CACHE) && !REFRESH) {
  message("[xr-path] using cached enrichment ", CACHE, " (--refresh to redo)"); readRDS(CACHE)
} else compute_enrichment()
long$class <- vapply(long$cell_type, class_of, character(1))

## ---- term selection, as the per-region panels --------------------------------
GO_SEM <- if (SIMPLIFY > 0)
  GOSemSim::godata("org.Hs.eg.db", ont = "BP", computeIC = FALSE) else NULL
simplify_ids <- function(ids, score, cutoff) {
  if (length(ids) < 2) return(ids)
  sim <- GOSemSim::termSim(ids, ids, semData = GO_SEM, method = "Wang")
  sim[is.na(sim)] <- 0
  remove <- character(0)
  for (i in seq_along(ids)) {
    ii <- which(sim[, i] > cutoff)
    if (length(ii) < 2) next
    best <- ii[score[ii] == min(score[ii])]
    remove <- union(remove, ids[setdiff(ii, best)])
  }
  setdiff(ids, remove)
}

sel_of <- function(dirn) {
  d  <- long[long$direction == dirn & long$n_genes >= MIN_DRAW_GENES, ]
  ds <- d[d$n_genes >= MIN_SEL_GENES & d$p.adjust < 0.05, ]
  if (!nrow(ds)) return(NULL)
  if (SIMPLIFY > 0) {
    sel <- ds %>% group_by(cell_type) %>% slice_min(pvalue, n = SIMPLIFY_POOL) %>%
      group_modify(function(g, k) dplyr::slice_min(g[g$ID %in% simplify_ids(g$ID, g$pvalue, SIMPLIFY), ],
                                                   pvalue, n = TOP_N)) %>% ungroup()
    keep <- simplify_ids(unique(sel$ID), sel$pvalue[match(unique(sel$ID), sel$ID)], SIMPLIFY)
  } else {
    sel  <- ds %>% group_by(cell_type) %>% slice_min(pvalue, n = TOP_N) %>% ungroup()
    keep <- unique(sel$ID)
  }
  ## drop terms only one cell type ever picked up, as the per-region panels do
  n_ct <- d %>% dplyr::filter(ID %in% keep, p.adjust < 0.05) %>%
    group_by(ID) %>% summarise(n = dplyr::n_distinct(cell_type), .groups = "drop")
  keep <- n_ct$ID[n_ct$n >= MIN_CT]
  d[d$ID %in% keep, ]
}
plotted <- bind_rows(lapply(c("up", "down"), sel_of))
if (!nrow(plotted)) stop("no terms selected")
DIR_LAB <- c(up = "Up in PWH", down = "Down in PWH")
plotted$direction <- factor(DIR_LAB[plotted$direction], levels = unname(DIR_LAB))
plotted$class     <- factor(plotted$class, levels = names(CLASS))
plotted$cell_type <- factor(plotted$cell_type, levels = intersect(CT_ORDER, unique(plotted$cell_type)))
ord <- plotted %>% group_by(direction, pathway) %>% summarise(best = min(p.adjust), .groups = "drop") %>%
  arrange(direction, best)
## same GO shorthand as revision/tools/fig3_updown_panels.R, so most labels stay on one line
abbreviate_go <- function(x) {
  rules <- c("positive regulation of " = "pos. reg. of ", "negative regulation of " = "neg. reg. of ",
             "regulation of " = "reg. of ", "biosynthetic process" = "biosynthesis",
             "catabolic process" = "catabolism", "metabolic process" = "metabolism",
             "signaling pathway" = "signaling", "signal transduction" = "signaling",
             "energy derivation by oxidation of organic compounds" = "energy derivation (oxidation)",
             "cellular response to " = "response to ")
  for (k in names(rules)) x <- gsub(k, rules[[k]], x, fixed = TRUE)
  x
}
plotted$label <- factor(stringr::str_wrap(abbreviate_go(plotted$pathway), 44),
                        levels = unique(stringr::str_wrap(abbreviate_go(ord$pathway), 44)))
utils::write.csv(plotted %>% dplyr::select(direction, class, cell_type, pathway, ID, p.adjust,
                                           pvalue, avr_logFC, n_genes),
                 file.path(OUT, sprintf("fig3_crossregion_pathways_top%d%s_plotted.csv", TOP_N, SRC_SFX)),
                 row.names = FALSE)
message("[xr-path] terms drawn: ",
        paste(sprintf("%s %d", levels(plotted$direction),
                      tapply(plotted$ID, plotted$direction, dplyr::n_distinct)), collapse = ", "),
        " over ", dplyr::n_distinct(plotted$cell_type), " cell types")

COL_LIM <- as.numeric(stats::quantile(abs(plotted$avr_logFC), 0.95))
message(sprintf("[xr-path] colour capped at +/-%.2f (max |mean logFC| %.2f)",
                COL_LIM, max(abs(plotted$avr_logFC))))
p <- ggplot(plotted, aes(x = cell_type, y = label, size = -log10(p.adjust), colour = avr_logFC)) +
  geom_point(na.rm = TRUE) +
  scale_size(range = c(0.9, 4.2), name = expression(-log[10]~adj.~p)) +
  ## a few terms carry a mean logFC of ~5 and flattened everything else to pale;
  ## the scale is capped at the 95th percentile of |mean logFC| and values beyond
  ## it are squished onto the end colour
  scale_colour_gradient2(low = "#2166AC", mid = "grey93", high = "#B2182B", midpoint = 0,
                         limits = c(-COL_LIM, COL_LIM), oob = scales::squish,
                         name = "mean logFC") +
  facet_grid(direction ~ class, scales = "free", space = "free", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_cowplot(font_size = 9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(colour = "grey78", fill = NA, linewidth = 0.5),
        panel.spacing.x = unit(0.35, "lines"), panel.spacing.y = unit(1.0, "lines"),
        strip.text.y.left = element_text(face = "bold", size = 9, angle = 90),
        strip.text.x = element_text(face = "bold", size = 9.5),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
        axis.text.y = element_text(size = 9, lineheight = 0.85),
        legend.position = "right", legend.title = element_text(size = 8.5, face = "bold"),
        legend.text = element_text(size = 7.5), legend.key.height = unit(0.75, "lines"),
        plot.margin = margin(6, 10, 4, 6))

n_path  <- dplyr::n_distinct(plotted$label)
n_cells <- dplyr::n_distinct(plotted$cell_type)
v_lines <- lengths(regmatches(levels(plotted$label), gregexpr("\n", levels(plotted$label)))) + 1
q_lines <- max(1, as.numeric(stats::quantile(v_lines, 0.9)))
w <- 2.5 + 1.4 + 0.30 * n_cells
h <- 1.05 + 0.185 * q_lines * n_path + 0.55
f <- file.path(OUT, sprintf("fig3_crossregion_pathways_top%d%s.png", TOP_N, SRC_SFX))
ggsave(f, p, width = w, height = h, dpi = 300, limitsize = FALSE, bg = "white")
ggsave(sub("\\.png$", ".pdf", f), p, width = w, height = h, device = grDevices::cairo_pdf, bg = "white")
message(sprintf("[xr-path] %d pathways x %d cell types -> %.1f x %.1f in  %s", n_path, n_cells, w, h, basename(f)))
message("[xr-path] done")
