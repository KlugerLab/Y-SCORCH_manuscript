# ==============================================================================
# Robustness of the Reactive-PWH astrocyte population to integration and
# clustering parameter choices (reviewer request).
#
# Method
#   Reactive-PWH was defined at resolution 1 as the union of four clusters (6, 12,
#   14, 19 in astrocyte_analysis.Rmd), so the share of its nuclei in one new cluster
#   penalises a state splitting into sub-clusters, which is not the question. The
#   question is whether the Reactive-PWH nuclei stay together, apart from the other
#   subtypes. At every setting:
#     majority clusters  new clusters in which > 50% of nuclei are Reactive-PWH
#     recall             Reactive-PWH nuclei inside those clusters / all Reactive-PWH
#     precision          Reactive-PWH nuclei inside them / all nuclei inside them
#     jaccard            |Reactive-PWH n union of majority clusters| / |union of both|
#   Because the claim is a PWH-enriched state, the recovered population (union of
#   majority clusters) goes through the same donor-level composition test as the
#   original label (revision/tools/astro_composition_donor.R: quasibinomial
#   ~ HIV + Age + Sex + Race + PMI, one row per donor). The single-cluster numbers
#   are also kept in the table (pct_in_best, best_jaccard). The same majority-cluster
#   Jaccard is computed for every subtype as a benchmark, and the per-cluster share
#   of the largest donor / region shows what clusters follow when a setting fails.
#
#   Settings (cells = the 7 astrocyte subtypes of astro_all_clusters_labeled_no1040vst.RDS,
#   Oligo and Endo-Astro removed; Louvain, default seed):
#   published        harmony_astro (Harmony over libraries), 30 PCs, k = 20, as in the Rmd
#   no_integration   pca, 30 PCs, k = 20
#   harmony_donor    Harmony re-run grouping by donor (not donor x region), 30 PCs, k = 20
#   pcs20 / pcs50    harmony_astro, 20 or 50 PCs, k = 20
#   k10 / k50        harmony_astro, 30 PCs, k = 10 or 50
#   each at resolution 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0
#
# Inputs   OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled_no1040vst.RDS
# Outputs  OUT_ROOT/06-astrocyte_analysis/figures/ask3_robustness/
#            clusters_<setting>.rds   cell x resolution assignments (metrics can be
#                                     recomputed without re-clustering)
#            robustness_metrics.csv, subtype_benchmark.csv, integration_probe.csv
#            ask3_robustness_Reactive_PWH.png / .pdf / _caption.txt
#
# Usage: revision/docker/run.sh revision/tools/astro_cluster_robustness.R [--replot]
#        (--replot redraws from the saved tables)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(harmony); library(ggplot2); library(dplyr)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
TOOLS <- file.path(SCORCH_REPO, "revision/tools")
source(file.path(TOOLS, "astro_pub_theme.R"))
source(file.path(TOOLS, "astro_composition_donor.R"))

args <- commandArgs(trailingOnly = TRUE)
REPLOT  <- "--replot" %in% args
TARGET  <- "Reactive-PWH"
LAB     <- "subcluster_labels"
DROP    <- c("Oligo", "Endo-Astro")
RES     <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0)
OUT <- file.path(OUT_ROOT, "06-astrocyte_analysis", "figures", "ask3_robustness")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time(); el <- function() format(round(difftime(Sys.time(), t0, units = "mins"), 1))

SETTINGS <- data.frame(
  setting   = c("published", "no_integration", "harmony_donor", "pcs20", "pcs50", "k10", "k50"),
  reduction = c("harmony_astro", "pca", "harmony_donor", "harmony_astro", "harmony_astro", "harmony_astro", "harmony_astro"),
  npcs      = c(30, 30, 30, 20, 50, 30, 30),
  k         = c(20, 20, 20, 20, 20, 10, 50),
  stringsAsFactors = FALSE)

## --replot redraws from the saved tables and does not load the object when they exist
MF <- file.path(OUT, "robustness_metrics.csv")
BF <- file.path(OUT, "subtype_benchmark.csv")
IF <- file.path(OUT, "integration_probe.csv")
FROM_TABLES <- REPLOT && all(file.exists(c(MF, BF, IF)))
if (!FROM_TABLES) {
  obj <- readRDS(file.path(OUT_ROOT, "06-astrocyte_analysis", "data", "astro_all_clusters_labeled_no1040vst.RDS"))
  Idents(obj) <- LAB
  obj <- subset(obj, idents = setdiff(unique(as.character(obj@meta.data[[LAB]])), DROP))
  lab <- setNames(as.character(obj@meta.data[[LAB]]), colnames(obj))
  target <- names(lab)[lab == TARGET]
  message(sprintf("[rob] %d nuclei in 7 subtypes, %d %s (%s min); reductions: %s", ncol(obj), length(target),
                  TARGET, el(), paste(names(obj@reductions), collapse = ", ")))
  for (r in c("harmony_astro", "pca"))
    if (!r %in% names(obj@reductions)) stop("reduction '", r, "' is not on the object")
}

## ---- cluster every setting (skipped for settings already saved) -------------------
if (!REPLOT) for (i in seq_len(nrow(SETTINGS))) {
  s <- SETTINGS[i, ]
  f <- file.path(OUT, sprintf("clusters_%s.rds", s$setting))
  if (file.exists(f)) { message("[rob] ", s$setting, ": saved - skipped"); next }
  if (s$reduction == "harmony_donor" && !"harmony_donor" %in% names(obj@reductions)) {
    set.seed(1)
    obj <- RunHarmony(obj, group.by.vars = "sample_ids", reduction.use = "pca",
                      dims.use = seq_len(ncol(Embeddings(obj, "pca"))),
                      reduction.save = "harmony_donor", verbose = FALSE)
    message("[rob] Harmony by donor done (", el(), " min)")
  }
  nd <- min(s$npcs, ncol(Embeddings(obj, s$reduction)))
  if (nd < s$npcs) message("[rob] ", s$reduction, " has only ", nd, " dims")
  g <- paste0("g_", s$setting)
  w <- FindNeighbors(obj, reduction = s$reduction, dims = seq_len(nd), k.param = s$k,
                     graph.name = c(paste0(g, "_nn"), paste0(g, "_snn")), verbose = FALSE)
  cl <- sapply(RES, function(r)
    as.character(FindClusters(w, graph.name = paste0(g, "_snn"), resolution = r, verbose = FALSE)$seurat_clusters))
  dimnames(cl) <- list(colnames(obj), paste0("res_", RES))
  saveRDS(cl, f)
  message(sprintf("[rob] %s: clustered at %d resolutions (%s min)", s$setting, length(RES), el()))
  rm(w); invisible(gc())
}

## ---- metrics ----------------------------------------------------------------------
if (FROM_TABLES) {
  M <- utils::read.csv(MF, stringsAsFactors = FALSE)
  B <- utils::read.csv(BF, stringsAsFactors = FALSE)
  I <- utils::read.csv(IF, stringsAsFactors = FALSE)
  message("[rob] --replot: metrics read from the saved tables")
} else {
comp_test <- function(recovered) {
  tmp <- obj
  tmp$recovered <- ifelse(colnames(tmp) %in% recovered, "recovered", "other")
  r <- tryCatch(astro_composition(tmp, label_col = "recovered", outdir = NULL)$results,
                error = function(e) { message("[rob]   composition test failed: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(c(or = NA, p = NA, pct_PWoH = NA, pct_PWH = NA))
  r <- r[r$subtype == "recovered", ]
  c(or = r$odds_ratio, p = r$p, pct_PWoH = r$mean_pct_PWoH, pct_PWH = r$mean_pct_PWH)
}
M <- do.call(rbind, lapply(SETTINGS$setting, function(st) {
  f <- file.path(OUT, sprintf("clusters_%s.rds", st))
  if (!file.exists(f)) { message("[rob] ", st, ": no clusters saved"); return(NULL) }
  cl <- readRDS(f)
  do.call(rbind, lapply(seq_along(RES), function(j) {
    v <- setNames(cl[, j], rownames(cl))
    frac <- tapply(names(v) %in% target, v, mean)
    maj  <- names(frac)[frac > 0.5]
    uni  <- names(v)[v %in% maj]
    hit  <- sum(uni %in% target)
    tb   <- table(v[target]); best <- names(tb)[which.max(tb)]
    ct   <- if (length(uni)) comp_test(uni) else c(or = NA, p = NA, pct_PWoH = NA, pct_PWH = NA)
    data.frame(setting = st, resolution = RES[j], n_clusters = length(unique(v)),
               n_majority_clusters = length(maj),
               recall = round(hit / length(target), 3),
               precision = if (length(uni)) round(hit / length(uni), 3) else NA,
               jaccard = round(hit / length(union(uni, target)), 3),
               pct_in_best = round(100 * max(tb) / length(target), 1),
               best_jaccard = round(max(tb) / (length(target) + sum(v == best) - max(tb)), 3),
               recovered_nuclei = length(uni),
               odds_ratio = unname(ct["or"]), p = unname(ct["p"]),
               mean_pct_PWoH = unname(ct["pct_PWoH"]), mean_pct_PWH = unname(ct["pct_PWH"]),
               stringsAsFactors = FALSE)
  }))
}))
utils::write.csv(M, MF, row.names = FALSE)
print(M[, c("setting", "resolution", "n_clusters", "recall", "precision", "jaccard", "odds_ratio", "p")], row.names = FALSE)

## BENCHMARK: the same majority-cluster Jaccard for every subtype,
##   so Reactive-PWH is judged against the established subtypes, not against 1.
## LIBRARY DOMINANCE: nucleus-weighted share of each new cluster coming from its
##   largest donor / region - shows what the clusters follow when a setting fails.
md <- obj@meta.data
B <- list(); I <- list()
for (st in SETTINGS$setting) {
  f <- file.path(OUT, sprintf("clusters_%s.rds", st)); if (!file.exists(f)) next
  cl <- readRDS(f)
  lb <- as.character(md[rownames(cl), LAB]); dn <- as.character(md[rownames(cl), "sample_ids"])
  rg <- toupper(as.character(md[rownames(cl), "region"]))
  for (j in seq_len(ncol(cl))) {
    v <- cl[, j]
    for (s in unique(lb)) {
      tg <- lb == s; frac <- tapply(tg, v, mean); uni <- v %in% names(frac)[frac > 0.5]
      B[[length(B) + 1]] <- data.frame(setting = st, resolution = RES[j], subtype = s, n = sum(tg),
                                       jaccard = round(sum(tg & uni) / sum(tg | uni), 3))
    }
    top <- function(g) { t <- table(v, g); sum(apply(t, 1, max)) / length(v) }
    I[[length(I) + 1]] <- data.frame(setting = st, resolution = RES[j], n_clusters = length(unique(v)),
      donor_top_share_weighted = round(top(dn), 3), region_top_share_weighted = round(top(rg), 3),
      clusters_one_donor_gt80pct = sum(apply(table(v, dn), 1, function(x) max(x) / sum(x)) > 0.8))
  }
}
B <- do.call(rbind, B); I <- do.call(rbind, I)
utils::write.csv(B, BF, row.names = FALSE); utils::write.csv(I, IF, row.names = FALSE)
}

## ---- figure -----------------------------------------------------------------------
## The no-integration odds ratio is not drawn: its "recovered" population holds ~10%
##   of the Reactive-PWH nuclei and the control rate is ~0, so the quasibinomial fit
##   degenerates (p = 1.000). Panel c is the benchmark heatmap of every subtype.
LABS <- c(published = "Published (Harmony by library, 30 PCs, k = 20)", pcs20 = "20 PCs", pcs50 = "50 PCs",
          k10 = "k = 10", k50 = "k = 50", harmony_donor = "Harmony by donor", no_integration = "No integration")
SHORT <- c(published = "Published", pcs20 = "20 PCs", pcs50 = "50 PCs", k10 = "k = 10", k50 = "k = 50",
           harmony_donor = "Harmony by donor", no_integration = "No integration")
COLS <- c(published = "black", pcs20 = "#67A9CF", pcs50 = "#2166AC", k10 = "#9E9AC8", k50 = "#54278F",
          harmony_donor = "#EF8A62", no_integration = "#B2182B")
P <- M; P$setting <- factor(P$setting, levels = names(LABS))
P$sig <- ifelse(is.finite(P$p) & P$p < 0.05, "p < 0.05", "p ≥ 0.05")
vline <- geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3)
p1 <- ggplot(P, aes(resolution, jaccard, colour = setting, group = setting)) + vline +
  geom_line(linewidth = 0.5) + geom_point(size = 1.3) +
  scale_colour_manual(values = COLS, labels = LABS, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = RES) +
  scale_y_continuous(limits = c(0, 1), name = "Recovery of Reactive-PWH\n(Jaccard)") + xlab(NULL) +
  theme_astro_pub() + guides(colour = guide_legend(ncol = 2))
Q <- P[P$setting != "no_integration" & is.finite(P$odds_ratio), ]
p2 <- ggplot(Q, aes(resolution, odds_ratio, colour = setting, group = setting)) + vline +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.5) + geom_point(aes(shape = sig), size = 1.5, fill = "white") +
  scale_colour_manual(values = COLS, guide = "none") +
  scale_shape_manual(values = c(`p < 0.05` = 16, `p ≥ 0.05` = 21), name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = RES) +
  scale_y_log10(name = "PWH vs PWoH odds ratio\n(donor-level)", limits = c(1, NA)) +
  xlab("Clustering resolution") + theme_astro_pub()
H <- B %>% group_by(subtype, setting) %>% summarise(jaccard = median(jaccard), n = n[1], .groups = "drop")
ord <- H %>% dplyr::filter(setting == "published") %>% arrange(jaccard)
H$subtype <- factor(H$subtype, levels = ord$subtype)
H$setting <- factor(H$setting, levels = names(SHORT))
p3 <- ggplot(H, aes(setting, subtype, fill = jaccard)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", jaccard), colour = jaccard > 0.55), size = 2.1, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08306B", limits = c(0, 1), breaks = c(0, 0.5, 1),
                      name = "Median Jaccard", guide = guide_colourbar(barwidth = 5, barheight = 0.5,
                                                                      title.position = "top")) +
  scale_x_discrete(labels = SHORT, position = "top") + labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") + theme_astro_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
        plot.margin = margin(5, 34, 5, 5),
        axis.text.y = element_text(face = ifelse(levels(H$subtype) == TARGET, "bold", "plain")))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is needed for the figure")
library(patchwork)
fig <- ((p1 / p2) | p3) + plot_layout(widths = c(1.45, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.box = "horizontal", plot.tag = element_text(face = "bold"))
fp <- file.path(OUT, "ask3_robustness_Reactive_PWH.png")
ggsave(fp, fig, width = 7.2, height = 5.6, dpi = 600, bg = "white")
ggsave(sub("\\.png$", ".pdf", fp), fig, width = 7.2, height = 5.6, device = grDevices::cairo_pdf, bg = "white")

## caption numbers come from the tables
pr <- P[P$setting %in% c("published", "pcs20", "pcs50", "k10", "k50"), ]
hd <- P[P$setting == "harmony_donor", ]
ni <- I[I$setting == "no_integration", ]
cap <- sprintf(paste(
  "The seven astrocyte subtypes were re-clustered from scratch under each integration and parameter setting, at",
  "seven resolutions (dashed line: resolution 1, at which the subtypes were defined as merges of marker-defined",
  "clusters). (a) Overlap (Jaccard) between the Reactive-PWH nuclei and the union of new clusters in which",
  "Reactive-PWH nuclei are the majority. (b) Odds ratio, PWH vs PWoH, for that recovered population from the",
  "donor-level composition model (quasibinomial, one observation per donor, adjusted for age, sex, race and PMI).",
  "Across the published setting and the PC and neighbourhood variants, Jaccard was %.2f-%.2f and the odds ratio",
  "%.1f-%.1f (p %.3f-%.3f). With Harmony run by donor rather than by library, recovery was partial (Jaccard",
  "%.2f-%.2f) and the recovered population remained PWH-enriched (odds ratio %.1f-%.1f, p <= %.3f). Without",
  "integration, clusters followed the library (%.0f-%.0f%% of a cluster's nuclei from one donor), so no activation",
  "state was recovered and the odds ratio is not shown. (c) The same recovery measure for every subtype (median over",
  "resolutions); lineage-mixed subtypes (Micro-Astro, Synaptic, Mural-Astro) survive without integration, activation",
  "states (Signaling, Reactive-PWH) do not. Tables: robustness_metrics.csv, subtype_benchmark.csv,",
  "integration_probe.csv."),
  min(pr$jaccard), max(pr$jaccard), min(pr$odds_ratio), max(pr$odds_ratio), min(pr$p), max(pr$p),
  min(hd$jaccard), max(hd$jaccard), min(hd$odds_ratio), max(hd$odds_ratio), max(hd$p),
  100 * min(ni$donor_top_share_weighted), 100 * max(ni$donor_top_share_weighted))
writeLines(strwrap(cap, 110), sub("\\.png$", "_caption.txt", fp))
message("[rob] figure -> ", fp, " (", el(), " min)")
message("[rob] done")
