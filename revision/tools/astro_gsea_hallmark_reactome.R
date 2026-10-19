# ==============================================================================
# GSEA of the Reactive-PWH astrocyte donor-level DE on MSigDB Hallmark and
# Reactome gene sets, and the Fig1E heatmap from the Hallmark leading edge.
#
# Method
#   Genes are ranked by sign(logFC) x -log10(PValue) (p floored at one tenth of the
#   smallest non-zero p; genes with negative QL F dropped) with no significance
#   threshold, and tested with clusterProfiler::GSEA (Hallmark) and
#   ReactomePA::gsePathway (Reactome), BH-adjusted p < 0.05. The same ranking is
#   used for GO BP; there, ciliary and developmental sets (motile cilium assembly,
#   spermatid development, neural tube patterning) fill the down direction of the
#   vs-Reactive contrasts, annotation artefacts in postmortem astrocytes. Hallmark
#   is 50 curated, non-redundant sets with no parent/child pairs; Reactome is
#   pathway-level and better defined than GO's process vocabulary.
#   Reactome needs Entrez ids, so symbols are mapped before ranking and results are
#   mapped back with setReadable, so leading-edge genes come out as symbols.
#   Three DE tables are tested: Reactive-PWH vs all other astrocytes (published
#   contrast), vs Reactive astrocytes (24 donors), and vs Reactive in PWH donors
#   only (20 donors). The NES figure shows the two vs-Reactive contrasts; vs_all is
#   written to csv for comparison with the published contrast. Overlapping sets
#   (leading-edge Jaccard >= 0.5) are pruned for the bar figure; top 8 per direction,
#   set size >= 15.
#
# Inputs   OUT_ROOT/06-astrocyte_analysis/astro_Reactive_PWH*_donorDE_no1040vst.csv
#          (and pwh_only/), OUT_ROOT/reference/h.all.v2024.1.Hs.symbols.gmt
# Outputs  OUT_ROOT/06-astrocyte_analysis/figures/Fig1E_gsea/
#            gsea_<table>_<db>.csv              full results per table and database
#            gsea_hallmark_reactome.png/.pdf    NES bars, up and down, per table
#            gsea_hallmark_reactome_plotted.csv
#          OUT_ROOT/06-astrocyte_analysis/figures/Fig1E_heatmaps/
#            Fig1E_hallmark_heatmap_no1040vst.png/.pdf  Fig1E layout from the Hallmark
#                                     leading edge of the vs-Reactive (24 donors) contrast
#
# Usage: revision/docker/run.sh revision/tools/astro_gsea_hallmark_reactome.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
TOOLS <- file.path(SCORCH_REPO, "revision/tools")
for (f in c("astro_pub_theme", "astro_cluster_features")) source(file.path(TOOLS, paste0(f, ".R")))
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
  library(ggplot2); library(dplyr); library(stringr); library(patchwork)
})

D    <- file.path(OUT_ROOT, "06-astrocyte_analysis")
OUT  <- file.path(D, "figures", "Fig1E_gsea")
HM   <- file.path(D, "figures", "Fig1E_heatmaps")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(HM, recursive = TRUE, showWarnings = FALSE)
GMT  <- file.path(OUT_ROOT, "reference", "h.all.v2024.1.Hs.symbols.gmt")
SFX  <- "_no1040vst"
N_DIR <- 8; MIN_SET <- 15
MARKERS <- c("SERPINA3", "OSMR", "CHI3L1")

TABLES <- list(
  vs_all          = file.path(D, "astro_Reactive_PWH_donorDE_no1040vst.csv"),
  vs_reactive_24  = file.path(D, "astro_Reactive_PWH_vs_Reactive_donorDE_no1040vst.csv"),
  vs_reactive_pwh = file.path(D, "pwh_only", "astro_Reactive_PWH_vs_Reactive_donorDE_no1040vst.csv"))
LABS <- c(vs_all = "vs all other astrocytes\n(published contrast)",
          vs_reactive_24 = "vs Reactive astrocytes\n(24 donors)",
          vs_reactive_pwh = "vs Reactive astrocytes\n(20 PWH donors)")

## HALLMARK_TNFA_SIGNALING_VIA_NFKB -> "TNFα signaling via NF-κB"
hallmark_label <- function(x) {
  x <- tolower(gsub("_", " ", sub("^HALLMARK_", "", x)))
  fix <- c(tnfa = "TNFα", nfkb = "NF-κB", il6 = "IL6", il2 = "IL2", jak = "JAK", stat3 = "STAT3",
           stat5 = "STAT5", kras = "KRAS", mtorc1 = "mTORC1", e2f = "E2F", myc = "MYC", g2m = "G2M",
           tgf = "TGF", beta = "β", uv = "UV", dn = "down", pi3k = "PI3K", akt = "AKT", mtor = "mTOR",
           dna = "DNA", wnt = "WNT", notch = "NOTCH", p53 = "p53", ifn = "IFN")
  vapply(strsplit(x, " "), function(w) {
    w <- ifelse(w %in% names(fix), fix[w], w)
    s <- paste(w, collapse = " "); paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  }, character(1))
}

## the ranking, identical to the one used for the GO BP GSEA
rank_vec <- function(de) {
  p <- de$PValue
  floor_p <- min(p[p > 0], na.rm = TRUE) / 10
  v <- sign(de$logFC) * -log10(pmax(p, floor_p))
  names(v) <- de$gene
  sort(v[is.finite(v) & !duplicated(names(v))], decreasing = TRUE)
}

h_sets <- clusterProfiler::read.gmt(GMT)
message("[hr] Hallmark: ", length(unique(h_sets$term)), " sets from ", basename(GMT))

res_all <- list()
for (tn in names(TABLES)) {
  f <- TABLES[[tn]]
  if (!file.exists(f)) { message("[hr] missing: ", f); next }
  de <- utils::read.csv(f, stringsAsFactors = FALSE)
  negF <- sum(de$F < 0, na.rm = TRUE)
  if (negF > 0) message("[hr] ", tn, ": dropping ", negF, " genes with negative F")
  de <- de[is.na(de$F) | de$F >= 0, ]
  v  <- rank_vec(de)
  message("\n[hr] === ", tn, ": ", length(v), " ranked genes ===")

  ## ---- Hallmark, on symbols ------------------------------------------------------------
  h <- try(clusterProfiler::GSEA(geneList = v, TERM2GENE = h_sets, pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH", verbose = FALSE, seed = TRUE), silent = TRUE)
  if (!inherits(h, "try-error") && nrow(as.data.frame(h))) {
    r <- as.data.frame(h)
    r$label <- hallmark_label(r$ID); r$table <- tn; r$db <- "Hallmark"
    utils::write.csv(r, file.path(OUT, sprintf("gsea_%s_hallmark.csv", tn)), row.names = FALSE)
    res_all[[paste(tn, "Hallmark")]] <- r
    message(sprintf("[hr] Hallmark  %3d sets | top up: %s (NES %.2f) | top down: %s (NES %.2f)",
                    nrow(r), r$label[which.max(r$NES)], max(r$NES),
                    r$label[which.min(r$NES)], min(r$NES)))
  } else message("[hr] Hallmark: nothing significant")

  ## ---- Reactome, on entrez, mapped back to symbols --------------------------------------
  ez <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = names(v), keytype = "SYMBOL",
                                               column = "ENTREZID", multiVals = "first"))
  ve <- v[!is.na(ez)]; names(ve) <- ez[!is.na(ez)]
  ve <- sort(ve[!duplicated(names(ve))], decreasing = TRUE)
  rp <- try({
    suppressPackageStartupMessages(library(ReactomePA))
    g <- ReactomePA::gsePathway(geneList = ve, organism = "human", pvalueCutoff = 0.05,
                                pAdjustMethod = "BH", verbose = FALSE, seed = TRUE)
    DOSE::setReadable(g, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }, silent = TRUE)
  if (!inherits(rp, "try-error") && nrow(as.data.frame(rp))) {
    r <- as.data.frame(rp)
    r$label <- r$Description; r$table <- tn; r$db <- "Reactome"
    utils::write.csv(r, file.path(OUT, sprintf("gsea_%s_reactome.csv", tn)), row.names = FALSE)
    res_all[[paste(tn, "Reactome")]] <- r
    message(sprintf("[hr] Reactome  %3d sets | top up: %s (NES %.2f) | top down: %s (NES %.2f)",
                    nrow(r), r$label[which.max(r$NES)], max(r$NES),
                    r$label[which.min(r$NES)], min(r$NES)))
  } else message("[hr] Reactome: nothing significant (or ReactomePA unavailable)")
}

res <- dplyr::bind_rows(res_all)
if (!nrow(res)) stop("no results from either database")

## ---- NES bars, same cleaning as the GO figure -------------------------------------------
prune_overlap <- function(r) {
  r  <- r[order(r$p.adjust), ]
  ce <- strsplit(r$core_enrichment, "/", fixed = TRUE)
  keep <- logical(nrow(r))
  for (i in seq_len(nrow(r)))
    keep[i] <- !any(vapply(which(keep), function(j)
      length(intersect(ce[[i]], ce[[j]])) / length(union(ce[[i]], ce[[j]])), numeric(1)) >= 0.5)
  r[keep, ]
}
## The figure shows the two Reactive-comparator panels only. vs_all is still
## computed and written to csv, to document what changed relative to the
## published contrast.
FIG_TABLES <- c("vs_reactive_24", "vs_reactive_pwh")
pf <- res %>% filter(setSize >= MIN_SET, table %in% FIG_TABLES) %>% group_split(table, db) %>%
  lapply(function(d) {
    d <- prune_overlap(d)
    dplyr::bind_rows(d %>% filter(NES > 0) %>% slice_max(NES, n = N_DIR),
                     d %>% filter(NES < 0) %>% slice_min(NES, n = N_DIR))
  }) %>% dplyr::bind_rows()

if (nrow(pf)) {
  pf <- pf %>%
    mutate(table = factor(table, levels = FIG_TABLES), db = factor(db, c("Hallmark", "Reactome")),
           direction = ifelse(NES > 0, "up in Reactive-PWH", "down in Reactive-PWH"),
           ## zero-width spaces keep a shared term orderable within each panel
           lab = paste0(str_wrap(label, 40),
                        strrep("​", as.integer(table) + 3 * as.integer(db)))) %>%
    arrange(db, table, NES) %>% mutate(lab = factor(lab, levels = unique(lab)))
  rng <- max(abs(pf$NES)) * 1.08
  p <- ggplot(pf, aes(NES, lab, fill = direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    facet_grid(db ~ table, scales = "free_y", space = "free_y",
               labeller = labeller(table = LABS)) +
    scale_fill_manual(values = c(`up in Reactive-PWH` = "#B2182B",
                                 `down in Reactive-PWH` = "#2166AC"), name = NULL) +
    scale_x_continuous(limits = c(-rng, rng)) +
    labs(x = "normalised enrichment score", y = NULL) +
    theme_astro_pub(7) +
    theme(legend.position = "top", legend.margin = margin(b = -4),
          axis.text.y = element_text(size = 6, lineheight = 0.9),
          panel.spacing = unit(0.8, "lines"), strip.text.y = element_text(face = "bold"))
  ggsave(file.path(OUT, "gsea_hallmark_reactome.png"), p, width = 8.4, height = 8, dpi = 600, bg = "white")
  ggsave(file.path(OUT, "gsea_hallmark_reactome.pdf"), p, width = 8.4, height = 8,
         device = grDevices::cairo_pdf, bg = "white")
  utils::write.csv(pf[, c("db", "table", "label", "setSize", "NES", "p.adjust")],
                   file.path(OUT, "gsea_hallmark_reactome_plotted.csv"), row.names = FALSE)
  message("\n[hr] plotted ", nrow(pf), " sets")
}

## ---- the Fig1E heatmap from the Hallmark leading edge -----------------------------------
hm <- res_all[["vs_reactive_24 Hallmark"]]
de <- utils::read.csv(TABLES$vs_reactive_24, stringsAsFactors = FALSE)
if (!is.null(hm)) {
  tr <- hm[hm$NES > 0, ]; tr <- utils::head(tr[order(tr$p.adjust), ], 10)
  in_terms <- unique(unlist(strsplit(tr$core_enrichment, "/", fixed = TRUE)))
  d_in <- de[de$gene %in% in_terms, ]
  show <- utils::head(d_in$gene[order(-d_in$logFC)], 40)
  show <- c(show, setdiff(intersect(MARKERS, in_terms), show))
  long <- do.call(rbind, lapply(seq_len(nrow(tr)), function(i) {
    g <- intersect(strsplit(tr$core_enrichment[i], "/", fixed = TRUE)[[1]], show)
    if (length(g)) data.frame(term = tr$label[i], gene = g, stringsAsFactors = FALSE) else NULL }))
  long$logFC <- de$logFC[match(long$gene, de$gene)]
  first <- tapply(match(long$term, tr$label), long$gene, min)
  lfc   <- tapply(long$logFC, long$gene, max)
  gene_order <- names(first)[order(first, -lfc[names(first)])]
  wrap <- function(x) vapply(x, function(s) paste(strwrap(s, 34), collapse = "\n"), character(1))
  term_lv <- rev(wrap(tr$label))
  long$term <- factor(wrap(long$term), levels = term_lv)
  long$gene <- factor(long$gene, levels = gene_order)
  bars <- data.frame(term = factor(wrap(tr$label), levels = term_lv), mlp = -log10(tr$p.adjust))
  p_hm <- ggplot(long, aes(gene, term, fill = logFC)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_gradient(low = "#FCDCCF", high = "#B2182B",
                        name = "log2 FC\n(Reactive-PWH vs\nReactive)") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) + theme_astro_pub(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.3, face = "italic"),
          axis.text.y = element_text(size = 7.6, lineheight = 0.85),
          panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.3),
          legend.position = "right", legend.key.height = unit(0.45, "cm"),
          legend.key.width = unit(0.25, "cm"), legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6.5))
  p_bar <- ggplot(bars, aes(mlp, term)) + geom_col(fill = "grey55", width = 0.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08)), n.breaks = 3) +
    labs(x = expression(-log[10]~adj.~p), y = NULL) + theme_astro_pub(base_size = 7) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.title.x = element_text(size = 6.5))
  f <- file.path(HM, paste0("Fig1E_hallmark_heatmap", SFX, ".png"))
  save_pub(p_hm + p_bar + plot_layout(widths = c(1, 0.12), guides = "collect"),
           f, width = 7.6, height = 1.1 + 0.3 * nrow(tr))
  astro_caption(f, paste0(
    "MSigDB Hallmark gene sets enriched in Reactive-PWH astrocytes relative to Reactive astrocytes ",
    "(24 donors; gene set enrichment analysis over all ", nrow(de), " tested genes ranked by ",
    "sign(log2 fold change) x -log10 p, no threshold; BH-adjusted p < 0.05). Tiles show leading-edge ",
    "genes for the top ", nrow(tr), " sets with a positive normalised enrichment score, limited to the ",
    "40 with the largest fold change plus SERPINA3, OSMR and CHI3L1 where present; bars, ",
    "-log10 adjusted p. Hallmark sets are curated to be non-overlapping, so no redundancy filter is applied."))
  message("[hr] wrote ", basename(f))
}
message("[hr] done: ", OUT)
