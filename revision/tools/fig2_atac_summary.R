# ==============================================================================
# Figure 2 (ATAC overview) summary panels B, C, E and the between-region
# correlation supplement, from the saved DE result tables only (no ATAC object
# is loaded).
#
#   B  CHROMATIN CHANGES ACROSS CELL TYPES
#      DA peaks per cell type x region, "more accessible in PWH" and "less
#      accessible in PWH" side by side. Same filters as the published count
#      strips (04-ATAC_DE/plotting/01-sum_statistics_atac.R): detected in > 5% of
#      nuclei in either group, |logFC| > 0.25, p_adj.glb < 0.05. Replaces the six
#      per-region strips with one panel. Blank = cell type not tested in that
#      region; 0 = tested, nothing significant.
#
#   C  RNA vs CHROMATIN, PER REGION
#      The published supplementary corr_logfc_heatmap (04-ATAC_DE/cross_region/03-corr_gene_act.R) correlated
#      gene-activity logFC BETWEEN regions. This correlates, WITHIN each region
#      and cell type, the RNA logFC (primary uncorrected DE) with the gene-activity
#      logFC of the same gene (ATAC reads over gene body + promoter), same tile +
#      number design. Two gene sets:
#        de   genes that are RNA DE (p_adj.glb < 0.05, |logFC| > 0.25)   - main panel
#        all  every gene tested in both, detected in > 5% of nuclei      - supplement
#      Cells need >= 20 genes to get a value. Sign agreement among the DE genes
#      is in the CSV.
#
#   S  RNA and gene-activity logFC correlation BETWEEN region pairs, per cell type
#      (the published supplementary heatmap, redrawn on the rerun results).
#
#   E  RNA vs gene-activity logFC scatter of genes significant in both modalities,
#      per region and all regions combined.
#
# INPUTS   OUT_ROOT/03-DE_analysis/<region>/results/<region>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#          OUT_ROOT/04-ATAC_DE/<region>/{data_ATAC,data_GeneACT}/<region>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
# OUTPUTS  OUT_ROOT/04-ATAC_DE/figures/fig2/
#   fig2B_DA_landscape.png/.pdf/_caption.txt + fig2B_DA_counts.csv
#   fig2C_rna_vs_geneactivity_{de,all}.png/.pdf/_caption.txt + fig2C_rna_vs_geneactivity.csv
#   figS_region_corr_{geneactivity,rna}.png/.pdf/_caption.txt/.csv
#   fig2E_rna_vs_geneactivity_scatter[_combined].png/.pdf/_caption.txt + fig2E_rna_vs_geneactivity_scatter.csv
#
# Usage: revision/docker/run.sh revision/tools/fig2_atac_summary.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

OUT <- file.path(OUT_ROOT, "04-ATAC_DE", "figures", "fig2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
REGIONS <- c("PFC", "INS", "VST")
REGION_COLS <- c(PFC = "#E69F00", INS = "#56B4E9", VST = "#CC79A7")
MIN_GENES <- 20

## cell-type classes and order, as the Fig3 pathway panels (MSNs with inhibitory)
MSN <- c("D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT",
         "D2-Striosome","D1/D2-Hybrid","MSNs")
CLASS <- list(
  Excitatory = c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT"),
  Inhibitory = c("Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip", MSN),
  `Glia & vascular` = c("Astro","Oligo","OPC","Micro","Macrophages","Immune","Other immune",
                        "Endo","PC","SMC","VLMC"))
class_of <- function(ct) { for (k in names(CLASS)) if (ct %in% CLASS[[k]]) return(k); "Other" }
ct_levels <- unlist(CLASS, use.names = FALSE)
relab <- function(x) ifelse(x == "Immune", "Other immune", x)

rd <- function(f) { x <- as.data.frame(readRDS(f)); x$cluster_id <- as.character(x$cluster_id); x }
tab <- function(stage, R, kind) {
  f <- if (stage == "rna") file.path(OUT_ROOT, "03-DE_analysis", R, "results",
                                     paste0(R, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))
       else file.path(OUT_ROOT, "04-ATAC_DE", R, kind, paste0(R, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))
  if (!file.exists(f)) stop("missing ", f)
  rd(f)
}
detected <- function(x) (x$avg_CTR_frq > 0.05) | (x$avg_HIV_frq > 0.05)

## ---- B: DA landscape ------------------------------------------------------------
b <- bind_rows(lapply(REGIONS, function(R) {
  x <- tab("atac", R, "data_ATAC")
  tested <- unique(x$cluster_id)
  s <- x[detected(x) & abs(x$logFC) > 0.25 & x$p_adj.glb < 0.05, ]
  expand.grid(region = R, cell_type = tested, direction = c("up", "down"), stringsAsFactors = FALSE) %>%
    left_join(s %>% mutate(direction = ifelse(logFC > 0, "up", "down")) %>%
                count(cell_type = cluster_id, direction, name = "n"), by = c("cell_type", "direction")) %>%
    mutate(n = ifelse(is.na(n), 0L, n))
}))
b$cell_type <- relab(b$cell_type)
b$class <- vapply(b$cell_type, class_of, character(1))
utils::write.csv(b, file.path(OUT, "fig2B_DA_counts.csv"), row.names = FALSE)
message("[fig2] B: DA counts - ", paste(capture.output(print(b %>% group_by(region, direction) %>%
        summarise(n = sum(n), .groups = "drop") %>% as.data.frame(), row.names = FALSE)), collapse = " | "))

lv <- intersect(ct_levels, unique(b$cell_type))
b$cell_type <- factor(b$cell_type, levels = rev(lv))
b$region <- factor(b$region, levels = REGIONS)
b$class <- factor(b$class, levels = names(CLASS))
panel_b <- function(dir, hi, title) {
  d <- b[b$direction == dir, ]
  ggplot(d, aes(x = region, y = cell_type, fill = log10(n + 1))) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = n, colour = log10(n + 1) > 0.6 * log10(max(b$n) + 1)), size = 2.3, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
    scale_fill_gradient(low = "grey96", high = hi, limits = c(0, log10(max(b$n) + 1)), guide = "none") +
    facet_grid(class ~ ., scales = "free_y", space = "free_y") +
    scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_astro_pub(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(face = "bold", size = 7.5, colour = REGION_COLS[REGIONS]),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 7),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          panel.spacing = unit(0.3, "lines"))
}
pB <- panel_b("up", "#B2182B", "More accessible in PWH") +
      theme(strip.text.y = element_blank()) +
      panel_b("down", "#2166AC", "Less accessible in PWH")
fB <- file.path(OUT, "fig2B_DA_landscape.png")
save_pub(pB, fB, width = 4.2, height = 0.5 + 0.15 * length(lv))
astro_caption(fB, paste(
  "Differentially accessible peaks per cell type and region (pseudobulk edgeR, HIV + age + sex + race + PMI;",
  "|log2FC| > 0.25, FDR < 0.05, peak detected in > 5% of nuclei). Left: more accessible in PWH;",
  "right: less accessible. Blank: cell type not tested in that region (too few nuclei / rank-deficient design).",
  "Tile shade is log-scaled."))

## ---- C: RNA vs gene activity, per region ------------------------------------------
cc <- bind_rows(lapply(REGIONS, function(R) {
  rna <- tab("rna", R); ga <- tab("atac", R, "data_GeneACT")
  rna <- rna[detected(rna), c("gene", "cluster_id", "logFC", "p_adj.glb")]
  ga  <- ga[detected(ga), c("gene", "cluster_id", "logFC", "p_adj.glb")]
  j <- inner_join(rna, ga, by = c("gene", "cluster_id"), suffix = c("_rna", "_ga"))
  j %>% group_by(cluster_id) %>% group_modify(function(g, k) {
    de <- g[g$p_adj.glb_rna < 0.05 & abs(g$logFC_rna) > 0.25, ]
    r_or_na <- function(x, y) if (length(x) >= MIN_GENES) suppressWarnings(cor(x, y)) else NA_real_
    data.frame(n_all = nrow(g), r_all = r_or_na(g$logFC_rna, g$logFC_ga),
               n_de = nrow(de), r_de = r_or_na(de$logFC_rna, de$logFC_ga),
               sign_agree_de = if (nrow(de)) mean(sign(de$logFC_rna) == sign(de$logFC_ga)) else NA_real_)
  }) %>% ungroup() %>% mutate(region = R)
}))
cc <- cc %>% rename(cell_type = cluster_id) %>% mutate(cell_type = relab(cell_type))
cc$class <- vapply(cc$cell_type, class_of, character(1))
utils::write.csv(cc, file.path(OUT, "fig2C_rna_vs_geneactivity.csv"), row.names = FALSE)
message("[fig2] C: ", nrow(cc), " region x cell type pairs; with >= ", MIN_GENES, " RNA DE genes: ",
        sum(!is.na(cc$r_de)))

draw_c <- function(col, nvar, stub, what) {
  d <- cc[!is.na(cc[[col]]), ]
  lvc <- intersect(ct_levels, unique(d$cell_type))
  d$cell_type <- factor(d$cell_type, levels = rev(lvc))
  d$region <- factor(d$region, levels = REGIONS)
  d$class <- factor(d$class, levels = names(CLASS))
  d$val <- d[[col]]; d$n <- d[[nvar]]
  p <- ggplot(d, aes(x = region, y = cell_type, fill = val)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", val)), size = 2.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-1, 1), name = "Pearson r") +
    facet_grid(class ~ ., scales = "free_y", space = "free_y") +
    scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_astro_pub(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(face = "bold", size = 7.5, colour = REGION_COLS[REGIONS]),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 7),
          legend.key.height = unit(0.45, "cm"), legend.key.width = unit(0.25, "cm"),
          panel.spacing = unit(0.3, "lines"))
  f <- file.path(OUT, paste0("fig2C_rna_vs_geneactivity_", stub, ".png"))
  save_pub(p, f, width = 3.2, height = 0.5 + 0.15 * length(lvc))
  astro_caption(f, paste0(
    "Agreement between RNA and chromatin changes in PWH, per region and cell type: Pearson correlation of ",
    "the RNA log2 fold change (pseudobulk edgeR, uncorrected counts) and the gene-activity log2 fold change ",
    "(ATAC reads over gene body and promoter, pseudobulk edgeR) for ", what, ". Cells with fewer than ",
    MIN_GENES, " genes are blank. Gene counts and sign agreement: fig2C_rna_vs_geneactivity.csv."))
}
draw_c("r_de", "n_de", "de", "genes differentially expressed in RNA (FDR < 0.05, |log2FC| > 0.25)")
draw_c("r_all", "n_all", "all", "all genes tested in both, detected in > 5% of nuclei")
## ---- S: between-region correlation, as the published supplementary heatmap ----
## Published: 04-ATAC_DE/cross_region/03-corr_gene_act.R - per cell type, Pearson
## r of the gene-activity logFC between two regions over every gene detected
## (> 5% of nuclei) in both. Drawn for gene activity (the published supplement)
## and for RNA (the RNA version is recomputed here so it uses the VST DE without
## donor 1040). Same tile style as panel C.
PAIRS <- list(c("INS", "PFC"), c("INS", "VST"), c("PFC", "VST"))
pair_corr <- function(stage) {
  tb <- lapply(setNames(REGIONS, REGIONS), function(R) {
    x <- if (stage == "rna") tab("rna", R) else tab("atac", R, "data_GeneACT")
    x[detected(x), c("gene", "cluster_id", "logFC")]
  })
  bind_rows(lapply(PAIRS, function(p) {
    j <- inner_join(tb[[p[1]]], tb[[p[2]]], by = c("gene", "cluster_id"))
    j %>% group_by(cluster_id) %>%
      summarise(n = dplyr::n(), r = if (dplyr::n() >= MIN_GENES) cor(logFC.x, logFC.y) else NA_real_,
                .groups = "drop") %>%
      mutate(pair = paste(p, collapse = "-"))
  })) %>% rename(cell_type = cluster_id) %>% mutate(cell_type = relab(cell_type))
}
for (stage in c("geneactivity", "rna")) {
  pc <- pair_corr(stage)
  pc$class <- vapply(pc$cell_type, class_of, character(1))
  utils::write.csv(pc, file.path(OUT, paste0("figS_region_corr_", stage, ".csv")), row.names = FALSE)
  d <- pc[!is.na(pc$r), ]
  lvp <- intersect(ct_levels, unique(d$cell_type))
  d$cell_type <- factor(d$cell_type, levels = rev(lvp))
  d$pair  <- factor(d$pair, levels = vapply(PAIRS, paste, character(1), collapse = "-"))
  d$class <- factor(d$class, levels = names(CLASS))
  p <- ggplot(d, aes(x = pair, y = cell_type, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-1, 1), name = "Pearson r") +
    facet_grid(class ~ ., scales = "free_y", space = "free_y") +
    ## pair names on two lines: "INS-PFC" etc. collided on one line at 3.4 in
    scale_x_discrete(position = "top", expand = c(0, 0), labels = function(x) sub("-", " vs\n", x)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_astro_pub(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(face = "bold", size = 7, lineheight = 0.9),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 7),
          legend.key.height = unit(0.45, "cm"), legend.key.width = unit(0.25, "cm"),
          panel.spacing = unit(0.3, "lines"))
  f <- file.path(OUT, paste0("figS_region_corr_", stage, ".png"))
  save_pub(p, f, width = 3.4, height = 0.5 + 0.15 * length(lvp))
  astro_caption(f, paste0(
    "Consistency of ", if (stage == "rna") "RNA" else "gene-activity (ATAC)", " changes in PWH between ",
    "regions: per cell type, Pearson correlation of the ", if (stage == "rna") "RNA" else "gene-activity",
    " log2 fold change between two regions over all genes detected (> 5% of nuclei) in both. Cell types ",
    "tested in only one region of a pair, or with fewer than ", MIN_GENES, " shared genes, are blank."))
  message("[fig2] region-pair correlation (", stage, "): ", nrow(d), " cell type x pair values")
  print(as.data.frame(pc %>% filter(cell_type %in% c("Astro", "Micro", "OPC", "Oligo")) %>%
          mutate(r = round(r, 2)) %>% dplyr::select(cell_type, pair, r) %>%
          tidyr::pivot_wider(names_from = pair, values_from = r)), row.names = FALSE)
}

## ---- E: RNA vs gene-activity logFC scatter, ONE PANEL PER REGION ------------------
## Published: 05-RNA_ATAC_DE/plotting/01-sum_statistics.R -> geom_logfc_gene_act_new.png:
## genes significant in BOTH the RNA DE and the gene-activity DE of the same cell
## type (each: > 5% detected, |logFC| > 0.25, p_adj.glb < 0.05), RNA logFC on x,
## gene-activity logFC on y, coloured by cell type, all three regions on one plot
## with region as shape. Here: one panel per region, then the combined plot below.
## Every point is kept; N_LABEL sets how many are labelled, ranked by distance from
## the origin (each axis scaled by its SD within the region); Inf = all. No forced
## labels. Colours: the 37-level cell-type palette of the published figure.
CT_COLS <- setNames(c(
  "indianred1","tomato1","sienna1","sienna3","maroon","#6A3D9A","#FB9A99","#CAB2D6","orchid1",
  "purple1","lightslateblue","orchid3","pink3","plum1","tomato1","tomato3","sienna2","burlywood3",
  "green4","palegreen2","darkseagreen","dodgerblue2","skyblue3","darkturquoise","skyblue1","steelblue4","deepskyblue",
  "#FF7F00","#FDBF6F","khaki2","gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4","yellow3"),
  c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
    "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
    "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
    "Astro","Oligo","OPC","Micro","Macrophages","Other immune","PC","Endo","SMC","VLMC"))
## All genes labelled, as in the published figure; no forced labels.
N_LABEL <- Inf; ALWAYS <- character(0)
sig <- function(x) x[detected(x) & abs(x$logFC) > 0.25 & x$p_adj.glb < 0.05, c("gene", "cluster_id", "logFC")]
sc <- bind_rows(lapply(REGIONS, function(R)
  inner_join(sig(tab("rna", R)), sig(tab("atac", R, "data_GeneACT")),
             by = c("gene", "cluster_id"), suffix = c("_rna", "_ga")) %>% mutate(region = R)))
sc$cell_type <- relab(sc$cluster_id)
sc <- sc %>% group_by(region) %>%
  mutate(dist = sqrt((logFC_rna / sd(logFC_rna))^2 + (logFC_ga / sd(logFC_ga))^2),
         lab = ifelse(rank(-dist, ties.method = "first") <= N_LABEL | gene %in% ALWAYS, gene, NA_character_),
         concordant = sign(logFC_rna) == sign(logFC_ga)) %>% ungroup()
utils::write.csv(sc %>% dplyr::select(region, cell_type, gene, logFC_rna, logFC_ga, concordant),
                 file.path(OUT, "fig2E_rna_vs_geneactivity_scatter.csv"), row.names = FALSE)
st <- sc %>% group_by(region) %>%
  summarise(n = dplyr::n(), pct = round(100 * mean(concordant)), r = cor(logFC_rna, logFC_ga), .groups = "drop")
message("[fig2] E: genes significant in both - ",
        paste(sprintf("%s n=%d (%d%% same direction, r=%.2f)", st$region, st$n, st$pct, st$r), collapse = "; "))
strip_lab <- setNames(sprintf("%s\nn = %s genes, %d%% same direction", st$region,
                              formatC(st$n, big.mark = ","), st$pct), st$region)
sc$region <- factor(sc$region, levels = REGIONS)
present <- intersect(names(CT_COLS), unique(sc$cell_type))
sc$cell_type <- factor(sc$cell_type, levels = present)
pE <- ggplot(sc, aes(x = logFC_rna, y = logFC_ga, colour = cell_type)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey40") +
  geom_point(size = 0.9, alpha = 0.8, stroke = 0) +
  ggrepel::geom_text_repel(aes(label = lab), size = 1.9, force = 2, colour = "black", na.rm = TRUE,
                           max.overlaps = Inf, box.padding = 0.25, min.segment.length = 0.1,
                           segment.size = 0.2, segment.colour = "grey50", seed = 1) +
  scale_colour_manual(values = CT_COLS[present], name = NULL) +
  facet_wrap(~ region, nrow = 1, scales = "free", labeller = as_labeller(strip_lab)) +
  labs(x = "log2 FC, gene expression (RNA)", y = "log2 FC, gene activity (ATAC)") +
  theme_astro_pub(base_size = 7) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.25, "cm"), strip.text = element_text(size = 7.5, lineheight = 0.9)) +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 2, alpha = 1)))
fE <- file.path(OUT, "fig2E_rna_vs_geneactivity_scatter.png")
save_pub(pE, fE, width = 7.2, height = 4.4)
astro_caption(fE, paste0(
  "Effect sizes of genes changed in PWH in both gene expression (x, RNA pseudobulk edgeR) and gene activity ",
  "(y, ATAC reads over gene body and promoter, pseudobulk edgeR) in the same cell type (each FDR < 0.05, ",
  "|log2FC| > 0.25, detected in > 5% of nuclei); one panel per region, coloured by cell type; every gene ",
  "is labelled. Strip: number of genes and the percentage changing in the same ",
  "direction in both modalities."))

## E, combined: all regions on one plot as published (VST alone has very few
## genes): region = shape, cell type = colour. Every gene labelled
## (N_LABEL = Inf); no forced labels.
sc2 <- sc %>% mutate(dist_all = sqrt((logFC_rna / sd(logFC_rna))^2 + (logFC_ga / sd(logFC_ga))^2),
                     lab_all = ifelse(rank(-dist_all, ties.method = "first") <= N_LABEL | gene %in% ALWAYS,
                                      gene, NA_character_))
REGION_SHAPES <- c(PFC = 16, INS = 17, VST = 15)
pE2 <- ggplot(sc2, aes(x = logFC_rna, y = logFC_ga, colour = cell_type, shape = region)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey40") +
  geom_point(size = 1.6, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = lab_all), size = 2.0, force = 2, colour = "black", na.rm = TRUE,
                           max.overlaps = Inf, box.padding = 0.3, min.segment.length = 0.1,
                           segment.size = 0.2, segment.colour = "grey50", seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = CT_COLS[present], name = "Cell type") +
  scale_shape_manual(values = REGION_SHAPES,
                     labels = sprintf("%s (n = %d)", st$region, st$n)[match(REGIONS, st$region)],
                     breaks = REGIONS, name = "Region") +
  labs(x = "log2 FC, gene expression (RNA)", y = "log2 FC, gene activity (ATAC)") +
  theme_astro_pub(base_size = 7) +
  theme(legend.position = "right", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm")) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 2)),
         colour = guide_legend(order = 2, ncol = 1, override.aes = list(size = 2, shape = 16)))
fE2 <- file.path(OUT, "fig2E_rna_vs_geneactivity_scatter_combined.png")
save_pub(pE2, fE2, width = 7.2, height = 5.6)
astro_caption(fE2, paste0(
  "Effect sizes of genes changed in PWH in both gene expression (x) and gene activity (y) in the same cell ",
  "type and region (each FDR < 0.05, |log2FC| > 0.25, detected in > 5% of nuclei); colour, cell type; shape, ",
  "region (PFC ", st$n[st$region == "PFC"], ", INS ", st$n[st$region == "INS"], ", VST ",
  st$n[st$region == "VST"], " genes). ", round(100 * mean(sc$concordant)), "% of genes change in the same ",
  "direction in both modalities. Every gene is labelled."))
message("[fig2] E combined: ", nrow(sc), " genes, ", round(100 * mean(sc$concordant)), "% same direction")

print(as.data.frame(cc %>% filter(cell_type %in% c("Astro", "Micro", "Oligo", "OPC")) %>%
        mutate(across(c(r_all, r_de, sign_agree_de), ~ round(.x, 2))) %>% arrange(cell_type, region)),
      row.names = FALSE)
message("[fig2] done - ", OUT)
