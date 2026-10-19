# ==============================================================================
# Shared helpers for the Figure 4 region-pair logFC scatters: the published
# cell-type palette, the region-pair data and the published label rule, defined
# once so every figure that draws these scatters uses the same definitions.
#
#   get(R)     regional edgeR DE, genes detected in > 5% of nuclei in either group
#   prep(p)    genes tested in the same cell type in both regions of a pair;
#              DE in a region = |logFC| > 0.25 and p_adj.glb < 0.05
#   pub_lab()  the published label rule (|logFC| > 1 and FDR < 0.05 in both
#              regions, or in one region with |logFC| < 0.25 and FDR > 0.05 in the
#              other), top n by |logFC.x| + |logFC.y|, one label per gene
#
# INPUT   OUT_ROOT/03-DE_analysis/<region>/results/<region>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
# Usage:  sourced, e.g. by revision/tools/fig4_cross_region_scatter.R.
#         Requires dplyr and OUT_ROOT (revision/config/paths.R).
# ==============================================================================

PAIRS <- list(c("INS", "PFC"), c("INS", "VST"), c("PFC", "VST"))

## palettes, as 03-DE_analysis/06-compare_logFC_cross_region.R
LEV <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
         "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
         "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
         "Astro","Oligo","OPC","Micro","Macrophages","Other immune","PC","Endo","SMC","VLMC")
cp_col <- setNames(c(
  "indianred1","tomato1","sienna1","sienna3","maroon","#6A3D9A","#FB9A99","#CAB2D6","orchid1",
  "purple1","lightslateblue","orchid3","pink3","plum1","tomato1","tomato3","sienna2","burlywood3",
  "green4","palegreen2","darkseagreen","dodgerblue2","skyblue3","darkturquoise","skyblue1","steelblue4","deepskyblue",
  "#FF7F00","#FDBF6F","khaki2","gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4","yellow3"), LEV)
cp_col2 <- c(Ext = "tomato2", MSNs = "royalblue3", Inh = "dodgerblue", Glial = "gold2")
COARSE <- setNames(c(rep("Ext", 9), rep("MSNs", 9), rep("Inh", 9), rep("Glial", 10)), LEV)

get <- function(R) {
  x <- readRDS(file.path(OUT_ROOT, "03-DE_analysis", R, "results",
                         paste0(R, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS")))
  x %>% ungroup() %>% dplyr::filter(avg_CTR_frq > 0.05 | avg_HIV_frq > 0.05) %>%
    dplyr::select(gene, cluster_id, logFC, p_adj.glb)
}

prep <- function(p) {
  d <- inner_join(REG[[p[1]]], REG[[p[2]]], by = c("gene", "cluster_id"), suffix = c(".x", ".y"))
  d$cluster_id <- as.character(d$cluster_id)
  d$cluster_id[d$cluster_id == "Immune"] <- "Other immune"
  d %>% mutate(
    pair  = paste(p, collapse = " vs "), x_reg = p[1], y_reg = p[2],
    DEx   = abs(logFC.x) > 0.25 & p_adj.glb.x < 0.05,
    DEy   = abs(logFC.y) > 0.25 & p_adj.glb.y < 0.05,
    state = factor(ifelse(DEx & DEy, "both", ifelse(DEx | DEy, "one", "neither")),
                   levels = c("neither", "one", "both")),
    cell_type   = factor(cluster_id, levels = LEV),
    coarse      = factor(COARSE[cluster_id], levels = names(cp_col2)),
    concordant  = sign(logFC.x) == sign(logFC.y),
    size_metric = abs(logFC.x) + abs(logFC.y))
}

## force / per_type: after the top n, each cell type in `force` that has
## no label yet gets its per_type strongest genes under the SAME published rule,
## skipping genes already labelled. A forced type with no rule-eligible gene gets
## nothing and is reported - it is not relabelled under a looser rule.
## genes = c("CD163:Micro", "B2M:L4 IT"): the examples the Results text discusses,
## labelled in every panel where that gene x cell type was tested in both regions,
## whether or not the published rule would have picked it (the text names them).
## updown = TRUE: n labels UP and n DOWN instead of n overall. Direction is
## the sign of logFC.x + logFC.y, i.e. of the larger change for opposite-sign genes.
## A gene is labelled once, at its strongest cell type, so a gene up in one cell
## type and down in another lands in whichever direction that cell type gives.
pub_lab <- function(d, n = N_PUB, force = character(0), per_type = 1, genes = character(0),
                    updown = FALSE) {
  topx <- abs(d$logFC.x) > 1 & d$p_adj.glb.x < 0.05
  topy <- abs(d$logFC.y) > 1 & d$p_adj.glb.y < 0.05
  nonx <- abs(d$logFC.x) < 0.25 & d$p_adj.glb.x > 0.05
  nony <- abs(d$logFC.y) < 0.25 & d$p_adj.glb.y > 0.05
  e <- d[(topx & topy) | (topx & nony) | (topy & nonx), ] %>% arrange(desc(size_metric))
  l <- e %>% distinct(gene, .keep_all = TRUE)
  l <- if (updown) {
    up <- (l$logFC.x + l$logFC.y) > 0
    bind_rows(head(l[up, ], n), head(l[!up, ], n))
  } else head(l, n)
  for (ct in setdiff(force, as.character(l$cell_type))) {
    add <- e %>% dplyr::filter(as.character(cell_type) == ct, !gene %in% l$gene) %>%
      distinct(gene, .keep_all = TRUE) %>% head(per_type)
    if (!nrow(add) && ct %in% as.character(d$cell_type))
      message("[fig4] ", d$pair[1], ": no gene of ", ct, " meets the label rule - not forced")
    l <- bind_rows(l, add)
  }
  for (gc in genes) {
    g <- sub(":.*$", "", gc); ct <- sub("^[^:]*:", "", gc)
    if (any(l$gene == g & as.character(l$cell_type) == ct)) next
    add <- d[d$gene == g & as.character(d$cell_type) == ct, ]
    if (!nrow(add)) { message("[fig4] ", d$pair[1], ": ", gc, " not tested in both regions"); next }
    ## the same gene labelled at another cell type would show twice - keep the named one
    l <- bind_rows(l[l$gene != g, ], add)
  }
  l
}
