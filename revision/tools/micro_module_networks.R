# ==============================================================================
# GO networks for the three microglia co-expression modules (Fig 5).
#
# Modules: fixed gene sets from Ruiqi's module analysis (selected_gene_modules.csv):
#   hiv_down_module = Module 1 (80 genes), hiv_up_module1 = Module 2 (48),
#   hiv_up_module2 = Module 3 (85), the naming used in the module figure.
#
# Method: as the published microglia network
#   (03-DE_analysis/09-microglia_cross_region_heatmap_pathways.R):
#   enrichGO BP, BH, p and q < 0.05, default universe (as the published module
#   GO), cnetplot of the top 5 terms, size_item 3 / size_category 1.5.
#   Gene colour = cross-region microglia log2FC from the cross-region table over
#   all tested genes (edgeR by default; SCORCH_XR_DE=donor_aware uses the
#   donor-aware model); a module gene not tested in microglia is drawn
#   uncoloured. Up modules use the published white-to-red ramp, the down module
#   the mirrored blue-to-white.
#   The module gene lists do not change with the reprocessing - only the
#   fold-change colours do, and the GO results only through annotation-database
#   versions.
#
# Input:  OUT_ROOT/microglia_module_analysis/selected_gene_modules.csv; the
#         cross-region DE table (xr_de_file() in revision/config/paths.R)
# Output: OUT_ROOT/microglia_module_analysis/networks[_donor_aware]/
#           module<k>_go_network.png          gene labels, no legend (as published)
#           module<k>_go_network_withlabel.png  term + gene labels, with legend
#           module<k>_go_network.pdf
#           module<k>_GO_BP.csv               full enrichment table
#
# Usage: revision/docker/run.sh revision/tools/micro_module_networks.R [--show 5]
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(clusterProfiler); library(enrichplot)
  library(org.Hs.eg.db); library(AnnotationDbi)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
SHOW <- as.integer(getarg("--show", "5"))
## verbatim (md5-checked) copy of the module analysis's gene-module table
MODULE_CSV <- file.path(OUT_ROOT, "microglia_module_analysis", "selected_gene_modules.csv")
## SCORCH_XR_DE=donor_aware (revision/config/paths.R) colours the nodes by the
## donor-aware cross-region logFC and writes to networks_donor_aware/
OUT <- file.path(OUT_ROOT, "microglia_module_analysis", paste0("networks", XR_OUT_SFX))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

mods <- utils::read.csv(MODULE_CSV, stringsAsFactors = FALSE)
MOD_NAME <- c(hiv_down_module = "Module 1", hiv_up_module1 = "Module 2", hiv_up_module2 = "Module 3")
stopifnot(all(unique(mods$module) %in% names(MOD_NAME)))

## cross-region microglia log2FC, every tested gene (not only DE genes: module
## membership is defined by co-expression, so most module genes are not DE)
de <- readRDS(xr_de_file()) %>%
  ungroup() %>% dplyr::filter(cluster_id == "Micro")
fc_all <- setNames(de$logFC, de$gene)

for (m in names(MOD_NAME)) {
  k <- sub("Module ", "", MOD_NAME[[m]])
  genes <- unique(mods$gene[mods$module == m])
  ent <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL",
                                                columns = "ENTREZID"))
  ent <- ent[!is.na(ent$ENTREZID), ]
  ego <- enrichGO(gene = unique(ent$ENTREZID), OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  res <- as.data.frame(ego)
  utils::write.csv(res, file.path(OUT, sprintf("module%s_GO_BP.csv", k)), row.names = FALSE)
  in_de <- sum(genes %in% names(fc_all))
  message(sprintf("[modnet] %s (%s): %d genes, %d with Entrez IDs, %d tested in microglia DE; %d GO BP terms at FDR < 0.05",
                  MOD_NAME[[m]], m, length(genes), nrow(ent), in_de, nrow(res)))
  if (!nrow(res)) { message("[modnet]   no significant terms - no network drawn"); next }
  message("[modnet]   top terms: ", paste(utils::head(res$Description, SHOW), collapse = " | "))

  fc <- sort(fc_all[intersect(genes, names(fc_all))], decreasing = TRUE)
  down <- grepl("down", m)
  ramp <- if (down) c("#2166AC", "white") else c("white", "red")

  draw <- function(label_mode) {
    p <- if (label_mode == "gene")
      cnetplot(ego, foldChange = fc, size_item = 3, size_category = 1.5,
               node_label = "gene", showCategory = SHOW)
    else
      cnetplot(ego, foldChange = fc, size_item = 3, size_category = 1.5, showCategory = SHOW)
    p + scale_color_gradientn(colours = ramp, name = "log2 FC")
  }
  ## Canvas sized to the network. The published microglia network has ~90 nodes
  ## on an 8.3-in canvas; a module has 15-25, so on the same canvas the nodes
  ## drift far apart. Node and label sizes are absolute, so a smaller canvas
  ## pulls them together without shrinking them.
  n_nodes <- length(unique(unlist(strsplit(utils::head(res$geneID, SHOW), "/")))) + min(SHOW, nrow(res))
  side <- max(4, min(8.33, 1.5 + 0.16 * n_nodes))
  message(sprintf("[modnet]   %d nodes -> %.1f in canvas", n_nodes, side))
  p1 <- draw("gene") + theme(legend.position = "none")
  f1 <- file.path(OUT, sprintf("module%s_go_network.png", k))
  grDevices::png(f1, width = side, height = side, units = "in", res = 600); print(p1); grDevices::dev.off()
  grDevices::cairo_pdf(sub("\\.png$", ".pdf", f1), width = side, height = side); print(p1); grDevices::dev.off()
  p2 <- draw("all")
  grDevices::png(file.path(OUT, sprintf("module%s_go_network_withlabel.png", k)),
                 width = side + 1.2, height = side, units = "in", res = 600); print(p2); grDevices::dev.off()
  message("[modnet]   wrote ", basename(f1))
}
message("[modnet] done -> ", OUT)
