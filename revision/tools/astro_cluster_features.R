# ==============================================================================
# Astrocyte subcluster characterisation: marker dotplots, data-derived subtype
# label suggestions, vascular-contamination check, and plotting/enrichment helpers
# (function library, sourced by 06-astrocyte_analysis/astrocyte_analysis.Rmd and
# the astrocyte figure tools).
#
# Method
#   The published astrocyte_analysis.Rmd maps cluster indices to subtype labels by
#   hand (e.g. "7" = "HIV-Reactive"). Those indices come from FindClusters on the
#   published object, and clusters renumber after ambient correction, so labels
#   are derived from the data instead: FindAllMarkers (Wilcoxon; min.pct 0.25,
#   logFC 0.25) per cluster, plus the fraction of each cluster expressing each
#   curated subtype panel (ASTRO_PANELS). A panel is scored by its best-detected
#   gene, scaled to the highest-scoring cluster for that panel, and the top panel
#   is suggested; a margin < 0.15 over the runner-up is marked not confident.
#   Nothing is renamed automatically.
#   Helpers cluster_gene_pct() and pwh_relabel() come from revision/config/paths.R.
#
# Outputs (astro_features / astro_save_plots)  OUT_ROOT/06-astrocyte_analysis/
#   astro_cluster_markers_all.csv, astro_cluster_markers_top.csv,
#   astro_cluster_suggested_labels.csv, astro_cluster_sizes.csv,
#   astro_marker_dotplot*.png, astro_top_marker_dotplot.png, astro_tsne_labeled*.png
#
# Usage
#     source(file.path(SCORCH_REPO, "revision/tools/astro_cluster_features.R"))
#     res <- astro_features(scorch_astrocytes)          # object already in session
#     res$dotplot                                       # marker dotplot
#     head(res$suggested)                               # data-derived labels
#     head(res$markers)                                 # FindAllMarkers table
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(ggplot2); library(dplyr)
})

## Literature panels for the subtypes the Rmd names. Kept lineage-specific:
## a marker shared by two subtypes cannot separate them.
ASTRO_PANELS <- list(
  Homeostatic  = c("SLC1A2","SLC1A3","GLUL","GJA1","AQP4","MFGE8","NTSR2"),
  Reactive     = c("GFAP","VIM","CD44","SERPINA3","C3","OSMR","CHI3L1","TNC"),
  ## Reactive-PWH follows the paper's definition of the population
  ##   (GFAP+/CD44+/OSMR+/SERPINA3+), not interferon genes. It deliberately
  ##   overlaps the Reactive panel: Reactive and Reactive-PWH differ in degree,
  ##   not in separate markers.
  `Reactive-PWH` = c("GFAP","CD44","OSMR","SERPINA3","CHI3L1","TNC"),
  Signaling    = c("NTRK2","GRM3","ADGRV1","PTPRZ1","LRP1B","NRG3"),
  Synaptic     = c("SPARCL1","GPC5","NRXN1","NLGN3","EPHB1","CADM2"),
  Ciliated     = c("FOXJ1","PIFO","TMEM212","DNAH11","SPAG17","HYDIN"),
  Hypoxic      = c("HIF1A","VEGFA","SLC2A1","ADM","NDRG1","BNIP3"),
  `Micro-Astro`= c("P2RY12","CSF1R","C1QB","CX3CR1","DOCK8"),
  ## contaminant / doublet populations, present as clusters but not astrocytes
  Oligo        = c("MBP","MOG","MOBP","ST18","CNP"),
  Endo         = c("CLDN5","FLT1","PTPRB","ITM2A","ABCB1")
)

#' @param obj      the astrocyte Seurat object (split or joined; both handled)
#' @param group    metadata column holding the clustering
#' @param outdir   where CSVs/PNGs go; NULL = don't write, just return
#' @param min_pct,logfc thresholds passed to FindAllMarkers
#' @param n_top    markers per cluster to keep in the summary table
astro_features <- function(obj,
                           group  = "harmony_astro_clusters",
                           outdir = file.path(OUT_ROOT, "06-astrocyte_analysis"),
                           min_pct = 0.25, logfc = 0.25, n_top = 20) {

  if (!group %in% colnames(obj@meta.data)) {
    alt <- intersect(c("seurat_clusters","harmony_astro_clusters"), colnames(obj@meta.data))
    if (!length(alt)) stop("no clustering column found; pass group=")
    group <- alt[1]
    message("[astro] using clustering column '", group, "'")
  }
  Idents(obj) <- group
  n_clust <- length(unique(as.character(obj@meta.data[[group]])))
  message("[astro] ", ncol(obj), " cells in ", n_clust, " clusters")

  ## ---- 1. per-cluster marker genes ------------------------------------------
  ## FindAllMarkers needs a single joined layer; a split assay errors. Join a
  ## working copy so the caller's object keeps its layers for re-integration.
  work <- obj
  if (length(Layers(work[["RNA"]], search = "counts")) > 1) {
    message("[astro] joining layers for marker detection (working copy only)")
    work <- JoinLayers(work)
  }
  message("[astro] FindAllMarkers (presto Wilcoxon) ...")
  markers <- FindAllMarkers(work, only.pos = TRUE,
                            min.pct = min_pct, logfc.threshold = logfc)
  top <- markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = n_top, with_ties = FALSE) %>% ungroup()

  ## ---- 2. fraction of each cluster expressing each panel gene ---------------
  ## cluster_pct_expressing() is layer-aware and drops absent genes with a
  ## warning; combine = "max" scores a panel by its best-detected marker, since
  ## averaging a near-universal marker with a sparse one measures the panel's
  ## composition as much as the biology.
  all_genes <- unique(unlist(ASTRO_PANELS))
  gene_pct  <- cluster_gene_pct(obj, all_genes, group = group)
  panel_pct <- sapply(ASTRO_PANELS, function(gs) {
    g <- intersect(gs, rownames(gene_pct))
    if (!length(g)) return(rep(NA_real_, ncol(gene_pct)))
    apply(gene_pct[g, , drop = FALSE], 2, max)
  })
  rownames(panel_pct) <- colnames(gene_pct)

  ## relative to the purest cluster of each subtype, because detection rate
  ## differs enormously between panels (interferon genes are sparse in nuclei;
  ## SLC1A2 is near-universal), so raw fractions are not comparable across them
  rel <- sweep(panel_pct, 2, apply(panel_pct, 2, max, na.rm = TRUE), "/")

  n_by <- table(as.character(obj@meta.data[[group]]))
  suggested <- data.frame(
    cluster      = rownames(rel),
    n_cells      = as.integer(n_by[rownames(rel)]),
    suggested    = colnames(rel)[apply(rel, 1, which.max)],
    best_rel     = round(apply(rel, 1, max, na.rm = TRUE), 3),
    second       = colnames(rel)[apply(rel, 1, function(v) order(v, decreasing = TRUE)[2])],
    margin       = round(apply(rel, 1, function(v) {
                      s <- sort(v, decreasing = TRUE); s[1] - s[2] }), 3),
    stringsAsFactors = FALSE
  )
  suggested <- cbind(suggested, round(panel_pct, 3))
  suggested <- suggested[order(-suggested$n_cells), ]
  ## a small margin means the call is ambiguous - flag rather than hide
  suggested$confident <- suggested$margin >= 0.15

  ## ---- 3. dotplot over the curated panel ------------------------------------
  present <- all_genes[all_genes %in% rownames(obj)]
  missing <- setdiff(all_genes, present)
  if (length(missing))
    message("[astro] panel genes absent from the object: ", paste(missing, collapse = ", "))
  dp <- DotPlot(work, features = present, group.by = group, cluster.idents = TRUE) +
    RotatedAxis() +
    theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 9)) +
    labs(x = NULL, y = "cluster",
         title = "Astrocyte subtype markers by cluster",
         subtitle = "panels: Homeostatic | Reactive | HIV-Reactive | Signaling | Synaptic | Ciliated | Hypoxic | Micro-Astro | Oligo")

  ## dotplot of each cluster's own top markers, for the labels the panels miss
  dp_data <- DotPlot(work, group.by = group,
                     features = unique(top$gene[top$cluster %in% levels(Idents(work))])) +
    RotatedAxis() + theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1))

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write.csv(markers,   file.path(outdir, "astro_cluster_markers_all.csv"), row.names = FALSE)
    write.csv(top,       file.path(outdir, "astro_cluster_markers_top.csv"), row.names = FALSE)
    write.csv(suggested, file.path(outdir, "astro_cluster_suggested_labels.csv"), row.names = FALSE)
    ggsave(file.path(outdir, "astro_marker_dotplot.png"), dp,
           width = 14, height = 2 + 0.3 * n_clust, dpi = 200, limitsize = FALSE)
    ggsave(file.path(outdir, "astro_top_marker_dotplot.png"), dp_data,
           width = 26, height = 2 + 0.3 * n_clust, dpi = 200, limitsize = FALSE)
    message("[astro] wrote 3 CSVs + 2 PNGs to ", outdir)
  }

  message("\n[astro] suggested labels (check `confident`; margin < 0.15 is ambiguous):")
  print(suggested[, c("cluster","n_cells","suggested","best_rel","second","margin","confident")],
        row.names = FALSE)
  n_amb <- sum(!suggested$confident)
  if (n_amb) message("[astro] ", n_amb, " cluster(s) are ambiguous - inspect the dotplot before labelling")

  invisible(list(markers = markers, top = top, suggested = suggested,
                 panel_pct = panel_pct, rel = rel, gene_pct = gene_pct,
                 dotplot = dp, top_dotplot = dp_data))
}

# ==============================================================================
# Plot helpers — a panel-annotated dotplot and a labelled tSNE.
# ==============================================================================

#' Dotplot with genes grouped and labelled by the subtype panel they belong to.
#'
#' A flat 50-gene axis is unreadable and, worse, invites mis-reading: without the
#' grouping you cannot see that a cluster is high on Reactive markers and low on
#' everything else. Facetting by panel makes the block structure the point.
#'
#' @param obj    astrocyte object
#' @param group  clustering or label column to put on the y axis
#' @param outfile PNG path; NULL returns the plot without writing
astro_dotplot <- function(obj,
                          group  = "harmony_astro_clusters",
                          panels = ASTRO_PANELS,
                          outfile = file.path(OUT_ROOT, "06-astrocyte_analysis",
                                              "astro_marker_dotplot_annotated.png"),
                          width = 20, height = NULL) {
  if (!group %in% colnames(obj@meta.data)) stop("no such column: ", group)
  work <- obj
  if (length(Layers(work[["RNA"]], search = "counts")) > 1) work <- JoinLayers(work)
  Idents(work) <- group

  ## gene -> panel, keeping panel order and dropping genes absent from the object
  g2p <- unlist(lapply(names(panels), function(p) setNames(rep(p, length(panels[[p]])), panels[[p]])))
  g2p <- g2p[names(g2p) %in% rownames(work)]
  g2p <- g2p[!duplicated(names(g2p))]
  if (!length(g2p)) stop("none of the panel genes are in the object")
  absent <- setdiff(unlist(panels), names(g2p))
  if (length(absent)) message("[astro] absent from object, dropped: ", paste(absent, collapse = ", "))

  d <- DotPlot(work, features = names(g2p), group.by = group)$data
  d$panel <- factor(g2p[as.character(d$features.plot)], levels = names(panels))
  d$features.plot <- factor(as.character(d$features.plot), levels = names(g2p))

  p <- ggplot(d, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
    facet_grid(~ panel, scales = "free_x", space = "free_x", switch = "x") +
    scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B",
                           midpoint = 0, name = "scaled\nexpression") +
    scale_size(range = c(0, 6), name = "% expressing") +
    labs(x = NULL, y = NULL,
         title = "Astrocyte subtype markers by cluster",
         subtitle = "genes grouped by the subtype panel they mark") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 9),
          panel.spacing = unit(2, "pt"),
          strip.placement = "outside",
          strip.background = element_rect(fill = "grey90", colour = NA),
          strip.text.x = element_text(size = 8, face = "bold", angle = 90))

  if (!is.null(outfile)) {
    if (is.null(height)) height <- 2.5 + 0.28 * length(unique(d$id))
    ggsave(outfile, p, width = width, height = height, dpi = 200, limitsize = FALSE)
    message("[astro] dotplot -> ", outfile)
  }
  p
}

#' Labelled tSNE (or UMAP) of the astrocyte subclusters.
#'
#' @param label_col what to colour and label by: cluster ids, or subtype labels
#'   once assigned ("subcluster_labels" uses the shared subtype palette)
astro_tsne <- function(obj,
                       label_col = "harmony_astro_clusters",
                       reduction = NULL,
                       outfile = file.path(OUT_ROOT, "06-astrocyte_analysis",
                                           "astro_tsne_labeled.png"),
                       split_by = NULL, pt.size = 0.4,
                       label = TRUE, label.size = 9, label.box = TRUE,
                       point.alpha = 0.45, width = 12, height = 10, cols = NULL) {
  ## Coloured by subtype: use the shared Fig1D palette (ASTRO_SUBTYPE_COLOURS,
  ##   pastel_blue) plus greys for the Oligo / Endo-Astro contamination clusters,
  ##   so the supplementals match Fig1D. Cluster-number colouring, or explicit
  ##   cols, is left alone.
  if (is.null(cols) && label_col == "subcluster_labels") {
    if (!exists("ASTRO_SUBTYPE_COLOURS"))
      source(file.path(SCORCH_REPO, "revision/tools/astro_reviewer_figures.R"))
    pal  <- c(ASTRO_SUBTYPE_COLOURS, ASTRO_EXTRA_COLOURS)
    labs_present <- unique(as.character(obj@meta.data[[label_col]]))
    cols <- pal[intersect(names(pal), labs_present)]
    if (length(setdiff(labs_present, names(pal)))) cols <- NULL   # unknown labels: fall back
    else obj@meta.data[[label_col]] <- factor(as.character(obj@meta.data[[label_col]]), levels = names(cols))
  }
  if (is.null(reduction)) {
    reduction <- intersect(c("tsne","umap"), names(obj@reductions))
    if (!length(reduction)) stop("no tsne/umap reduction on the object")
    reduction <- reduction[1]
  }
  if (!label_col %in% colnames(obj@meta.data)) stop("no such column: ", label_col)

  ## PWH / PWoH for display. condition_new stores "HIV"/"CTR" and the DE design
  ## depends on those values, so only the shown labels change (pwh_relabel()).
  if (!is.null(split_by) && split_by == "condition_new") {
    obj$.pwh <- pwh_relabel(obj@meta.data[[split_by]])
    split_by <- ".pwh"
  }

  n <- table(as.character(obj@meta.data[[label_col]]))
  message("[astro] ", reduction, " coloured by ", label_col, " (",
          length(n), " groups, ", sum(n), " cells)")

  ## Readable labels on a dense tSNE: DimPlot's own label geom is small and thin,
  ## so labels are drawn separately with LabelClusters(), which takes fontface
  ## and a background box. Points are dimmed so the text carries the figure;
  ## alpha is set on the layer afterwards because DimPlot() does not expose it.
  p <- DimPlot(obj, reduction = reduction, group.by = label_col,
               label = FALSE, repel = TRUE,
               raster = FALSE, shuffle = TRUE, pt.size = pt.size,
               split.by = split_by) +
    labs(title = paste0("Astrocyte subclusters (", reduction, ")"),
         subtitle = paste0(sum(n), " cells, ", length(n), " groups")) +
    theme(plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 12),
          legend.text = element_text(size = 13),
          legend.key.size = grid::unit(1.1, "lines"))

  if (label) {
    if (!is.null(point.alpha) && length(p$layers))
      p$layers[[1]]$aes_params$alpha <- point.alpha
    p <- LabelClusters(p, id = label_col, repel = TRUE,
                       size = label.size, fontface = "bold",
                       colour = "black",
                       box = label.box,
                       fill = if (label.box) scales::alpha("white", 0.75) else NA,
                       segment.colour = "grey30", segment.size = 0.3,
                       max.overlaps = Inf)
  }

  if (!is.null(outfile)) {
    w <- if (is.null(split_by)) width else width * length(unique(obj@meta.data[[split_by]]))
    ggsave(outfile, p, width = w, height = height, dpi = 300, limitsize = FALSE)
    message("[astro] tSNE -> ", outfile)
  }
  p
}

#' Write both plots plus a per-cluster count table in one call.
astro_save_plots <- function(obj, group = "harmony_astro_clusters",
                             outdir = file.path(OUT_ROOT, "06-astrocyte_analysis")) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dp <- astro_dotplot(obj, group = group,
                      outfile = file.path(outdir, "astro_marker_dotplot_annotated.png"))
  tp <- astro_tsne(obj, label_col = group,
                   outfile = file.path(outdir, "astro_tsne_labeled.png"))
  tc <- if ("condition_new" %in% colnames(obj@meta.data))
    astro_tsne(obj, label_col = group, split_by = "condition_new",
               outfile = file.path(outdir, "astro_tsne_labeled_by_condition.png")) else NULL
  cnt <- as.data.frame(table(cluster = as.character(obj@meta.data[[group]])))
  names(cnt)[2] <- "n_cells"
  write.csv(cnt, file.path(outdir, "astro_cluster_sizes.csv"), row.names = FALSE)
  invisible(list(dotplot = dp, tsne = tp, tsne_by_condition = tc, counts = cnt))
}

# ==============================================================================
# Endfeet vs endothelial contamination vs doublet.
#
# A vascular-looking cluster inside an astrocyte subset has three explanations,
# and they need different actions:
#
#   endothelial contamination  CLDN5/FLT1 high, AQP4/SLC1A2 LOW    -> drop
#   perivascular / endfeet     CLDN5/FLT1 low,  AQP4 HIGH          -> KEEP, real astrocytes
#   astro-endo doublet         both high, and nCount/nFeature inflated -> drop
#
# Endfeet are AQP4-enriched (AQP4 is polarised to the endfoot membrane and
# anchored by the dystrophin complex - DTNA, SNTA1, DAG1), so a genuine endfeet
# population scores HIGH on astrocyte markers, not low. Markers alone cannot
# separate endfeet from doublets, since both co-express; library size is what
# distinguishes them, so it is reported here too.
VASC_PANELS <- list(
  Endothelial = c("CLDN5","FLT1","PECAM1","VWF","ABCB1","PTPRB"),
  Endfeet     = c("AQP4","MLC1","DTNA","SNTA1","DAG1","ATP1A2","GJA1"),
  Astrocyte   = c("SLC1A2","SLC1A3","GLUL","NTSR2"),
  Mural       = c("PDGFRB","RGS5","ACTA2","NOTCH3")   # pericyte/SMC, a 4th option
)

#' @param clusters which clusters to examine; default = all
astro_vascular_check <- function(obj, clusters = NULL,
                                 group = "harmony_astro_clusters") {
  if (!group %in% colnames(obj@meta.data)) stop("no such column: ", group)
  g <- as.character(obj@meta.data[[group]])
  if (is.null(clusters)) clusters <- sort(unique(g))
  clusters <- as.character(clusters)

  genes <- unique(unlist(VASC_PANELS))
  gp <- cluster_gene_pct(obj, genes, group = group)

  score <- sapply(VASC_PANELS, function(gs) {
    gg <- intersect(gs, rownames(gp))
    if (!length(gg)) return(rep(NA_real_, ncol(gp)))
    apply(gp[gg, , drop = FALSE], 2, max)
  })
  rownames(score) <- colnames(gp)

  ## library size — inflated counts are the doublet signature
  qc <- data.frame(
    cluster  = names(tapply(obj$nCount_RNA, g, median)),
    n_cells  = as.integer(table(g)[names(tapply(obj$nCount_RNA, g, median))]),
    med_nCount   = as.numeric(tapply(obj$nCount_RNA, g, median)),
    med_nFeature = as.numeric(tapply(obj$nFeature_RNA, g, median)),
    stringsAsFactors = FALSE
  )
  overall_nCount <- median(obj$nCount_RNA)
  qc$nCount_vs_median <- round(qc$med_nCount / overall_nCount, 2)

  out <- merge(qc, data.frame(cluster = rownames(score), round(score, 3),
                              check.names = FALSE), by = "cluster")
  out <- out[out$cluster %in% clusters, ]
  out$call <- with(out, ifelse(
    Endothelial > 0.30 & Astrocyte < 0.30, "endothelial - drop",
    ifelse(Endothelial > 0.30 & Astrocyte >= 0.30 & nCount_vs_median > 1.4, "likely doublet - drop",
    ifelse(Endothelial <= 0.30 & Endfeet >= 0.50, "endfeet/perivascular astrocyte - KEEP",
    ifelse(Mural > 0.30 & Astrocyte < 0.30, "mural (pericyte/SMC) - drop", "ambiguous - inspect")))))
  out <- out[order(-out$Endothelial), ]

  message("[vasc] median nCount across all astrocytes: ", round(overall_nCount))
  message("[vasc] thresholds are heuristics - read the numbers, not the call:")
  print(out[, c("cluster","n_cells","Endothelial","Endfeet","Astrocyte","Mural",
                "med_nCount","nCount_vs_median","call")], row.names = FALSE)
  message("\n[vasc] per-gene detail:")
  print(round(t(gp[, clusters, drop = FALSE]), 3))
  invisible(out)
}

# ==============================================================================
# Colour schemes for FeaturePlots of sparse nuclear expression.
#
# The governing constraint is zero-inflation: in snRNA-seq most cells do not
# detect a given gene, so whatever colour sits at the low end covers most of the
# panel. Plain viridis is perceptually uniform but its low end is dark purple,
# which reads as signal - the zeros end up looking like the strongest feature of
# the plot. Every scheme here keeps the floor near-white for that reason.
#
#   feature_scale("grey_red")    Seurat-like but warmer; safe default
#   feature_scale("heat")        white -> amber -> crimson -> plum; good dynamic range
#   feature_scale("blues")       cool, high contrast, prints well in greyscale
#   feature_scale("viridis_hi")  viridis with the dark end trimmed off
#   feature_scale("magma_hi")    magma trimmed; strong for a single hero gene
#
# USE:
#   FeaturePlot(obj, features = g, reduction = "tsne", order = TRUE,
#               combine = FALSE) |> lapply(\(p) p + feature_scale("heat"))
#
# order = TRUE matters as much as the palette: it draws expressing cells on top
# of the zeros instead of letting them be buried.
FEATURE_PALETTES <- list(
  grey_red   = c("grey88", "#FDBB84", "#E34A33", "#B30000"),
  heat       = c("grey95", "#FFC300", "#FF5733", "#C70039", "#581845"),
  blues      = c("grey92", "#9ECAE1", "#4292C6", "#08519C", "#08306B"),
  viridis_hi = c("grey93", "#3EBC96", "#35B779", "#B5DE2B", "#FDE725"),
  magma_hi   = c("grey93", "#FEB078", "#F1605D", "#B63679", "#721F81")
)

feature_scale <- function(palette = "heat", limits = NULL, name = "expr") {
  cols <- FEATURE_PALETTES[[palette]]
  if (is.null(cols))
    stop("unknown palette '", palette, "'. Options: ",
         paste(names(FEATURE_PALETTES), collapse = ", "))
  ggplot2::scale_colour_gradientn(colours = cols, limits = limits, name = name,
                                  na.value = "grey93")
}

#' Draw the same gene under every palette, so the choice is made by looking
#' rather than by description.
astro_palette_preview <- function(obj, feature, reduction = "tsne",
                                  outfile = file.path(OUT_ROOT, "06-astrocyte_analysis",
                                                      "figures", "feature_palette_preview.png")) {
  ps <- lapply(names(FEATURE_PALETTES), function(nm) {
    FeaturePlot(obj, features = feature, reduction = reduction, order = TRUE,
                pt.size = 0.3, raster = FALSE) +
      feature_scale(nm) + ggtitle(nm) +
      theme_void() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
                           legend.position = "none")
  })
  p <- patchwork::wrap_plots(ps, ncol = 3) +
    patchwork::plot_annotation(title = paste0("Palette options — ", feature))
  if (!is.null(outfile)) {
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    ggsave(outfile, p, width = 15, height = 10, dpi = 200)
    message("[fig] palette preview -> ", outfile)
  }
  p
}

# ==============================================================================
# Restrict an enrichResult to a subset of genes FOR DISPLAY only.
#
# A heatplot of 900 genes is unreadable, but shrinking the input list to make the
# picture legible changes the statistics: enrichment p-values depend on the size
# and composition of the gene set. So compute enrichment on the full list, then
# trim what is drawn.
#
#   ego  <- enrichGO(gene = gene, universe = sde$de$gene, ...)   # all 900
#   show <- head(gene[order(-lfc[gene])], 60)                    # top 60 by logFC
#   heatplot(trim_enrich(ego, show), foldChange = geneList, showCategory = 10)
#
# The p-values and term selection are unchanged; only the gene axis is thinned.
# Terms left with no displayed genes are dropped so the panel has no empty rows.
#' @param n_terms keep only the n most significant terms. Do this HERE rather
#'   than relying on heatplot's showCategory: trimming genes can empty a term,
#'   which makes showCategory fall through to the next term down and quietly
#'   change which pathways are displayed. Capping the terms first makes the panel
#'   deterministic.
#' @param keep_genes genes to draw. If NULL, every gene of the kept terms.
trim_enrich <- function(ego, keep_genes = NULL, n_terms = NULL, drop_empty = TRUE) {
  r <- ego@result
  r <- r[order(r$p.adjust), , drop = FALSE]
  if (!is.null(n_terms)) r <- utils::head(r, n_terms)
  if (!is.null(keep_genes)) {
    r$geneID <- vapply(strsplit(r$geneID, "/", fixed = TRUE), function(g)
      paste(intersect(g, keep_genes), collapse = "/"), character(1))
    r$Count <- lengths(strsplit(r$geneID, "/", fixed = TRUE))
    if (drop_empty) r <- r[nzchar(r$geneID), , drop = FALSE]
  }
  ego@result <- r
  ego
}

#' Genes belonging to the top n terms of an enrichResult, most significant first.
#' Selecting display genes from the terms actually shown keeps the heatmap dense
#' - picking globally-top genes leaves many of them in no displayed pathway.
genes_in_top_terms <- function(ego, n_terms = 10, de = NULL, n_genes = NULL) {
  r <- ego@result[order(ego@result$p.adjust), , drop = FALSE]
  g <- unique(unlist(strsplit(utils::head(r$geneID, n_terms), "/", fixed = TRUE)))
  if (!is.null(de)) {
    d <- de[de$gene %in% g, , drop = FALSE]
    d <- d[order(d$FDR, -d$logFC), ]
    g <- d$gene
  }
  if (!is.null(n_genes)) g <- utils::head(g, n_genes)
  g
}

#' Top n genes of a DE table for DISPLAY.
#'
#' by = "significance" (default) ranks on FDR; by = "logFC" ranks on effect size.
#'
#' Significance is the better default and the reason is concrete: ranking the
#' Reactive-PWH markers by logFC pushed the complement genes out of the top 60
#' (C3 rank 259, C4A 392, SERPING1 455, C1R 494) even though SERPING1 has an FDR
#' of 1.5e-09 - among the most reliably changed genes in the test. Effect-size
#' ranking favours noisy genes with large fold-changes over consistently
#' regulated ones, which is the opposite of what a marker panel should show.
top_marker_genes <- function(de, n = 60, fdr = 0.05, min_lfc = 0.5,
                             by = c("significance", "logFC")) {
  by <- match.arg(by)
  d <- de[!is.na(de$FDR) & de$FDR < fdr & de$logFC > min_lfc, ]
  d <- if (by == "significance") d[order(d$FDR, -d$logFC), ] else d[order(-d$logFC), ]
  utils::head(d$gene, n)
}

# ==============================================================================
# Marker figure grouped by DATABASE-DERIVED terms, not by hand.
#
# Genes are assigned to GO Biological Process terms from org.Hs.eg.db via
# clusterProfiler::enrichGO. Which terms appear, and which genes belong to them,
# come from the ontology and the enrichment statistics - not from the analyst.
# Nothing here invents a grouping or a label.
#
#   ego <- enrichGO(gene = up_genes, universe = sde$de$gene, OrgDb = org.Hs.eg.db,
#                   keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
#                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
#   astro_term_module_plot(ego, sde$de, n_terms = 8)
#
# Each facet is a GO term with its adjusted p-value; genes are those the ontology
# assigns to it, restricted to the significant ones and ranked by fold-change.
astro_term_module_plot <- function(ego, de, n_terms = 8, genes_per_term = 12,
                                   fdr = 0.05, min_lfc = 0,
                                   title = NULL,
                                   outfile = file.path(OUT_ROOT, "06-astrocyte_analysis",
                                                       "figures", "marker_terms.png")) {
  r <- ego@result[order(ego@result$p.adjust), , drop = FALSE]
  r <- utils::head(r, n_terms)
  if (!nrow(r)) stop("no enriched terms")

  d <- de[!is.na(de$FDR) & de$FDR < fdr & de$logFC > min_lfc, c("gene","logFC","FDR")]

  long <- do.call(rbind, lapply(seq_len(nrow(r)), function(i) {
    g <- strsplit(r$geneID[i], "/", fixed = TRUE)[[1]]
    x <- d[d$gene %in% g, , drop = FALSE]
    if (!nrow(x)) return(NULL)
    x <- x[order(-x$logFC), ]
    x <- utils::head(x, genes_per_term)
    ## label carries the term's own adjusted p-value, so the panel is auditable
    x$term <- sprintf("%s\n(p.adj %.1e)", r$Description[i], r$p.adjust[i])
    x$term_order <- i
    x
  }))
  if (is.null(long)) stop("no significant genes in the top terms")

  long$term <- factor(long$term, levels = unique(long$term[order(long$term_order)]))
  long <- long[order(long$term, long$logFC), ]
  long$gene_lab <- factor(paste(long$term, long$gene), levels = unique(paste(long$term, long$gene)))
  long$neglog10FDR <- pmin(-log10(long$FDR), 20)

  p <- ggplot(long, aes(x = logFC, y = gene_lab)) +
    geom_segment(aes(x = 0, xend = logFC, yend = gene_lab), colour = "grey75", linewidth = 0.4) +
    geom_point(aes(colour = neglog10FDR, size = neglog10FDR)) +
    scale_y_discrete(labels = function(x) sub("^.*\\s", "", x)) +
    facet_grid(term ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_gradientn(colours = c("#FDBB84", "#E34A33", "#7F0000"),
                           name = expression(-log[10]~FDR)) +
    scale_size(range = c(2, 6), guide = "none") +
    labs(x = expression(log[2]~"fold-change vs the same donor's other astrocytes"),
         y = NULL,
         title = title %||% "Subtype markers by enriched GO Biological Process",
         subtitle = "terms and gene assignments from GO via clusterProfiler; genes shown are significant at the stated FDR") +
    theme_bw(base_size = 11) +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 8, hjust = 0),
          strip.background = element_rect(fill = "grey93", colour = NA),
          panel.grid.major.y = element_blank())

  if (!is.null(outfile)) {
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    h <- 2 + 0.2 * nrow(long)
    ggsave(outfile, p, width = 9, height = h, dpi = 300, limitsize = FALSE)
    ggsave(sub("\\.png$", ".pdf", outfile), p, width = 9, height = h, limitsize = FALSE)
    message("[term] ", nrow(long), " genes across ", nlevels(long$term), " GO terms -> ", outfile)
  }
  p
}

# ==============================================================================
# Reactome enrichment.
#
# Reactome is a curated pathway database with a different structure to GO: fewer,
# more specific terms and far less of GO BP's near-duplicate redundancy, which
# usually makes the resulting figure easier to read.
#
# It is keyed on Entrez ids, not symbols, so the mapping is done here and the
# genes that fail to map are REPORTED rather than silently dropped - symbol-to-
# Entrez loss is a common and invisible source of shrinkage in these analyses.
#
#   rp <- astro_reactome(gene, universe = sde$de$gene)
#   dotplot(rp$enrich, showCategory = 10)
#   astro_term_module_plot(rp$enrich, sde$de, n_terms = 8)
astro_reactome <- function(genes, universe = NULL, pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05, readable = TRUE) {
  suppressPackageStartupMessages({
    library(ReactomePA); library(org.Hs.eg.db); library(AnnotationDbi)
  })
  map <- function(x) {
    e <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(x),
           keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"))
    e[!is.na(e)]
  }
  ge <- map(genes)
  message("[reactome] gene set: ", length(ge), " of ", length(unique(genes)),
          " symbols mapped to Entrez",
          if (length(ge) < length(unique(genes)))
            paste0(" (", length(unique(genes)) - length(ge), " unmapped)") else "")
  ue <- if (!is.null(universe)) map(universe) else NULL
  if (!is.null(ue)) message("[reactome] universe: ", length(ue), " of ",
                            length(unique(universe)), " mapped")

  res <- enrichPathway(gene = unname(ge), universe = if (is.null(ue)) NULL else unname(ue),
                       organism = "human", pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff, readable = readable)
  n <- if (is.null(res)) 0 else sum(res@result$p.adjust < pvalueCutoff, na.rm = TRUE)
  message("[reactome] enriched pathways at p.adj < ", pvalueCutoff, ": ", n)
  list(enrich = res, entrez = ge, universe_entrez = ue)
}
