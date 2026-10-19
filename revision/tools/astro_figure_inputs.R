# ==============================================================================
# Load the inputs that the figure chunks of 06-astrocyte_analysis/astrocyte_analysis.Rmd
# expect, plus figure-saving and palette helpers (function library).
#
# Method
#   The published figure chunks reference objects that were session state and are
#   never created by the Rmd itself. This provides them:
#   de_edgeR       pseudobulk DE results (cluster_id, logFC, p_adj.glb, avg_CTR_frq,
#                  avg_HIV_frq); by default the cross-region table returned by
#                  xr_de_file() (revision/config/paths.R), on uncorrected counts
#   genes          the enrichGO() universe: genes actually tested (or detected in
#                  >5% of nuclei in either group), never the whole genome
#   scorch_meta    combined-object metadata for all cell types, needed for each
#                  sample's astrocyte fraction
#   heatmap_list   empty accumulator the S15 loop writes into
#   deg_results    per-cell-type tables in FindMarkers column names, rebuilt from
#                  the donor-level pseudobulk DE (deg_results_from_de) instead of
#                  the cell-level Wilcoxon tests used in the published analysis
#   signac_ins is not loaded here.
#
# Inputs   OUT_ROOT/03-DE_analysis/... (cross-region or <REGION>/results DE tables)
#          OUT_ROOT/01-Regional_data_processing/combined/data/combined_v3_cell_metadata.csv
#          (written by revision/tools/dump_combined_metadata.R)
# Outputs  figures via save_fig()/with_fig() under OUT_ROOT/06-astrocyte_analysis/figures/
#
# Usage
#   source("revision/tools/astro_figure_inputs.R")
#   inp <- astro_figure_inputs()
#   de_edgeR <- inp$de_edgeR; genes <- inp$genes
#   scorch_meta <- inp$scorch_meta; heatmap_list <- inp$heatmap_list
#   deg_results <- inp$deg_results
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr) })

#' @param source_de "cross_region" (default) or "INS" / "PFC" / "VST"
#'
#'   The astrocyte object spans all three regions, so the cross-region result is
#'   the matching contrast. Per-region results are available if a figure is meant
#'   to be region-specific - state which in the legend either way, because
#'   "astrocyte DE" means different things across the two.
#'
#' @param counts "uncorrected" (default, the primary analysis) or "corrected"
astro_figure_inputs <- function(source_de = "cross_region",
                                counts = c("uncorrected", "corrected"),
                                universe = c("tested", "detected")) {
  counts   <- match.arg(counts)
  universe <- match.arg(universe)
  sfx <- if (counts == "corrected") "_corrected" else ""
  DE <- file.path(OUT_ROOT, "03-DE_analysis")

  ## xr_de_file(): the edgeR table, or the donor-aware one under SCORCH_XR_DE=donor_aware
  ## (revision/config/paths.R). Callers add XR_OUT_SFX to their output names.
  p <- if (identical(source_de, "cross_region"))
    xr_de_file(sfx)
  else
    file.path(DE, source_de, "results",
              sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq%s.RDS", source_de, sfx))
  if (!file.exists(p)) stop("missing DE result: ", p)

  de_edgeR <- as.data.frame(readRDS(p))
  need <- c("gene","cluster_id","logFC","logCPM","p_adj.glb","avg_CTR_frq","avg_HIV_frq")
  miss <- setdiff(need, colnames(de_edgeR))
  if (length(miss)) stop("DE result lacks: ", paste(miss, collapse = ", "))
  message("[inp] de_edgeR: ", source_de, " (", counts, ") - ", nrow(de_edgeR),
          " rows, ", length(unique(de_edgeR$cluster_id)), " cell types")

  ## ---- enrichGO universe ---------------------------------------------------
  ## The genes that were actually TESTED, not every gene in the annotation.
  ## Using a genome-wide universe against a DE-derived gene list inflates
  ## enrichment, because the background includes genes that could never have
  ## been called.
  genes <- if (universe == "tested") unique(de_edgeR$gene) else
    unique(de_edgeR$gene[de_edgeR$avg_CTR_frq > 0.05 | de_edgeR$avg_HIV_frq > 0.05])
  message("[inp] genes (enrichGO universe, '", universe, "'): ", length(genes))

  ## ---- scorch_meta ---------------------------------------------------------
  mp <- file.path(OUT_ROOT, "01-Regional_data_processing", "combined", "data",
                  "combined_v3_cell_metadata.csv")
  scorch_meta <- NULL
  if (file.exists(mp)) {
    scorch_meta <- utils::read.csv(mp, row.names = 1, stringsAsFactors = FALSE)
    message("[inp] scorch_meta: ", nrow(scorch_meta), " cells x ",
            ncol(scorch_meta), " columns")
  } else {
    warning("combined_v3_cell_metadata.csv not found - run ",
            "revision/docker/run.sh revision/tools/dump_combined_metadata.R first. ",
            "scorch_meta returned as NULL.")
  }

  list(de_edgeR = de_edgeR, genes = genes, scorch_meta = scorch_meta,
       heatmap_list = list(), de_path = p,
       ## same shape as the old cell-level object, donor-level statistics inside
       deg_results = deg_results_from_de(de_edgeR))
}

# ==============================================================================
# save_fig() — one call to write a figure with consistent naming and size.
#
# Most figure chunks in the astrocyte Rmd end by printing a plot, so nothing
# reaches disk. This writes both
# a PNG (for looking at) and a PDF (vector, for the manuscript) under
# 06-astrocyte_analysis/figures/, and returns the plot so it still displays
# inline in the notebook.
#
#   p <- ggplot(...) + ...
#   save_fig(p, "Fig1B_GO_PWoH_astrocytes")
#
# For a base-graphics figure (ComplexHeatmap, cnetplot via print) use
#   with_fig("S15_pathway_heatmap", { print(p) })
# ==============================================================================

ASTRO_FIG_DIR <- function() file.path(OUT_ROOT, "06-astrocyte_analysis", "figures")

save_fig <- function(p, name, width = 9, height = 7, dpi = 300, pdf = TRUE) {
  d <- ASTRO_FIG_DIR()
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  name <- gsub("[^A-Za-z0-9._-]+", "_", name)
  ggplot2::ggsave(file.path(d, paste0(name, ".png")), p,
                  width = width, height = height, dpi = dpi, limitsize = FALSE)
  if (pdf)
    ggplot2::ggsave(file.path(d, paste0(name, ".pdf")), p,
                    width = width, height = height, limitsize = FALSE)
  message("[fig] ", name, " -> ", d)
  invisible(p)
}

#' Base-graphics equivalent: opens a device, evaluates the expression, closes it.
with_fig <- function(name, expr, width = 2400, height = 2000, res = 300) {
  d <- ASTRO_FIG_DIR()
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  name <- gsub("[^A-Za-z0-9._-]+", "_", name)
  grDevices::png(file.path(d, paste0(name, ".png")), width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  message("[fig] ", name, " -> ", d)
}

#' Save every element of a named list of plots (e.g. the S15 heatmap_list).
save_fig_list <- function(lst, prefix, ...) {
  if (!length(lst)) { message("[fig] ", prefix, ": empty list, nothing saved"); return(invisible(NULL)) }
  nm <- names(lst)
  if (is.null(nm)) nm <- as.character(seq_along(lst))
  for (i in seq_along(lst)) {
    p <- lst[[i]]
    if (inherits(p, "ggplot")) save_fig(p, paste0(prefix, "_", nm[i]), ...)
    else message("[fig] ", prefix, "_", nm[i], ": not a ggplot, skipped")
  }
  invisible(lst)
}

# ==============================================================================
# deg_results, rebuilt from donor-level pseudobulk.
#
# The S15 loops expect a named list (one element per cell type) of data frames
# with Seurat's FindMarkers column names and genes as rownames:
#     p_val_adj, avg_log2FC, pct.1
# The published object came from cell-level Wilcoxon tests, which are not
# donor-aware and have inflated type I error. This builds the same shape from the
# pseudobulk edgeR results, so the loop code runs unchanged while the statistics
# underneath are donor-level (for subtype-vs-rest comparisons, see
# revision/tools/astro_subtype_de_donor.R):
#
#     p_val_adj  <- p_adj.glb        (donor-level FDR)
#     avg_log2FC <- logFC            (pseudobulk log2 fold-change)
#     pct.1      <- avg_HIV_frq      (fraction expressing in PWH)
#     pct.2      <- avg_CTR_frq      (fraction expressing in PWoH)
#
# The filters in those loops (p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.05)
# keep their meaning, but the logFC scale differs between cell-level and
# pseudobulk tests, so a cutoff tuned on one is not equivalent on the other.
deg_results_from_de <- function(de_edgeR, min_frq = NULL) {
  stopifnot(all(c("gene","cluster_id","logFC","p_adj.glb") %in% colnames(de_edgeR)))
  d <- as.data.frame(de_edgeR)
  if (!is.null(min_frq))
    d <- d[(d$avg_CTR_frq > min_frq) | (d$avg_HIV_frq > min_frq), , drop = FALSE]
  cts <- sort(unique(as.character(d$cluster_id)))
  out <- lapply(cts, function(ct) {
    x <- d[d$cluster_id == ct, , drop = FALSE]
    x <- x[!duplicated(x$gene), , drop = FALSE]
    r <- data.frame(
      avg_log2FC = x$logFC,
      p_val_adj  = x$p_adj.glb,
      pct.1      = if ("avg_HIV_frq" %in% colnames(x)) x$avg_HIV_frq else NA_real_,
      pct.2      = if ("avg_CTR_frq" %in% colnames(x)) x$avg_CTR_frq else NA_real_,
      row.names  = x$gene, stringsAsFactors = FALSE)
    r[order(r$p_val_adj), , drop = FALSE]
  })
  names(out) <- cts
  message("[inp] deg_results rebuilt from donor-level DE: ", length(out), " cell types (",
          paste(utils::head(cts, 6), collapse = ", "),
          if (length(cts) > 6) ", ..." else "", ")")
  out
}

# ==============================================================================
# Donor palette: shades of blue for PWoH, shades of orange for PWH.
#
# A stacked bar coloured by donor is unreadable with a categorical rainbow - you
# cannot see condition at all. Ramping within each condition keeps every donor
# distinguishable while making the PWoH/PWH split visible at a glance, and it
# matches PWH_COLOURS used elsewhere so the figures agree.
#
#   color_map <- donor_condition_palette(metadata$sample_ids, metadata$condition_new)
#   ... + scale_fill_manual(values = color_map)
donor_condition_palette <- function(donors, conditions,
                                    pwoh = c("#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                                    pwh  = c("#FDD0A2", "#FD8D3C", "#D94801", "#7F2704")) {
  d <- data.frame(donor = as.character(donors),
                  cond  = as.character(pwh_relabel(conditions)),
                  stringsAsFactors = FALSE)
  d <- d[!duplicated(d$donor), ]
  d <- d[order(d$cond, d$donor), ]
  out <- character(0)
  for (cc in unique(d$cond)) {
    ids <- d$donor[d$cond == cc]
    ramp <- if (cc == "PWH") pwh else pwoh
    cols <- grDevices::colorRampPalette(ramp)(max(length(ids), 2))[seq_along(ids)]
    out <- c(out, setNames(cols, ids))
  }
  message("[pal] ", sum(d$cond == "PWoH"), " PWoH (blues), ",
          sum(d$cond == "PWH"), " PWH (oranges)")
  out
}
