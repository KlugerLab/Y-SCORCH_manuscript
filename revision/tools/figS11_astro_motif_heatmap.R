# ==============================================================================
# Figure S11: astrocyte motif activity heatmaps, rebuilt on the rerun chromVAR results.
#
# STRUCTURE AS PUBLISHED (04-ATAC_DE/regional/02-motif_analysis_new.R):
#   a  motifs over-represented in DA peaks UP in PWH
#   b  motifs over-represented in DA peaks DOWN in PWH
#   columns  region x condition (6), rows z-scored across those columns,
#   rows k-means clustered (row_km = 5).
#
# REBUILT FROM THE RERUN, not from the published objects:
#   motif choice  FindMotifs over-representation (p.adjust < 0.05) on each region's
#                 DA peaks, kept if over-represented in >= --min-regions regions
#                 (default all 3) - tables already written by
#                 revision/tools/atac_motif_chromvar.R as <CT>_findmotifs_{up,down}.csv
#   activity      per-cell chromVAR z-scores from <CT>_chromvar_cells.rds, averaged
#                 within region x condition. The published run used Signac's
#                 RunChromVAR scale.data; Signac 1.17 dropped RunChromVAR, so the
#                 z-scores come from chromVAR directly (same pipeline, see
#                 atac_motif_chromvar.R).
#
# ADDED  a right-hand marker for motifs that are significant in the DONOR-AWARE
#   pooled model (revision/tools/atac_motif_crossregion.R, centered variant, FDR < 0.05).
#   The published stars came from FindMotifs over-representation, which is not a
#   PWH-vs-PWoH test and is not donor-aware.
#
# INPUT   OUT_ROOT/04-ATAC_DE/motif_chromvar/<R>/<CT>_chromvar_cells.rds, <CT>_findmotifs_{up,down}.csv
#         OUT_ROOT/04-ATAC_DE/motif_chromvar/crossregion/<CT>_chromvar_crossregion.csv
# OUTPUT  OUT_ROOT/04-ATAC_DE/figures/figS11/
#           figS11_<celltype>_motif_activity_{up,down}.png/.pdf
#           figS11_<celltype>_motif_activity_{up,down}.csv   every value drawn
#
# Usage: revision/docker/run.sh revision/tools/figS11_astro_motif_heatmap.R
#          [--celltypes Astro] [--km 5] [--min-regions 3] [--max-labels 25]
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ComplexHeatmap); library(circlize) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
CTS <- trimws(strsplit(getarg("--celltypes", "Astro"), ",")[[1]])
KM  <- as.integer(getarg("--km", "5"))
## Requiring the motif to be over-represented in EVERY region (not the published
## union) cuts Astro up 547 -> 102 and down 721 -> 41. With DA peak sets this large,
## FindMotifs calls nearly all 746 JASPAR motifs enriched somewhere, so the union
## is close to unselective.
MIN_REG  <- as.integer(getarg("--min-regions", "3"))
MAX_LAB  <- as.integer(getarg("--max-labels", "25"))
REGIONS <- c("INS", "PFC", "VST")
BASE <- file.path(OUT_ROOT, "04-ATAC_DE", "motif_chromvar")
OUT  <- file.path(OUT_ROOT, "04-ATAC_DE", "figures", "figS11")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
safe <- function(x) gsub("[^A-Za-z0-9_.+-]", "_", x)

for (ct in CTS) {
  ## per-cell chromVAR z-scores -> mean per region x condition
  cells <- lapply(REGIONS, function(R) {
    f <- file.path(BASE, R, paste0(safe(ct), "_chromvar_cells.rds"))
    if (!file.exists(f)) { message("[S11] missing ", f); return(NULL) }
    x <- readRDS(f)
    cond <- as.character(pwh_relabel(x$meta$condition_new))
    sapply(PWH_LEVELS, function(cl) rowMeans(x$z[, cond == cl, drop = FALSE]))
  })
  names(cells) <- REGIONS
  cells <- cells[!vapply(cells, is.null, logical(1))]
  if (!length(cells)) { message("[S11] no chromVAR for ", ct); next }
  motif_names <- utils::read.csv(file.path(BASE, names(cells)[1], "motif_names.csv"),
                                 stringsAsFactors = FALSE)
  nm <- setNames(motif_names$name, motif_names$motif)

  act <- do.call(cbind, lapply(names(cells), function(R) {
    m <- cells[[R]]; colnames(m) <- paste(R, colnames(m), sep = "-"); m }))
  ## published column order: region blocks, PWoH then PWH inside each
  col_orders <- as.vector(t(outer(names(cells), PWH_LEVELS, paste, sep = "-")))
  act <- act[, col_orders, drop = FALSE]

  ## donor-aware pooled result, for the added marker
  xrf <- file.path(BASE, "crossregion", paste0(safe(ct), "_chromvar_crossregion.csv"))
  xr <- if (file.exists(xrf)) utils::read.csv(xrf, stringsAsFactors = FALSE) else NULL

  for (dirn in c("up", "down")) {
    ## motifs over-represented in that direction's DA peaks in >= MIN_REG regions
    per_region <- lapply(names(cells), function(R) {
      f <- file.path(BASE, R, sprintf("%s_findmotifs_%s.csv", safe(ct), dirn))
      if (!file.exists(f)) return(character(0))
      d <- utils::read.csv(f, stringsAsFactors = FALSE)
      if (!nrow(d) || !"p.adjust" %in% names(d)) return(character(0))
      d$motif[d$p.adjust < 0.05]
    })
    tab <- table(unlist(per_region))
    sel <- intersect(names(tab)[tab >= MIN_REG], rownames(act))
    message(sprintf("[S11] %s %s: %d motifs over-represented in >= %d of %d regions (union would be %d)",
                    ct, dirn, length(sel), MIN_REG, length(per_region), length(tab)))
    if (length(sel) < 3) { message("[S11]   too few motifs - no heatmap"); next }

    mat <- t(scale(t(act[sel, , drop = FALSE])))          # z-score across the 6 columns
    rownames(mat) <- nm[sel]
    keep <- !apply(mat, 1, function(r) any(!is.finite(r)))
    mat <- mat[keep, , drop = FALSE]
    utils::write.csv(data.frame(motif = sel[keep], name = rownames(mat), mat, check.names = FALSE),
                     file.path(OUT, sprintf("figS11_%s_motif_activity_%s.csv", tolower(ct), dirn)),
                     row.names = FALSE)

    donor_sig <- if (is.null(xr)) rep(FALSE, nrow(mat)) else
      sel[keep] %in% xr$motif[xr$FDR < 0.05]
    ## Too many rows to name one by one: label the donor-aware significant motifs,
    ## topped up with the largest-swing ones, and hide the rest of the row names.
    rng <- apply(mat, 1, function(r) diff(range(r)))
    lab_idx <- which(donor_sig)
    if (length(lab_idx) > MAX_LAB) lab_idx <- lab_idx[order(-rng[lab_idx])][seq_len(MAX_LAB)]
    if (length(lab_idx) < MAX_LAB) {
      extra <- setdiff(order(-rng), lab_idx)
      lab_idx <- sort(c(lab_idx, head(extra, MAX_LAB - length(lab_idx))))
    }
    ## the mark labels are only needed when the row names are hidden, otherwise
    ## every motif would be named twice
    show_names <- nrow(mat) <= 45
    ra <- if (show_names)
      rowAnnotation(`donor-aware\nFDR < 0.05` =
                      anno_simple(ifelse(donor_sig, "yes", "no"),
                                  col = c(yes = "black", no = "grey92")),
                    annotation_name_gp = gpar(fontsize = 8))
    else
      rowAnnotation(`donor-aware\nFDR < 0.05` =
                      anno_simple(ifelse(donor_sig, "yes", "no"),
                                  col = c(yes = "black", no = "grey92")),
                    mark = anno_mark(at = lab_idx, labels = rownames(mat)[lab_idx],
                                     labels_gp = gpar(fontsize = 6.5)),
                    annotation_name_gp = gpar(fontsize = 8))
    ha <- HeatmapAnnotation(
      Region = sub("-.*", "", colnames(mat)),
      Condition = sub(".*-", "", colnames(mat)),
      col = list(Region = c(INS = "#E69F00", PFC = "#56B4E9", VST = "#CC79A7"),
                 Condition = setNames(c("#0072B2", "darkorange1"), PWH_LEVELS)),
      annotation_name_gp = gpar(fontsize = 8))
    p <- Heatmap(mat, name = "z-scored\nmotif activity",
                 col = colorRamp2(c(-2, 0, 2), c("steelblue", "white", "red1")),
                 top_annotation = ha, right_annotation = ra,
                 row_km = KM, cluster_columns = FALSE, show_column_names = TRUE,
                 show_row_names = show_names,
                 row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 9),
                 column_title = sprintf("%s: motifs enriched in peaks %s in PWH", ct, dirn))
    f <- file.path(OUT, sprintf("figS11_%s_motif_activity_%s.png", tolower(ct), dirn))
    h <- min(9.5, max(4.5, 0.055 * nrow(mat) + 2.2))   # fits a journal page
    png(f, width = 7, height = h, units = "in", res = 300); draw(p); dev.off()
    grDevices::cairo_pdf(sub("\\.png$", ".pdf", f), width = 7, height = h); draw(p); dev.off()
    message(sprintf("[S11]   %d motifs drawn (%d donor-aware significant) -> %s",
                    nrow(mat), sum(donor_sig), basename(f)))
  }
}
message("[S11] done -> ", OUT)
