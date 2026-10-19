# ==============================================================================
# Per-gene, per-cell-type ambient RNA fraction, and functions that flag regional DE
# calls with it (function library; the driver is revision/tools/run_ambient_fraction.R).
#
# Method
#   Primary DE runs on uncorrected counts, because CellBender over-corrects the
#   reactive-astrocyte panel and creates artifactual complete separation. This
#   library quantifies, for each reported call, how much of the gene's signal was
#   ambient. It removes no genes and changes no test; it attaches a number to each
#   call so the caveat can be stated quantitatively.
#
#   For gene g in cell type k:
#       ambient_frac(g,k) = 1 - sum(corrected[g, cells in k]) / sum(RAW[g, cells in k])
#   i.e. the fraction of the observed signal that CellBender attributed to ambient
#   RNA, over exactly the cells that formed that cell type's pseudobulk, so it joins
#   onto the edgeR results one-to-one on gene x cluster_id. It is CellBender's
#   estimate, not ground truth: a high value is a reason for caution, not proof
#   that the biology is wrong. Pairs with raw_sum < 10 are set to NA.
#
#   Companion columns:
#   ambient_frac_donor_median  same quantity within each donor, then the median
#                              across donors. Far below ambient_frac means one
#                              sample drives it.
#   ambient_frac_other         the same fraction in all other cell types pooled.
#                              Separates ambient signal present everywhere (high
#                              in both) from ambient signal specific to this cell
#                              type (high here, low elsewhere); the latter means
#                              the gene's presence in this cell type is the part
#                              the model called contamination.
#
#   Rounding: the corrected DE run rounds CellBender's fractional posterior counts
#   to integers for the QL F-test. This uses the unrounded posterior, the better
#   estimate of how much signal was removed; the two differ by <<1% at the
#   pseudobulk sum level.
#
# Inputs   OUT_ROOT/01-Regional_data_processing/<REGION>/data/.../integrated_seurat_obj_v3.RDS
#          (corrected counts) plus the uncorrected counts via uncorrected_counts()
#          from revision/config/paths.R
# Outputs  OUT_ROOT/03-DE_analysis/<REGION>/results/<REGION>_ambient_fraction_by_celltype.csv
#
# Usage
#   source("revision/config/paths.R"); source("revision/tools/ambient_fraction.R")
#   amb <- ambient_fraction_table("INS")        # slow: loads the region's objects
#   de  <- readRDS(".../INS_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
#   de  <- annotate_de_with_ambient(de, amb, threshold = 0.10)
#   ambient_flag_summary(de)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr)
})

AMBIENT_CT_COL     <- "consensus_annotation_v3"
AMBIENT_SAMPLE_COL <- "sample_ids"

#' Path of the cached table for a region.
ambient_table_path <- function(region) {
  region <- toupper(region)
  file.path(OUT_ROOT, "03-DE_analysis", region, "results",
            paste0(region, "_ambient_fraction_by_celltype.csv"))
}

#' Regional object used by the DE scripts, so the cells match exactly.
#'
#' These paths are copied from the readRDS() lines of 03-DE_analysis/
#' 01/02/03-*_edgeR.R and must stay in step with them. VST reads from a
#' different stage directory than INS and PFC - using the wrong object would
#' give a cell set that does not correspond to the DE, silently.
.ambient_region_object <- function(region) {
  region <- toupper(region)
  sub <- switch(region,
                INS = "03-cell_integration_vis",
                PFC = "03-cell_integration_vis",
                VST = "04-VST-tsne_marker_plotting",
                stop("unknown region '", region, "' - add its object path here"))
  file.path(OUT_ROOT, "01-Regional_data_processing", region, "data", sub,
            "integrated_seurat_obj_v3.RDS")
}

# ------------------------------------------------------------------------------
#' Build the gene x cell-type ambient fraction table for one region.
#'
#' @param region "INS", "PFC" or "VST"
#' @param obj    optional, an already-loaded object (skips the read)
#' @param write  cache to ambient_table_path(region)
#'
#' @return data.frame: region, cluster_id, gene, n_cells, n_donors, raw_sum,
#'         corr_sum, ambient_frac, ambient_frac_donor_median, ambient_frac_other
ambient_fraction_table <- function(region, obj = NULL, write = TRUE,
                                   min_raw_sum = 10) {
  region <- toupper(region)

  if (is.null(obj)) {
    p <- .ambient_region_object(region)
    if (!file.exists(p)) stop("missing regional object: ", p)
    message("[amb] loading ", p)
    obj <- readRDS(p)
  }
  stopifnot(AMBIENT_CT_COL %in% colnames(obj@meta.data),
            AMBIENT_SAMPLE_COL %in% colnames(obj@meta.data))

  ## ---- the two matrices, same cells, same genes, same order ----------------
  corr <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  colnames(corr) <- colnames(obj); rownames(corr) <- rownames(obj)
  message("[amb] corrected: ", nrow(corr), " x ", ncol(corr),
          if (any(corr@x != floor(corr@x))) "  (fractional posterior, kept unrounded)" else "")

  raw <- uncorrected_counts(obj, region)
  stopifnot(identical(dimnames(raw), dimnames(corr)))

  ## Genes that uncorrected_counts() had to zero-fill (absent from a sample's
  ## RAW assay) would look like 100% ambient, which is meaningless. A gene with
  ## corrected signal but zero RAW signal is exactly that case; frac() returns NA.
  ct  <- as.character(obj@meta.data[[AMBIENT_CT_COL]])
  sid <- as.character(obj@meta.data[[AMBIENT_SAMPLE_COL]])
  cts <- sort(unique(ct[!is.na(ct)]))
  message("[amb] ", length(cts), " cell types, ", length(unique(sid)), " donors")

  ## global sums, for the "all other cell types" column
  raw_tot  <- Matrix::rowSums(raw)
  corr_tot <- Matrix::rowSums(corr)

  frac <- function(r, c) {
    f <- 1 - (c / r)
    f[!is.finite(f) | r <= 0] <- NA_real_
    pmin(pmax(f, 0), 1)                     # clamp to [0, 1]
  }

  out <- lapply(cts, function(k) {
    idx <- which(ct == k)
    r <- Matrix::rowSums(raw[, idx, drop = FALSE])
    c <- Matrix::rowSums(corr[, idx, drop = FALSE])

    ## per-donor, then median across donors that actually have this cell type
    sk <- sid[idx]
    dn <- unique(sk)
    dm <- vapply(rownames(raw), function(g) NA_real_, numeric(1))
    if (length(dn) > 1L) {
      per <- vapply(dn, function(d) {
        j <- idx[sk == d]
        rr <- Matrix::rowSums(raw[, j, drop = FALSE])
        cc <- Matrix::rowSums(corr[, j, drop = FALSE])
        ff <- frac(rr, cc)
        ff[rr < min_raw_sum] <- NA_real_     # unstable in a tiny donor x ct cell
        ff
      }, numeric(nrow(raw)))
      dm <- apply(per, 1, function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE))
    }

    data.frame(region = region, cluster_id = k, gene = rownames(raw),
               n_cells = length(idx), n_donors = length(dn),
               raw_sum = as.numeric(r), corr_sum = as.numeric(c),
               ambient_frac = frac(r, c),
               ambient_frac_donor_median = as.numeric(dm),
               ambient_frac_other = frac(raw_tot - r, corr_tot - c),
               stringsAsFactors = FALSE, row.names = NULL)
  })

  res <- dplyr::bind_rows(out)
  ## Below min_raw_sum the ratio is noise, not an estimate, so it is set to NA.
  n_thin <- sum(res$raw_sum < min_raw_sum)
  res$ambient_frac[res$raw_sum < min_raw_sum] <- NA_real_
  message("[amb] ", nrow(res), " gene x cell-type pairs; ", n_thin,
          " (", round(100 * n_thin / nrow(res), 1), "%) below raw_sum ",
          min_raw_sum, " -> NA")
  message("[amb] median ambient_frac ",
          round(stats::median(res$ambient_frac, na.rm = TRUE), 3),
          "; >10% in ", round(100 * mean(res$ambient_frac > 0.10, na.rm = TRUE), 1), "% of pairs")

  if (write) {
    f <- ambient_table_path(region)
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, f, row.names = FALSE)
    message("[amb] wrote ", f)
  }
  res
}

# ------------------------------------------------------------------------------
#' Join the ambient fraction onto a DE result and flag calls above a threshold.
#'
#' Adds: ambient_frac, ambient_frac_donor_median, ambient_frac_other,
#'       ambient_flag (logical), ambient_specific (high here, low elsewhere).
#'
#' Nothing is removed. `ambient_flag` is a REPORTING annotation.
annotate_de_with_ambient <- function(de, amb, threshold = 0.10,
                                     specific_ratio = 2) {
  de <- as.data.frame(de)
  stopifnot(all(c("gene", "cluster_id") %in% colnames(de)))
  a <- amb[, c("cluster_id", "gene", "ambient_frac",
               "ambient_frac_donor_median", "ambient_frac_other")]
  a$cluster_id <- as.character(a$cluster_id)
  de$cluster_id <- as.character(de$cluster_id)

  n0 <- nrow(de)
  de <- dplyr::left_join(de, a, by = c("cluster_id", "gene"))
  stopifnot(nrow(de) == n0)

  miss <- sum(is.na(de$ambient_frac))
  if (miss)
    message("[amb] ", miss, " of ", n0, " DE rows have no ambient estimate ",
            "(gene x cell type below the raw-count floor) - flag is NA, not FALSE")

  de$ambient_flag <- de$ambient_frac > threshold
  ## "soup specifically here": high in this cell type, comparatively low elsewhere
  de$ambient_specific <- de$ambient_flag &
    (de$ambient_frac > specific_ratio * pmax(de$ambient_frac_other, 1e-6))

  message("[amb] flagged (>", threshold, "): ",
          sum(de$ambient_flag, na.rm = TRUE), " of ", n0, " tested rows")
  de
}

# ------------------------------------------------------------------------------
#' Flag rate among tested vs among significant rows, at several thresholds. If
#' the rate is the same in both, ambient content is not what drives the
#' significant calls.
ambient_flag_summary <- function(de, thresholds = c(0.05, 0.10, 0.25, 0.50),
                                 fdr = 0.05) {
  stopifnot("ambient_frac" %in% colnames(de))
  sig <- !is.na(de$p_adj.glb) & de$p_adj.glb < fdr
  do.call(rbind, lapply(thresholds, function(t) {
    f <- de$ambient_frac > t
    data.frame(threshold = t,
               pct_tested = 100 * mean(f, na.rm = TRUE),
               pct_significant = 100 * mean(f[sig], na.rm = TRUE),
               n_significant_flagged = sum(f[sig], na.rm = TRUE),
               n_significant = sum(sig),
               enrichment = mean(f[sig], na.rm = TRUE) / mean(f, na.rm = TRUE),
               row.names = NULL)
  }))
}

#' Per-cell-type version of the same table, at one threshold.
ambient_flag_by_celltype <- function(de, threshold = 0.10, fdr = 0.05) {
  de %>%
    dplyr::mutate(sig = !is.na(p_adj.glb) & p_adj.glb < fdr,
                  flag = ambient_frac > threshold) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n_tested = dplyr::n(),
                     pct_tested_flagged = 100 * mean(flag, na.rm = TRUE),
                     n_sig = sum(sig),
                     n_sig_flagged = sum(flag & sig, na.rm = TRUE),
                     pct_sig_flagged = 100 * mean(flag[sig], na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(pct_sig_flagged))
}

#' The flagged significant genes, highest ambient fraction first (supplementary
#' table).
ambient_flagged_genes <- function(de, threshold = 0.10, fdr = 0.05, n = Inf) {
  d <- de[!is.na(de$ambient_frac) & de$ambient_frac > threshold &
            !is.na(de$p_adj.glb) & de$p_adj.glb < fdr, , drop = FALSE]
  keep <- intersect(c("cluster_id", "gene", "logFC", "logCPM", "p_adj.glb",
                      "avg_CTR_frq", "avg_HIV_frq", "ambient_frac",
                      "ambient_frac_donor_median", "ambient_frac_other",
                      "ambient_specific"), colnames(d))
  d <- d[order(-d$ambient_frac), keep, drop = FALSE]
  utils::head(d, n)
}
