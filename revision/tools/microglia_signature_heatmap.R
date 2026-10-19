# ==============================================================================
# Heatmaps of the microglial state gene lists across regions, donor-aware.
#
# Two panels:
#   A  gene x region, filled by the regional pseudobulk edgeR logFC: does the
#      gene move the same way in every region? Donor-aware because the DE is
#      pseudobulk (n = donors, never cells), but it does not show the donors.
#   B  gene x donor, filled by z-scored logCPM (within region), split by region
#      and condition: is the signature carried by the whole group or by a few
#      donors? Every column is one person, so a few-donor effect is visible.
#
# Colour: logFC and the z-score are diverging, so the ramp has two hues with a
# neutral grey midpoint. The poles are the paper's condition colours
# (PWH_COLOURS: blue #0072B2 = PWoH, orange = PWH; Okabe-Ito, colourblind-safe),
# so orange means higher in PWH as everywhere else in the paper.
#
# Genes in several signatures are repeated in each block so each list reads as
# a list; non-significant cells are not blanked (the asterisk carries
# significance, the fill carries magnitude); genes absent from a region's DE are
# drawn as explicit grey NA cells.
#
# Inputs: OUT_ROOT/03-DE_analysis/<R>/results/<R>_HIV-CTR_res_edgeR_ambientFlagged.csv;
#         microglia pseudobulk from revision/tools/microglia_signature_score.R
# Output: OUT_ROOT/03-DE_analysis/microglia_signatures/figures/
#         microglia_signatures_by_region.{png,pdf}, microglia_signatures_by_donor.{png,pdf}
#
# Usage: revision/docker/run.sh revision/tools/microglia_signature_heatmap.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/microglia_signature_score.R"))
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(edgeR) })

FIGDIR  <- file.path(MICRO_OUT, "figures")
REGIONS <- c("INS", "PFC", "VST")

## Diverging ramp: blue <- neutral grey -> orange. Grey, not white, so a
## zero-effect cell is visibly a measured zero rather than blank paper.
DIV_LOW  <- unname(PWH_COLOURS["PWoH"])   # #0072B2
DIV_MID  <- "#E8E6E3"
DIV_HIGH <- "darkorange1"
NA_FILL  <- "#BFBFBF"

# ------------------------------------------------------------------------------
#' Per-gene microglia DE across regions, read by COLUMN NAME.
#'
#' The `<donor>.frq` columns differ in number between regions, so column
#' positions differ too - indexing by number silently reads a different field in
#' another region (e.g. VST's F statistic read as a fold change of 33, which
#' looks like complete separation). Always look up by name.
microglia_de_long <- function(regions = REGIONS, cell_type = MICRO_LABEL) {
  do.call(rbind, lapply(regions, function(r) {
    f <- file.path(OUT_ROOT, "03-DE_analysis", r, "results",
                   paste0(r, "_HIV-CTR_res_edgeR_ambientFlagged.csv"))
    if (!file.exists(f)) { warning("no DE table for ", r); return(NULL) }
    d <- data.table::fread(f, showProgress = FALSE)
    need <- c("gene", "cluster_id", "logFC", "p_adj.loc")
    if (!all(need %in% names(d))) stop("missing columns in ", f)
    d <- d[d$cluster_id == cell_type, ]
    data.frame(region = r, gene = d$gene, logFC = d$logFC,
               p_adj = d$`p_adj.loc`,
               ambient = if ("ambient_frac" %in% names(d)) d$ambient_frac else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

#' Signature membership as a long table; genes repeat across signatures.
.sig_long <- function(sigs = MICROGLIA_SIGNATURES, include_exploratory = FALSE) {
  keep <- names(sigs)[include_exploratory |
                      !vapply(sigs, function(s) isTRUE(s$exploratory), logical(1))]
  do.call(rbind, lapply(keep, function(n)
    data.frame(signature = n, gene = sigs[[n]]$genes,
               direction = sigs[[n]]$direction, stringsAsFactors = FALSE)))
}

# ------------------------------------------------------------------------------
#' PANEL A - gene x region, filled by logFC.
microglia_heatmap_region <- function(de = NULL, sigs = MICROGLIA_SIGNATURES,
                                     cap = 6, include_exploratory = FALSE) {
  if (is.null(de)) de <- microglia_de_long()
  sl <- .sig_long(sigs, include_exploratory)

  ## Cross-join signature x gene against ALL regions BEFORE joining the DE, so a
  ## gene missing from one region's DE still gets a tile and is drawn in the NA
  ## colour. A plain left join drops those rows and leaves white gaps in the
  ## grid, which a reader parses as "no change" when it means "not measured".
  d <- merge(sl, data.frame(region = REGIONS, stringsAsFactors = FALSE), by = NULL)
  d <- merge(d, de, by = c("gene", "region"), all.x = TRUE)
  d$region <- factor(d$region, levels = REGIONS)

  ## order genes within each signature by mean effect, so blocks read top-down
  ord <- d %>% group_by(signature, gene) %>%
    summarise(m = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    arrange(signature, m)
  d$row <- factor(paste(d$signature, d$gene), levels = paste(ord$signature, ord$gene))
  d$star <- ifelse(!is.na(d$p_adj) & d$p_adj < 0.05, "*", "")
  d$fill <- pmax(pmin(d$logFC, cap), -cap)

  ggplot(d, aes(x = region, y = row, fill = fill)) +
    geom_tile(colour = "white", linewidth = 0.6) +          # 2px-equivalent gap
    geom_text(aes(label = star), size = 4, vjust = 0.78, colour = "grey15") +
    facet_grid(signature ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH,
                         midpoint = 0, limits = c(-cap, cap),
                         na.value = NA_FILL,
                         name = expression(log[2]~FC),
                         breaks = c(-cap, 0, cap),
                         labels = c(paste0("≤ -", cap), "0", paste0("≥ ", cap))) +
    scale_y_discrete(labels = function(x) sub("^[^ ]+ ", "", x)) +
    labs(title = "Microglial state signatures across regions",
         subtitle = paste0("pseudobulk edgeR, n = donors (28 / 29 / 27)\n",
                           "* p_adj < 0.05   |   orange = higher in PWH\n",
                           "genes repeat across signatures; grey = not measured"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
          strip.placement = "outside",
          strip.background = element_rect(fill = "grey93", colour = NA),
          panel.grid = element_blank(),
          axis.text.y = element_text(size = 6.5),
          legend.position = "right")
}

# ------------------------------------------------------------------------------
#' PANEL B - gene x donor, z-scored logCPM, split by region and condition.
#'
#' @param which_sigs signatures to draw. Defaults to the cross-region hits plus
#'   DAM_Haney2024 as a deliberate NEGATIVE contrast - showing a flat block
#'   beside the moving ones is what stops a reader concluding that everything
#'   moved.
microglia_heatmap_donor <- function(which_sigs = c("LDAM_Haney2024",
                                                   "Homeostatic",
                                                   "LDAM_Marschallinger2020_down",
                                                   "DAM_Haney2024"),
                                    sigs = MICROGLIA_SIGNATURES,
                                    regions = REGIONS, cap = 2.5) {
  sl <- .sig_long(sigs[which_sigs], include_exploratory = TRUE)

  long <- do.call(rbind, lapply(regions, function(r) {
    pb <- microglia_pseudobulk(r)                 # cached after the first build
    y  <- calcNormFactors(DGEList(counts = pb$counts))
    lc <- edgeR::cpm(y, log = TRUE, prior.count = 1)
    g  <- intersect(unique(sl$gene), rownames(lc))
    z  <- t(scale(t(lc[g, , drop = FALSE])))      # per gene, WITHIN region
    z  <- z[is.finite(rowSums(z)), , drop = FALSE]
    data.frame(region = r,
               gene   = rep(rownames(z), ncol(z)),
               donor  = rep(colnames(z), each = nrow(z)),
               z      = as.vector(z),
               condition = pb$meta$condition[match(rep(colnames(z), each = nrow(z)),
                                                   pb$meta$donor)],
               stringsAsFactors = FALSE)
  }))

  d <- merge(sl, long, by = "gene")
  d$region    <- factor(d$region, levels = regions)
  d$condition <- factor(as.character(d$condition), levels = c("PWoH", "PWH"))
  d$signature <- factor(d$signature, levels = which_sigs)
  d$row  <- factor(paste(d$signature, d$gene))
  d$fill <- pmax(pmin(d$z, cap), -cap)

  ## donors ordered within region x condition, so blocks are visually contiguous
  du <- unique(d[, c("region", "condition", "donor")])
  du <- du[order(du$region, du$condition, du$donor), ]
  d$donor <- factor(d$donor, levels = unique(du$donor))

  ggplot(d, aes(x = donor, y = row, fill = fill)) +
    geom_tile() +
    facet_grid(signature ~ region + condition,
               scales = "free", space = "free", switch = "y") +
    scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH,
                         midpoint = 0, limits = c(-cap, cap), na.value = NA_FILL,
                         name = "z (within region)",
                         breaks = c(-cap, 0, cap),
                         labels = c(paste0("≤ -", cap), "0", paste0("≥ ", cap))) +
    scale_y_discrete(labels = function(x) sub("^[^ ]+ ", "", x)) +
    labs(title = "Donor-level expression of the signature genes",
         subtitle = paste("every column is one person; z within region, so",
                          "regions are comparable on CONTRAST, not baseline.",
                          "\nDAM_Haney2024 is included as a negative contrast."),
         x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
          strip.text.x = element_text(size = 6.5),
          strip.placement = "outside",
          strip.background = element_rect(fill = "grey93", colour = NA),
          panel.grid = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6.5))
}

# ------------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

  message("[mgfig] panel A: gene x region")
  de <- microglia_de_long()
  message("[mgfig] DE rows for ", MICRO_LABEL, ": ",
          paste(sprintf("%s=%d", REGIONS,
                        vapply(REGIONS, function(r) sum(de$region == r), integer(1))),
                collapse = "  "))
  pA <- microglia_heatmap_region(de)
  ggsave(file.path(FIGDIR, "microglia_signatures_by_region.png"), pA,
         width = 6.2, height = 13, dpi = 300, limitsize = FALSE)
  ggsave(file.path(FIGDIR, "microglia_signatures_by_region.pdf"), pA,
         width = 6.2, height = 13, limitsize = FALSE)
  message("[mgfig] wrote panel A")

  message("[mgfig] panel B: gene x donor (builds/loads pseudobulk)")
  pB <- microglia_heatmap_donor()
  ggsave(file.path(FIGDIR, "microglia_signatures_by_donor.png"), pB,
         width = 13, height = 6.5, dpi = 300, limitsize = FALSE)
  ggsave(file.path(FIGDIR, "microglia_signatures_by_donor.pdf"), pB,
         width = 13, height = 6.5, limitsize = FALSE)
  message("[mgfig] wrote panel B")
  message("[mgfig] done - ", FIGDIR)
}
