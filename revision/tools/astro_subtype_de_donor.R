# Donor-aware marker DE for an astrocyte subtype: paired pseudobulk, within donor.
#
# Method: replaces the cell-level Wilcoxon FindMarkers test (subtype vs other
#   astrocytes), whose p-values are anticonservative because the effective n is
#   the number of donors, not nuclei (reviewer request). Each donor is pseudobulked
#   twice, over its cells of the target subtype and over its remaining astrocytes
#   (or over a named comparator subtype), and edgeR QL tests
#       ~ donor + half
#   so donor absorbs between-person variation. Because the comparison is internal
#   to each donor no control arm is needed: PWoH donors contribute almost no
#   Reactive-PWH cells, and that scarcity is the HIV-enrichment result (see
#   revision/tools/astro_composition_donor.R). A donor enters only if both halves
#   have >= 10 cells; the test needs >= 5 such donors. Markers: FDR < 0.05 and
#   |logFC| > 0.25. HIV vs PWoH within the subtype is not tested: with few PWoH
#   donors contributing cells that contrast is not estimable.
#
# Input:  a Seurat astrocyte object with subtype labels (passed by the caller).
# Output: OUT_ROOT/06-astrocyte_analysis/astro_<subtype>[_vs_<comparator>]_donorDE.csv,
#         _donor_counts.csv and _donor_heatmap.png.
#
# Usage: defines astro_subtype_de(); sourced by revision/tools/run_astro_reviewer_figures.R
#   and other astrocyte figure scripts, e.g.
#     res <- astro_subtype_de(clean_scorch_astro, subtype = "Reactive-PWH")

## Optional output-name suffix (e.g. "_no1040vst"), set by
## run_astro_reviewer_figures.R --suffix, so each file names the sample set it
## came from. Empty by default.
if (!exists("ASTRO_OUT_SUFFIX")) ASTRO_OUT_SUFFIX <- ""
.sfxname <- function(name) {
  if (!nzchar(ASTRO_OUT_SUFFIX)) return(name)
  sub("(\\.[A-Za-z0-9]+)$", paste0(ASTRO_OUT_SUFFIX, "\\1"), name)
}

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(edgeR); library(limma)
  library(ggplot2); library(dplyr)
})

MIN_CELLS_PER_HALF <- 10   # cells required in BOTH halves, so no donor contributes a noise-dominated profile
MIN_DONORS_PAIRED  <- 5    # below this the paired test is not worth running

#' @param obj        astrocyte object with subtype labels
#' @param subtype    the target subtype, e.g. "Reactive-PWH"
#' @param label_col  metadata column holding subtype labels
#' @param donor_col  donor id
#' @param restrict_to_condition  optionally "HIV" to run within PWH donors only
astro_subtype_de <- function(obj,
                             subtype,
                             ## vs_subtype: NULL keeps the published contrast
                             ## (subtype vs every other astrocyte). Naming a subtype
                             ## (e.g. the reactive population, present in both PWoH
                             ## and PWH) restricts the comparator half to it, which
                             ## asks what distinguishes this state from ordinary
                             ## reactivity rather than from astrocytes in general.
                             ## Cells in neither half are dropped.
                             vs_subtype = NULL,
                             label_col = "subcluster_labels",
                             donor_col = "sample_ids",
                             cond_col  = "condition_new",
                             restrict_to_condition = NULL,
                             outdir = file.path(OUT_ROOT, "06-astrocyte_analysis"),
                             fdr = 0.05, min_lfc = 0.25, n_heat = 30,
                             assay = "RNA") {
  ## assay: lets the same paired test run on snATAC gene-activity counts
  ## (Fig S14b). The default "RNA" is the Fig 1E analysis.

  md <- obj@meta.data
  for (cc in c(label_col, donor_col, cond_col))
    if (!cc %in% colnames(md)) stop("missing metadata column: ", cc)
  lab <- as.character(md[[label_col]])
  if (!subtype %in% lab) stop("subtype '", subtype, "' not present in ", label_col)

  don <- as.character(md[[donor_col]])
  cnd <- as.character(md[[cond_col]])
  if (!is.null(restrict_to_condition)) {
    keep <- cnd == restrict_to_condition
    message("[sde] restricted to ", restrict_to_condition, ": ", sum(keep), " of ", length(keep), " cells")
    obj <- obj[, keep]; md <- obj@meta.data
    lab <- lab[keep]; don <- don[keep]; cnd <- cnd[keep]
  }

  ## ---- which donors can contribute a matched pair --------------------------
  if (!is.null(vs_subtype)) {
    miss <- setdiff(vs_subtype, lab)
    if (length(miss)) stop("comparator subtype not present in ", label_col, ": ",
                           paste(miss, collapse = ", "))
    message("[sde] comparator half restricted to: ", paste(vs_subtype, collapse = ", "))
  }
  in_other <- if (is.null(vs_subtype)) lab != subtype else lab %in% vs_subtype
  n_in  <- table(don[lab == subtype])
  n_out <- table(don[in_other])
  all_d <- union(names(n_in), names(n_out))
  cnt <- data.frame(
    donor    = all_d,
    n_subtype = as.integer(n_in[all_d]),
    n_other   = as.integer(n_out[all_d]),
    stringsAsFactors = FALSE
  )
  cnt$n_subtype[is.na(cnt$n_subtype)] <- 0L
  cnt$n_other[is.na(cnt$n_other)]     <- 0L
  cnt$condition <- vapply(cnt$donor, function(d) unique(cnd[don == d])[1], character(1))
  cnt$condition_label <- as.character(pwh_relabel(cnt$condition))
  cnt$usable <- cnt$n_subtype >= MIN_CELLS_PER_HALF & cnt$n_other >= MIN_CELLS_PER_HALF
  cnt <- cnt[order(-cnt$n_subtype), ]

  use <- cnt$donor[cnt$usable]
  message("[sde] donors contributing a matched pair (>= ", MIN_CELLS_PER_HALF,
          " cells in BOTH halves): ", length(use), " of ", nrow(cnt))
  message("[sde]   by condition: ",
          paste(sprintf("%s=%d", names(table(cnt$condition[cnt$usable])),
                        table(cnt$condition[cnt$usable])), collapse = "  "))
  if (length(use) < MIN_DONORS_PAIRED)
    stop("only ", length(use), " donors can contribute a pair; the paired test is ",
         "not worth running. Report the composition result instead.")

  ## ---- pseudobulk each donor twice -----------------------------------------
  work <- obj
  if (length(Layers(work[[assay]], search = "counts")) > 1) work[[assay]] <- JoinLayers(work[[assay]])
  cts <- LayerData(work, assay = assay, layer = "counts")
  grp <- ifelse(lab == subtype, "subtype", ifelse(in_other, "other", NA_character_))

  keycols <- paste(don, grp, sep = "__")
  wanted  <- as.vector(outer(use, c("subtype","other"), paste, sep = "__"))
  idx <- keycols %in% wanted
  cts <- cts[, idx, drop = FALSE]; keycols <- keycols[idx]

  ## sum counts within each donor x half
  f <- factor(keycols, levels = wanted)
  pb <- as.matrix(cts %*% Matrix::sparse.model.matrix(~ 0 + f))
  colnames(pb) <- levels(f)
  message("[sde] pseudobulk: ", nrow(pb), " genes x ", ncol(pb),
          " profiles (", length(use), " donors x 2)")

  meta <- data.frame(
    sample  = colnames(pb),
    donor   = sub("__.*$", "", colnames(pb)),
    half    = factor(sub("^.*__", "", colnames(pb)), levels = c("other","subtype")),
    stringsAsFactors = FALSE
  )
  stopifnot(identical(meta$sample, colnames(pb)))

  ## ---- paired edgeR --------------------------------------------------------
  ## donor as a blocking factor: every donor contributes both halves, so the
  ## subtype coefficient is estimated WITHIN donor and between-person variation
  ## is absorbed rather than inflating the comparison.
  y <- DGEList(counts = pb, samples = meta)
  design <- model.matrix(~ donor + half, data = meta)
  keep <- filterByExpr(y, design)
  message("[sde] genes kept by filterByExpr: ", sum(keep), " of ", nrow(y))
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  coef_name <- grep("^halfsubtype$", colnames(design), value = TRUE)
  stopifnot(length(coef_name) == 1)
  res <- glmQLFTest(fit, coef = coef_name)
  tt <- topTags(res, n = Inf)$table
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene","logFC","logCPM","F","PValue","FDR")]

  n_sig <- sum(tt$FDR < fdr & abs(tt$logFC) > min_lfc, na.rm = TRUE)
  message("[sde] markers at FDR < ", fdr, " & |logFC| > ", min_lfc, ": ", n_sig)
  message("[sde] BCV: ", round(sqrt(y$common.dispersion), 3),
          "  (donor-level biological variation)")

  ## ---- per-donor heatmap of the top markers --------------------------------
  ## Expression per person (reviewer request), so a marker driven
  ## by one or two donors is visible rather than hidden behind a p-value.
  top <- head(tt$gene[tt$FDR < fdr & tt$logFC > min_lfc], n_heat)
  hm <- NULL
  if (length(top) >= 3) {
    cpm_mat <- edgeR::cpm(y, log = TRUE)[top, , drop = FALSE]
    long <- data.frame(
      gene  = rep(rownames(cpm_mat), ncol(cpm_mat)),
      sample = rep(colnames(cpm_mat), each = nrow(cpm_mat)),
      value = as.vector(cpm_mat), stringsAsFactors = FALSE)
    long$donor <- sub("__.*$", "", long$sample)
    long$half  <- sub("^.*__", "", long$sample)
    long$condition <- pwh_relabel(cnt$condition[match(long$donor, cnt$donor)])
    long$gene  <- factor(long$gene, levels = rev(top))
    long$half  <- factor(long$half, levels = c("other","subtype"))
    hm <- ggplot(long, aes(x = donor, y = gene, fill = value)) +
      geom_tile() +
      facet_grid(~ half + condition, scales = "free_x", space = "free_x") +
      scale_fill_viridis_c(name = "log2 CPM") +
      labs(x = NULL, y = NULL,
           title = paste0(subtype, " markers, per donor"),
           subtitle = "each column is one person; paired halves of the same donor side by side") +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            strip.text = element_text(size = 7))
  }

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    sfx <- gsub("[^A-Za-z0-9]+", "_", subtype)
    ## a restricted comparator writes to its own files, so the published
    ## subtype-vs-rest tables are never overwritten
    if (!is.null(vs_subtype))
      sfx <- paste0(sfx, "_vs_", gsub("[^A-Za-z0-9]+", "_", paste(vs_subtype, collapse = "_")))
    write.csv(tt,  file.path(outdir, .sfxname(paste0("astro_", sfx, "_donorDE.csv"))), row.names = FALSE)
    write.csv(cnt, file.path(outdir, .sfxname(paste0("astro_", sfx, "_donor_counts.csv"))), row.names = FALSE)
    if (!is.null(hm))
      ggsave(file.path(outdir, .sfxname(paste0("astro_", sfx, "_donor_heatmap.png"))), hm,
             width = 2 + 0.28 * length(use) * 2, height = 2 + 0.18 * length(top),
             dpi = 200, limitsize = FALSE)
    message("[sde] wrote results to ", outdir)
  }

  message("\n[sde] top markers:")
  print(utils::head(tt, 15), row.names = FALSE, digits = 3)
  message("\n[sde] donor contribution (usable = both halves >= ", MIN_CELLS_PER_HALF, "):")
  print(utils::head(cnt, 15), row.names = FALSE)

  invisible(list(de = tt, donors = cnt, dge = y, design = design, heatmap = hm))
}
