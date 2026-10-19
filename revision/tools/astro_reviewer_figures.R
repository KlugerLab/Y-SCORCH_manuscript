# ==============================================================================
# Figures answering the reviewer's astrocyte requests: per-donor contribution to
# each subtype, donor-level abundance, robustness to clustering parameters, and
# QC/doublet metrics (function library; driver revision/tools/run_astro_reviewer_figures.R).
#
# Method
#   ask 1  per-donor contribution      -> astro_fig_donor_contribution()
#   ask 2  donor-level abundance       -> astro_composition() in
#                                         revision/tools/astro_composition_donor.R
#   ask 3  robustness to resolution    -> astro_fig_resolution_stability() (the
#                                         multi-setting analysis is
#                                         revision/tools/astro_cluster_robustness.R)
#   ask 4  QC / doublet metrics        -> astro_fig_qc_metrics()
#   The question is whether a subtype is a real state or a few donors' cells, so
#   every figure has donor on an axis or as a point. Also defines the shared
#   astrocyte subtype palette (ASTRO_SUBTYPE_COLOURS) used by Fig1D and the
#   supplementals.
#
# Inputs   labelled astrocyte object (in session); scDblFinder calls from
#          OUT_ROOT/01-Regional_data_processing/<REGION>/data/sample_metadata_withMT/
# Outputs  OUT_ROOT/06-astrocyte_analysis/figures/
#            ask1A_subtype_composition_per_donor<sfx>.png, ask1B_donor_origin_per_subtype<sfx>.png,
#            ask1_donor_dominance<sfx>.csv, ask3_resolution_stability_<subtype><sfx>.png/.csv,
#            ask4_qc_doublet_by_subtype<sfx>.png, ask4_qc_medians<sfx>.csv (+ .pdf, _caption.txt)
#
# Usage: source(file.path(SCORCH_REPO, "revision/tools/astro_reviewer_figures.R"))
# ==============================================================================

## Output-name suffix, so a file carries the sample set it came from in its own
##   name. Set by run_astro_reviewer_figures.R --suffix (e.g. "_no1040vst" when
##   donor 1040's VST astrocytes are excluded). Empty by default.
if (!exists("ASTRO_OUT_SUFFIX")) ASTRO_OUT_SUFFIX <- ""
.sfxname <- function(name) {
  if (!nzchar(ASTRO_OUT_SUFFIX)) return(name)
  sub("(\\.[A-Za-z0-9]+)$", paste0(ASTRO_OUT_SUFFIX, "\\1"), name)
}

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(ggplot2); library(dplyr)
})

# ------------------------------------------------------------------------------
# ASK 1 — per-donor contribution to each cluster.
#
# Two complementary views:
#   A. per donor: what fraction of this person's astrocytes is each subtype?
#      -> is the subtype present across people, or only in a few?
#   B. per subtype: what fraction of this subtype's cells came from each donor?
#      -> is the subtype built out of one donor's cells?
# B answers whether a subtype is driven by a subset of cases.
# ------------------------------------------------------------------------------
## Subtype colours ("pastel_blue") = the manuscript astrocyte tSNE (Fig1D), used
##   for Fig1D (tSNE + UMAP), the reviewer figures and the labelled supplementals.
##   Based on ggplot's default hue palette (scales::hue_pal()(7)) in the level
##   order below, except Signaling -> soft blue hcl(250, 100, 65) and Micro-Astro
##   -> teal hcl(190, 100, 65), because the default Signaling teal #00C094 sits
##   too close to Reactive-PWH green #53B400 where the two mix in the UMAP. Same
##   hcl family as the rest, so the palette stays pastel.
ASTRO_SUBTYPE_COLOURS <- setNames(
  c("#F8766D", "#C49A00", "#53B400", "#3FA1FF", "#00C0BC", "#A58AFF", "#FB61D7"),
  c("Homeostatic", "Reactive", "Reactive-PWH", "Signaling", "Micro-Astro", "Synaptic", "Mural-Astro"))
## the two contamination clusters, shown only in the 9-group supplementals: greys,
## so they read as "not an astrocyte state"
ASTRO_EXTRA_COLOURS <- c(Oligo = "#9E9E9E", `Endo-Astro` = "#5E5E5E")

astro_fig_donor_contribution <- function(obj,
                                         label_col = "subcluster_labels",
                                         donor_col = "sample_ids",
                                         cond_col  = "condition_new",
                                         outdir = file.path(OUT_ROOT, "06-astrocyte_analysis/figures")) {
  md <- obj@meta.data
  d <- data.frame(donor = as.character(md[[donor_col]]),
                  lab   = as.character(md[[label_col]]),
                  cond  = as.character(md[[cond_col]]), stringsAsFactors = FALSE)
  d$cond <- pwh_relabel(d$cond)   # PWH / PWoH for display

  ## subtype order and colours follow the tSNE; anything unexpected goes last, grey
  present   <- unique(d$lab)
  lab_order <- c(intersect(names(ASTRO_SUBTYPE_COLOURS), present),
                 sort(setdiff(present, names(ASTRO_SUBTYPE_COLOURS))))
  sub_cols  <- setNames(ifelse(lab_order %in% names(ASTRO_SUBTYPE_COLOURS),
                               ASTRO_SUBTYPE_COLOURS[lab_order], "grey60"), lab_order)
  if (any(sub_cols == "grey60"))
    warning("subtypes with no tSNE colour: ", paste(lab_order[sub_cols == "grey60"], collapse = ", "))
  d$lab <- factor(d$lab, levels = lab_order)
  ## formatC, not format(): format() pads a vector to a common width ("n =    357")
  fmt <- function(x) formatC(as.integer(x), format = "d", big.mark = ",")

  ## donor order: condition, then size, so PWH/PWoH blocks read cleanly
  dtot <- d %>% count(donor, cond, name = "n_total")
  dtot <- dtot[order(dtot$cond, -dtot$n_total), ]
  d$donor <- factor(d$donor, levels = dtot$donor)
  dtot$donor <- factor(dtot$donor, levels = dtot$donor)

  ## n per subtype (legend) and per condition (facet strips / legend)
  n_lab   <- table(d$lab)
  lab_key <- sprintf("%s (n = %s)", lab_order, fmt(n_lab[lab_order]))
  csum    <- dtot %>% group_by(cond) %>%
    summarise(n_donors = dplyr::n(), n_cells = sum(n_total), .groups = "drop")
  cond_key <- setNames(sprintf("%s (%d donors, n = %s astrocytes)", csum$cond,
                               csum$n_donors, fmt(csum$n_cells)), as.character(csum$cond))

  ## --- A: composition within each donor ---
  a <- d %>% count(donor, cond, lab, name = "n") %>%
    group_by(donor) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
  ## publication styling: shared theme, no subtitle - the reading instructions
  ## go to the caption file written at the end
  if (!exists("theme_astro_pub"))
    source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
  pA <- ggplot(a, aes(x = donor, y = pct, fill = lab)) +
    geom_col(width = 0.85) +
    ## astrocytes per donor, above each bar
    geom_text(data = dtot, aes(x = donor, y = 101.5, label = fmt(n_total)),
              inherit.aes = FALSE, angle = 90, hjust = 0, size = 2.1) +
    ## two-line strip: the PWoH facet is narrower than a one-line label and clipped it
    facet_grid(~ cond, scales = "free_x", space = "free_x",
               labeller = as_labeller(setNames(
                 sprintf("%s\n%d donors, %s nuclei", csum$cond, csum$n_donors, fmt(csum$n_cells)),
                 as.character(csum$cond)))) +
    scale_fill_manual(values = sub_cols, breaks = lab_order, labels = lab_key) +
    scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 116),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Donor", y = "Donor's astrocytes (%)", fill = "Subtype",
         title = "Astrocyte subtype composition per donor") +
    theme_astro_pub() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
          legend.position = "right")

  ## --- B: where each subtype's cells came from ---
  b <- d %>% count(lab, donor, cond, name = "n") %>%
    group_by(lab) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
  bn <- b %>% group_by(lab) %>%
    summarise(n_cells = sum(n), n_donors = sum(n > 0), .groups = "drop")
  ## per-subtype donor dominance: share of the subtype's nuclei from its top donor(s)
  dom <- b %>% group_by(lab) %>%
    summarise(n_cells = sum(n),
              n_donors = sum(n > 0),
              n_donors_ge10 = sum(n >= 10),
              top_donor = donor[which.max(pct)],
              top_donor_pct = round(max(pct), 1),
              top3_pct = round(sum(sort(pct, decreasing = TRUE)[1:3], na.rm = TRUE), 1),
              .groups = "drop") %>%
    arrange(-top_donor_pct)

  ## Fill stays PWoH/PWH - that is what this panel shows. The subtype's tSNE
  ## colour is a strip under each bar, and n cells / n donors sits above it.
  pB <- ggplot(b, aes(x = lab, y = pct, fill = cond)) +
    geom_col(aes(group = donor), colour = "white", linewidth = 0.12, width = 0.8) +
    geom_tile(data = bn, aes(x = lab, y = -3.5), inherit.aes = FALSE,
              height = 4, width = 0.8, fill = unname(sub_cols[as.character(bn$lab)])) +
    geom_text(data = bn, inherit.aes = FALSE, vjust = 0, size = 2.4, lineheight = 0.9,
              aes(x = lab, y = 101.5,
                  label = sprintf("n = %s\n%d donors", fmt(n_cells), n_donors))) +
    scale_fill_manual(values = PWH_COLOURS, breaks = names(cond_key),
                      labels = unname(cond_key), name = NULL) +
    scale_y_continuous(breaks = seq(0, 100, 25), limits = c(-6, 114),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "Subtype's nuclei (%)",
         title = "Donor origin of each astrocyte subtype") +
    theme_astro_pub() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "bottom", legend.direction = "vertical")

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    fA <- file.path(outdir, .sfxname("ask1A_subtype_composition_per_donor.png"))
    fB <- file.path(outdir, .sfxname("ask1B_donor_origin_per_subtype.png"))
    save_pub(pA, fA, width = 7.2, height = 3.8)
    save_pub(pB, fB, width = 4.8, height = 4.2)
    astro_caption(fA, paste(
      "Each bar is one donor; the number above it is that donor's astrocyte nuclei.",
      "Legend n = nuclei per subtype; facet labels give donors and nuclei per condition.",
      "Seven astrocyte subtypes (Oligo and Endo-Astro contamination clusters excluded);",
      "colours match the astrocyte subtype tSNE (Fig. 1D)."))
    astro_caption(fB, paste(
      "Each segment is one donor, filled by condition; a subtype dominated by a few",
      "segments would be donor-driven. Above each bar: nuclei and contributing donors.",
      "The strip under each bar is the subtype colour from the tSNE (Fig. 1D)."))
    write.csv(dom, file.path(outdir, .sfxname("ask1_donor_dominance.csv")), row.names = FALSE)
  }
  message("[fig] donor dominance per subtype (top_donor_pct is the key number):")
  print(as.data.frame(dom), row.names = FALSE)
  invisible(list(per_donor = pA, per_subtype = pB, dominance = dom))
}

# ------------------------------------------------------------------------------
# ASK 4 — QC and doublet metrics per cluster.
#
# scDblFinder was run per sample upstream and its calls are in the stage-01
# metadata CSVs, but the columns did not survive into the combined object, so
# they are joined back here by barcode. Every retained cell is a called singlet
# (the doublets were removed upstream), which is shown explicitly - the score
# still varies and an elevated score in one cluster is informative even when no
# cell was called.
# ------------------------------------------------------------------------------
astro_join_doublet_scores <- function(obj, donor_col = "sample_ids") {
  md <- obj@meta.data
  if ("scDblFinder.score" %in% colnames(md)) {
    message("[qc] scDblFinder.score already present"); return(obj)
  }
  if (!"region" %in% colnames(md)) stop("need a 'region' column to locate the metadata CSVs")
  cells <- colnames(obj)
  sids  <- as.character(md[[donor_col]])
  regs  <- toupper(as.character(md$region))
  ## cell names are "<barcode>-<sample_ids>[-<REGION>]"
  suf_long  <- paste0("-", sids, "-", regs); suf_short <- paste0("-", sids)
  ew <- function(x, s) substr(x, nchar(x) - nchar(s) + 1L, nchar(x)) == s
  ul <- ew(cells, suf_long); us <- !ul & ew(cells, suf_short)
  bcs <- substr(cells, 1L, nchar(cells) - ifelse(ul, nchar(suf_long), nchar(suf_short)))

  score <- rep(NA_real_, length(cells)); klass <- rep(NA_character_, length(cells))
  for (rg in unique(regs)) {
    dir <- file.path(OUT_ROOT, "01-Regional_data_processing", rg, "data", "sample_metadata_withMT")
    if (!dir.exists(dir)) { warning("no metadata dir for ", rg); next }
    for (sid in unique(sids[regs == rg])) {
      hit <- list.files(dir, pattern = paste0("^", sid, "_.*\\.csv$"), full.names = TRUE)
      if (length(hit) != 1) next
      m <- utils::read.csv(hit[1], row.names = 1)
      if (!"scDblFinder.score" %in% colnames(m)) next
      i <- which(regs == rg & sids == sid)
      j <- match(bcs[i], rownames(m))
      score[i] <- m$scDblFinder.score[j]
      klass[i] <- as.character(m$scDblFinder.class[j])
    }
  }
  message("[qc] joined scDblFinder for ", sum(!is.na(score)), " of ", length(cells), " cells")
  obj$scDblFinder.score <- score
  obj$scDblFinder.class <- klass
  obj
}

astro_fig_qc_metrics <- function(obj,
                                 label_col = "subcluster_labels",
                                 outdir = file.path(OUT_ROOT, "06-astrocyte_analysis/figures")) {
  if (!"scDblFinder.score" %in% colnames(obj@meta.data))
    obj <- astro_join_doublet_scores(obj)
  md <- obj@meta.data
  vars <- intersect(c("nCount_RNA","nFeature_RNA","percent.mt","scDblFinder.score"),
                    colnames(md))
  if (!length(vars)) stop("no QC columns found")
  long <- do.call(rbind, lapply(vars, function(v)
    data.frame(lab = as.character(md[[label_col]]), metric = v,
               value = as.numeric(md[[v]]), stringsAsFactors = FALSE)))
  long <- long[!is.na(long$value), ]

  ## medians on the RAW values, before any display transform
  tab <- long %>% group_by(lab, metric) %>%
    summarise(median = round(median(value), 3), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = metric, values_from = median)

  ## publication styling: shared theme, no subtitle, one row, tSNE subtype colours
  ## and order, readable metric names. UMI and gene counts are drawn on log10 -
  ## on a linear scale the long right tail flattens every violin.
  if (!exists("theme_astro_pub"))
    source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
  sub_cols <- astro_subtype_palette(long$lab)
  metric_names <- c(nCount_RNA = "UMIs per nucleus (log10)",
                    nFeature_RNA = "Genes per nucleus (log10)",
                    percent.mt = "Mitochondrial reads (%)",
                    scDblFinder.score = "scDblFinder doublet score")
  pl <- long
  logm <- pl$metric %in% c("nCount_RNA", "nFeature_RNA")
  pl$value[logm] <- log10(pl$value[logm])
  pl$lab    <- factor(pl$lab, levels = names(sub_cols))
  pl$metric <- factor(metric_names[pl$metric], levels = metric_names[vars])

  p <- ggplot(pl, aes(x = lab, y = value, fill = lab)) +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey25") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = sub_cols) +
    labs(x = NULL, y = NULL, title = "QC and doublet metrics by astrocyte subtype") +
    theme_astro_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    f4 <- file.path(outdir, .sfxname("ask4_qc_doublet_by_subtype.png"))
    save_pub(p, f4, width = 2.4 * length(vars) + 0.4, height = 2.9)
    astro_caption(f4, paste(
      "All retained nuclei are scDblFinder singlets (doublets were removed upstream);",
      "the continuous doublet score is shown because it still varies between subtypes.",
      "UMI and gene counts on a log10 scale; boxes show median and interquartile range.",
      "Medians per subtype are in ask4_qc_medians.csv."))
    write.csv(tab, file.path(outdir, .sfxname("ask4_qc_medians.csv")), row.names = FALSE)
  }
  message("[fig] median QC per subtype:"); print(as.data.frame(tab), row.names = FALSE)
  invisible(list(plot = p, table = tab, obj = obj))
}

# ------------------------------------------------------------------------------
# ASK 3 — robustness to clustering / parameter choices.
#
# Re-cluster at a range of resolutions and ask whether the cells currently
# labelled as the target subtype stay together. If the population is real it
# should survive as a recognisable cluster across parameters; if it is an
# artifact of one resolution it will dissolve.
#
# Reported per resolution:
#   n_clusters        how finely the data was cut
#   best_jaccard      overlap between the target set and its best-matching new cluster
#   pct_in_best       what fraction of the target cells landed in that one cluster
# A population that keeps pct_in_best high across resolutions is robust.
#
# This re-clusters, so it is slow. It does not modify the caller's object.
# ------------------------------------------------------------------------------
astro_fig_resolution_stability <- function(obj,
                                           subtype,
                                           label_col = "subcluster_labels",
                                           reduction = "harmony_astro",
                                           dims = 1:30,
                                           resolutions = c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0),
                                           outdir = file.path(OUT_ROOT, "06-astrocyte_analysis/figures")) {
  target <- colnames(obj)[as.character(obj@meta.data[[label_col]]) == subtype]
  if (!length(target)) stop("no cells labelled ", subtype)
  message("[stab] target: ", length(target), " ", subtype, " cells")
  if (!reduction %in% names(obj@reductions))
    stop("reduction '", reduction, "' not on the object; have: ",
         paste(names(obj@reductions), collapse = ", "))

  work <- FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE)
  out <- do.call(rbind, lapply(resolutions, function(r) {
    w <- FindClusters(work, resolution = r, verbose = FALSE)
    cl <- as.character(Idents(w))
    names(cl) <- colnames(w)
    tc <- cl[target]
    tab <- table(tc)
    best <- names(tab)[which.max(tab)]
    in_best <- as.integer(max(tab))
    size_best <- sum(cl == best)
    data.frame(resolution = r,
               n_clusters = length(unique(cl)),
               best_cluster = best,
               pct_in_best = round(100 * in_best / length(target), 1),
               best_jaccard = round(in_best / (length(target) + size_best - in_best), 3),
               n_clusters_spanned = sum(tab >= 0.05 * length(target)),
               stringsAsFactors = FALSE)
  }))

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    sfx <- gsub("[^A-Za-z0-9]+", "_", subtype)
    write.csv(out, file.path(outdir, .sfxname(paste0("ask3_resolution_stability_", sfx, ".csv"))),
              row.names = FALSE)
  }
  p <- astro_plot_resolution_stability(out, subtype, outdir = outdir)
  message("[stab] stability across resolutions:"); print(out, row.names = FALSE)
  invisible(list(plot = p, table = out))
}

## Drawing is separate from the re-clustering so the figure can be redrawn from
## the saved table: run_astro_reviewer_figures.R --only ask3plot.
## Cluster counts sit under the resolution on the x axis, clear of the Jaccard line.
astro_plot_resolution_stability <- function(out, subtype,
                                            outdir = file.path(OUT_ROOT, "06-astrocyte_analysis/figures")) {
  if (!exists("theme_astro_pub"))
    source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
  out <- out[order(out$resolution), ]
  out$res_lab <- factor(sprintf("%s\n(%d)", format(out$resolution), out$n_clusters),
                        levels = sprintf("%s\n(%d)", format(out$resolution), out$n_clusters))
  BAR <- "#0072B2"; LINE <- "darkorange1"
  p <- ggplot(out, aes(x = res_lab)) +
    geom_col(aes(y = pct_in_best), fill = BAR, width = 0.65) +
    geom_line(aes(y = best_jaccard * 100, group = 1), colour = LINE, linewidth = 0.6) +
    geom_point(aes(y = best_jaccard * 100), colour = LINE, size = 1.4) +
    scale_y_continuous(name = paste0(subtype, " nuclei in one cluster (%)"),
                       limits = c(0, 100), expand = expansion(mult = c(0, 0.02)),
                       sec.axis = sec_axis(~ . / 100, name = "Jaccard index with best cluster")) +
    labs(x = "Clustering resolution (number of clusters)",
         title = paste0("Robustness of ", subtype, " to clustering resolution")) +
    theme_astro_pub() +
    theme(axis.title.y.left  = element_text(colour = BAR),
          axis.title.y.right = element_text(colour = LINE))

  if (!is.null(outdir)) {
    sfx <- gsub("[^A-Za-z0-9]+", "_", subtype)
    f3 <- file.path(outdir, .sfxname(paste0("ask3_resolution_stability_", sfx, ".png")))
    save_pub(p, f3, width = 4.6, height = 3.3)
    astro_caption(f3, paste0(
      "The astrocytes were re-clustered at each resolution. Bars (left axis): percentage of ",
      subtype, " nuclei that fall in a single cluster. Line (right axis): Jaccard overlap ",
      "between the ", subtype, " nuclei and their best-matching cluster. Numbers in ",
      "parentheses: clusters at that resolution."))
  }
  p
}
