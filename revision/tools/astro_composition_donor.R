# ==============================================================================
# Donor-level composition test for astrocyte subtypes: is each subtype's share of
# a donor's astrocytes different between PWH and PWoH? (function library, used by
# revision/tools/run_astro_reviewer_figures.R and astro_cluster_robustness.R)
#
# Method
#   Cell-level tests (FindMarkers over cells) treat every nucleus as independent,
#   so a donor contributing 4,000 cells counts 4,000 times. Here each subtype's
#   fraction of a donor's astrocytes is modelled at the donor level with the same
#   adjustment set as the main DE:
#       cbind(n_subtype, n_other) ~ HIV + Age + Sex + Race + PMI
#   quasibinomial GLM (overdispersion-aware), one observation per donor, including
#   donors with zero cells of that subtype. A subtype almost absent from controls
#   is the case where within-subtype DE is not answerable but the compositional
#   difference is. BH FDR across subtypes. Race levels held by a single donor
#   (Native Alaskan/American Indian) are collapsed to "Other".
#   Also reports per-donor counts per subtype and whether each subtype can support
#   a donor-level pseudobulk DE (>= MIN_DONORS donors per group with >= MIN_CELLS
#   cells); where it cannot, it is reported as not estimable.
#   (The regional edgeR DE is already donor-aware: pseudobulk by sample_id with
#   donor-level covariates.)
#
# Inputs   astrocyte Seurat object (in session); SAMPLE_META from revision/config/paths.R
# Outputs  OUT_ROOT/06-astrocyte_analysis/
#            astro_composition_donor_test<sfx>.csv, astro_within_subtype_feasibility<sfx>.csv,
#            astro_per_donor_counts<sfx>.csv, astro_composition_per_donor<sfx>.png/.pdf
#
# Usage
#   source(file.path(SCORCH_REPO, "revision/tools/astro_composition_donor.R"))
#   res <- astro_composition(clean_scorch_astro, label_col = "subcluster_labels")
#   res$results          # per-subtype effect, p, FDR
#   res$feasibility      # which subtypes can support within-subtype DE
#   res$plot             # per-donor proportions by condition
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
  library(Seurat); library(ggplot2); library(dplyr)
})

MIN_CELLS  <- 10   # cells of a subtype for a donor to count toward that subtype
MIN_DONORS <- 4    # donors per group needed before within-subtype DE is sensible

#' @param obj        astrocyte object carrying subcluster labels
#' @param label_col  metadata column with the subtype labels
#' @param donor_col  donor identifier (one row per donor in the model)
#' @param cond_col   HIV/CTR
#' @param outdir     NULL to skip writing
astro_composition <- function(obj,
                              label_col = "subcluster_labels",
                              donor_col = "sample_ids",
                              cond_col  = "condition_new",
                              outdir = file.path(OUT_ROOT, "06-astrocyte_analysis"),
                              sample_meta = SAMPLE_META) {

  md <- obj@meta.data
  for (cc in c(label_col, donor_col, cond_col))
    if (!cc %in% colnames(md)) stop("missing metadata column: ", cc)

  md[[label_col]] <- as.character(md[[label_col]])
  md[[donor_col]] <- as.character(md[[donor_col]])

  ## ---- per-donor x subtype counts -------------------------------------------
  tab <- as.data.frame.matrix(table(md[[donor_col]], md[[label_col]]))
  tab$total <- rowSums(tab)
  subtypes <- setdiff(colnames(tab), "total")
  message("[comp] ", nrow(tab), " donors, ", length(subtypes), " subtypes, ",
          sum(tab$total), " astrocytes")

  ## donor-level condition, taken from the cells (and checked for consistency)
  dcond <- tapply(as.character(md[[cond_col]]), md[[donor_col]],
                  function(v) { u <- unique(v); if (length(u) != 1)
                    stop("donor with >1 condition") ; u })
  dcond <- dcond[rownames(tab)]

  ## ---- donor-level covariates, same adjustment set as the main DE ----------
  meta <- utils::read.csv(sample_meta, stringsAsFactors = FALSE)
  meta <- meta[!duplicated(meta$sample_ids), ]
  rownames(meta) <- as.character(meta$sample_ids)
  keep <- rownames(tab)[rownames(tab) %in% rownames(meta)]
  if (length(keep) < nrow(tab))
    warning(nrow(tab) - length(keep), " donor(s) absent from the sample sheet: ",
            paste(setdiff(rownames(tab), keep), collapse = ", "))
  tab   <- tab[keep, , drop = FALSE]
  dcond <- dcond[keep]

  df <- data.frame(
    donor = keep,
    ## modelling keeps CTR/HIV as the factor levels so the design and the
    ## "HIVHIV" coefficient name stay valid; display uses PWH/PWoH via pwh_relabel()
    HIV   = factor(ifelse(dcond == "HIV", "HIV", "CTR"), levels = c("CTR","HIV")),
    Age   = as.numeric(meta[keep, "age_value_v2"]),
    Sex   = factor(make.names(meta[keep, "sex_v2"])),
    Race  = factor(make.names(meta[keep, "race_v2"])),
    PMI   = as.numeric(meta[keep, "postmortem_interval_v2"]),
    total = tab$total,
    stringsAsFactors = FALSE
  )
  ## As in the main DE: a race level held by one donor is not identified and can
  ## make the design rank-deficient, so the rare Native Alaskan/American Indian
  ## levels are collapsed.
  rare <- grepl("Native", as.character(df$Race))
  if (any(rare)) {
    df$Race <- as.character(df$Race); df$Race[rare] <- "Other"; df$Race <- factor(df$Race)
    message("[comp] collapsed ", sum(rare), " donor(s) into Race='Other' (single-donor levels are not identified)")
  }
  if (anyNA(df$Age) || anyNA(df$PMI))
    stop("NA in Age/PMI after the sample-sheet join; the model would drop those donors")
  message("[comp] donors: ", sum(df$HIV == "HIV"), " PWH / ", sum(df$HIV == "CTR"), " PWoH")

  ## ---- feasibility of a within-subtype pseudobulk DE ------------------------
  feas <- do.call(rbind, lapply(subtypes, function(s) {
    n <- tab[[s]]
    data.frame(subtype = s,
               cells_total = sum(n),
               donors_PWH_ok  = sum(n[df$HIV == "HIV"] >= MIN_CELLS),
               donors_PWoH_ok = sum(n[df$HIV == "CTR"] >= MIN_CELLS),
               stringsAsFactors = FALSE)
  }))
  feas$within_subtype_DE <- ifelse(
    feas$donors_PWH_ok >= MIN_DONORS & feas$donors_PWoH_ok >= MIN_DONORS,
    "feasible", "NOT ESTIMABLE - report composition only")
  feas <- feas[order(-feas$cells_total), ]

  ## ---- composition test, per subtype ---------------------------------------
  res <- do.call(rbind, lapply(subtypes, function(s) {
    y <- cbind(tab[[s]], df$total - tab[[s]])
    fit <- try(glm(y ~ HIV + Age + Sex + Race + PMI, data = df,
                   family = quasibinomial()), silent = TRUE)
    if (inherits(fit, "try-error"))
      return(data.frame(subtype = s, est = NA, se = NA, p = NA,
                        note = "model failed", stringsAsFactors = FALSE))
    co <- summary(fit)$coefficients
    if (!"HIVHIV" %in% rownames(co))
      return(data.frame(subtype = s, est = NA, se = NA, p = NA,
                        note = "HIV term dropped (rank deficient)", stringsAsFactors = FALSE))
    data.frame(subtype = s,
               est = unname(co["HIVHIV","Estimate"]),
               se  = unname(co["HIVHIV","Std. Error"]),
               p   = unname(co["HIVHIV",4]),
               note = "", stringsAsFactors = FALSE)
  }))
  res$odds_ratio <- round(exp(res$est), 3)
  res$FDR <- p.adjust(res$p, method = "BH")
  ## descriptive means alongside the model, so a significant p with a trivial
  ## difference is visible rather than implied
  res$mean_pct_PWoH <- round(100 * sapply(res$subtype, function(s)
    mean(tab[df$HIV == "CTR", s] / df$total[df$HIV == "CTR"])), 2)
  res$mean_pct_PWH <- round(100 * sapply(res$subtype, function(s)
    mean(tab[df$HIV == "HIV", s] / df$total[df$HIV == "HIV"])), 2)
  res <- res[order(res$p), c("subtype","mean_pct_PWoH","mean_pct_PWH",
                             "odds_ratio","est","se","p","FDR","note")]

  ## ---- per-donor proportions, plotted --------------------------------------
  long <- do.call(rbind, lapply(subtypes, function(s)
    data.frame(donor = df$donor, HIV = df$HIV, subtype = s,
               pct = 100 * tab[[s]] / df$total, stringsAsFactors = FALSE)))
  long$HIV <- factor(as.character(pwh_relabel(long$HIV)), levels = PWH_LEVELS)
  ## publication styling: shared theme, no subtitle, subtypes in the Fig1D tSNE order
  if (!exists("theme_astro_pub"))
    source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
  long$subtype <- factor(long$subtype, levels = names(astro_subtype_palette(long$subtype)))
  p <- ggplot(long, aes(x = HIV, y = pct, fill = HIV)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6, linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.9, alpha = 0.85, stroke = 0) +
    facet_wrap(~ subtype, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = PWH_COLOURS) +
    labs(y = "Donor's astrocytes (%)", x = NULL,
         title = "Astrocyte subtype abundance per donor") +
    theme_astro_pub() + theme(legend.position = "none")

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write.csv(res,  file.path(outdir, .sfxname("astro_composition_donor_test.csv")), row.names = FALSE)
    write.csv(feas, file.path(outdir, .sfxname("astro_within_subtype_feasibility.csv")), row.names = FALSE)
    write.csv(cbind(donor = rownames(tab), HIV = as.character(df$HIV), tab),
              file.path(outdir, .sfxname("astro_per_donor_counts.csv")), row.names = FALSE)
    fp <- file.path(outdir, .sfxname("astro_composition_per_donor.png"))
    save_pub(p, fp, width = 7.2, height = 4.2)
    astro_caption(fp, paste(
      "One point per donor, the unit of analysis; boxes show median and interquartile range.",
      "The donor-level composition test (quasibinomial, adjusted) is in",
      "astro_composition_donor_test.csv."))
    message("[comp] wrote 3 CSVs + PNG/PDF to ", outdir)
  }

  message("\n[comp] composition test (donor-level, quasibinomial, adjusted):")
  print(res, row.names = FALSE, digits = 3)
  message("\n[comp] can a within-subtype DE be run at all? (>= ", MIN_DONORS,
          " donors/group with >= ", MIN_CELLS, " cells)")
  print(feas, row.names = FALSE)

  invisible(list(results = res, feasibility = feas, counts = tab,
                 donor_data = df, plot = p))
}
