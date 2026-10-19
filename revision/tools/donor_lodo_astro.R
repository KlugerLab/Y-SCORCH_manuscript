# Leave-one-donor-out (LODO) sensitivity of the astrocyte subtype composition result (SERPINA3/OSMR reactive program).
#
# Method: reviewer request to show the major glial findings hold when outlier or
#   pathology cases are excluded, including leave-one-donor-out analyses (the
#   microglial findings are in revision/tools/donor_lodo_microglia.R). The claim
#   tested is composition: the Reactive-PWH subtype is 1.2% of astrocytes in PWoH
#   and 12.9% in PWH (OR 12.0, FDR 0.023).
#   Per subtype, a quasibinomial GLM (subtype count vs other astrocytes)
#   ~ HIV + Age + Sex + Race + PMI on a donor x subtype count table, as in
#   astro_composition(); BH across subtypes. The table comes from metadata alone,
#   so each LODO iteration is a row drop plus a refit. Fits are run on the full
#   cohort, with each donor dropped in turn, and with pathology donors 1039
#   (HIVE + PML) and 1181 (multi-pathology) dropped together.
#   - Removing a donor can make the GLM rank-deficient (mainly a Race level held
#     by one person); rare Race levels are collapsed as in astro_composition(),
#     and fits that lose the HIV term are reported, not silently dropped.
#   - Percentages are computed per donor and then averaged, not pooled over
#     cells, so a donor with many astrocytes cannot dominate them.
#   Oligo and Endo-Astro contamination clusters are excluded.
#
# Inputs:  OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled_no1040vst.RDS
#          (or --in <file>; falls back, with a warning, to astro_all_clusters_labeled.RDS)
#          SAMPLE_META (donor covariates)
# Outputs: OUT_ROOT/03-DE_analysis/donor_sensitivity/
#          astro_composition_lodo_allfits[_no1040vst].csv, astro_composition_lodo_summary[_no1040vst].csv
#
# Usage: revision/docker/run.sh revision/tools/donor_lodo_astro.R [--in <astrocyte object file>]

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

suppressPackageStartupMessages({library(Seurat)})

LODO_OUT <- file.path(OUT_ROOT, "03-DE_analysis", "donor_sensitivity")
dir.create(LODO_OUT, recursive = TRUE, showWarnings = FALSE)

PATHOLOGY_DONORS <- c("1039", "1181")   # HIVE+PML, and multi-pathology
FDR_THRESHOLD    <- 0.05
DROP_LABELS      <- c("Oligo", "Endo-Astro")   # contamination clusters, as in 06-astrocyte_analysis/astrocyte_analysis.Rmd

## --in selects the object; otherwise the object without donor 1040's VST
##   astrocytes is preferred. Falling back to the unfiltered object is announced,
##   not silent, because the LODO numbers must match the figures they accompany.
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }

.astro_dir <- file.path(OUT_ROOT, "06-astrocyte_analysis", "data")
.cands <- c(
  if (nzchar(getarg("--in", ""))) file.path(.astro_dir, getarg("--in", "")) else NA_character_,
  file.path(.astro_dir, "astro_all_clusters_labeled_no1040vst.RDS"),
  file.path(.astro_dir, "astro_all_clusters_labeled.RDS")
)
IN <- .cands[!is.na(.cands) & file.exists(.cands)][1]
if (is.na(IN)) stop("no astrocyte object found in ", .astro_dir)
if (!grepl("no1040", IN))
  warning("[lodo-astro] using ", basename(IN), " - still contains donor 1040's ",
          "VST astrocytes", immediate. = TRUE)

## output suffix follows the input, so LODO results cannot be confused with the
## pre-exclusion ones
SFX <- if (grepl("no1040", IN)) "_no1040vst" else ""

message("[lodo-astro] loading ", IN, " ...")
obj <- readRDS(IN)
LAB <- "subcluster_labels"
stopifnot(LAB %in% colnames(obj@meta.data))

md <- obj@meta.data
md[[LAB]] <- as.character(md[[LAB]])
md$sample_ids <- as.character(md$sample_ids)
md <- md[!(md[[LAB]] %in% DROP_LABELS), , drop = FALSE]
message("[lodo-astro] ", nrow(md), " astrocytes after dropping ",
        paste(DROP_LABELS, collapse = "/"))

## ---- donor x subtype table, built once ---------------------------------------
tab_full <- as.data.frame.matrix(table(md$sample_ids, md[[LAB]]))
subtypes <- colnames(tab_full)

dcond <- tapply(as.character(md$condition_new), md$sample_ids, function(v) {
  u <- unique(v); if (length(u) != 1) stop("donor with >1 condition"); u })

sm <- utils::read.csv(SAMPLE_META, stringsAsFactors = FALSE)
sm <- sm[!duplicated(sm$sample_ids), ]
rownames(sm) <- as.character(sm$sample_ids)

donors_all <- rownames(tab_full)[rownames(tab_full) %in% rownames(sm)]
message("[lodo-astro] ", length(donors_all), " donors, ", length(subtypes), " subtypes")

## ---- one composition fit on a donor subset -----------------------------------
.fit_comp <- function(drop = character(0)) {
  keep <- setdiff(donors_all, drop)
  if (length(keep) < 8) return(NULL)
  tab <- tab_full[keep, , drop = FALSE]
  total <- rowSums(tab)

  df <- data.frame(
    donor = keep,
    HIV   = factor(ifelse(dcond[keep] == "HIV", "HIV", "CTR"), levels = c("CTR","HIV")),
    Age   = as.numeric(sm[keep, "age_value_v2"]),
    Sex   = factor(make.names(sm[keep, "sex_v2"])),
    Race  = factor(make.names(sm[keep, "race_v2"])),
    PMI   = as.numeric(sm[keep, "postmortem_interval_v2"]),
    total = total, stringsAsFactors = FALSE)

  ## same rare-level collapse as astro_composition(): a Race level held by one
  ## donor is not identified and makes the design rank-deficient
  rare <- grepl("Native", as.character(df$Race))
  if (any(rare)) {
    df$Race <- as.character(df$Race); df$Race[rare] <- "Other"; df$Race <- factor(df$Race)
  }
  if (anyNA(df$Age) || anyNA(df$PMI)) return(NULL)
  if (nlevels(droplevels(df$HIV)) < 2) return(NULL)

  do.call(rbind, lapply(subtypes, function(s) {
    y <- cbind(tab[[s]], df$total - tab[[s]])
    fit <- try(stats::glm(y ~ HIV + Age + Sex + Race + PMI, data = df,
                          family = stats::quasibinomial()), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    co <- summary(fit)$coefficients
    if (!"HIVHIV" %in% rownames(co)) return(NULL)
    ## descriptive percentages: per donor, then averaged by group
    pct <- 100 * tab[[s]] / df$total
    data.frame(subtype = s,
               est = unname(co["HIVHIV","Estimate"]),
               p   = unname(co["HIVHIV",4]),
               odds_ratio   = unname(exp(co["HIVHIV","Estimate"])),
               mean_pct_CTR = mean(pct[df$HIV == "CTR"]),
               mean_pct_HIV = mean(pct[df$HIV == "HIV"]),
               n_donors = nrow(df), stringsAsFactors = FALSE)
  }))
}

.add_fdr <- function(x) { if (is.null(x)) return(NULL); x$FDR <- p.adjust(x$p, "BH"); x }

## ---- full cohort, LODO, and the targeted exclusion ---------------------------
full <- .add_fdr(.fit_comp())
if (is.null(full)) stop("full-cohort composition fit failed")

fits <- lapply(donors_all, function(d) {
  f <- .add_fdr(.fit_comp(drop = d))
  if (!is.null(f)) f$dropped <- d
  f
})
ok <- !vapply(fits, is.null, logical(1))
if (any(!ok))
  message("[lodo-astro] ", sum(!ok), " refit(s) failed (rank deficiency or ",
          "too few donors): ", paste(donors_all[!ok], collapse = ", "))
lodo <- do.call(rbind, fits[ok])

present <- intersect(PATHOLOGY_DONORS, donors_all)
excl <- if (length(present)) .add_fdr(.fit_comp(drop = present)) else NULL
message("[lodo-astro] pathology donors present in the astrocyte object: ",
        if (length(present)) paste(present, collapse = ", ") else "NONE")

## ---- collapse ----------------------------------------------------------------
res <- do.call(rbind, lapply(subtypes, function(s) {
  fu <- full[full$subtype == s, ]
  lo <- lodo[lodo$subtype == s, ]
  if (!nrow(fu) || !nrow(lo)) return(NULL)
  ex <- if (!is.null(excl)) excl[excl$subtype == s, ] else NULL
  worst <- which.max(lo$FDR)
  data.frame(
    subtype          = s,
    n_donors_full    = fu$n_donors,
    mean_pct_CTR     = round(fu$mean_pct_CTR, 2),
    mean_pct_HIV     = round(fu$mean_pct_HIV, 2),
    odds_ratio_full  = round(fu$odds_ratio, 3),
    fdr_full         = fu$FDR,
    lodo_fdr_min     = min(lo$FDR),
    lodo_fdr_max     = max(lo$FDR),
    lodo_fdr_median  = stats::median(lo$FDR),
    lodo_or_min      = round(min(lo$odds_ratio), 3),
    lodo_or_max      = round(max(lo$odds_ratio), 3),
    lodo_n_fits      = nrow(lo),
    lodo_n_sig       = sum(lo$FDR < FDR_THRESHOLD),
    lodo_all_sig     = all(lo$FDR < FDR_THRESHOLD),
    ## does the effect ever reverse? OR crossing 1 means the direction flipped
    direction_stable = all(sign(lo$est) == sign(fu$est)),
    most_influential = lo$dropped[worst],
    fdr_excl_pathology = if (!is.null(ex) && nrow(ex)) ex$FDR else NA_real_,
    or_excl_pathology  = if (!is.null(ex) && nrow(ex)) round(ex$odds_ratio, 3) else NA_real_,
    excluded_donors  = paste(present, collapse = ","),
    stringsAsFactors = FALSE)
}))
res <- res[order(res$fdr_full), ]

utils::write.csv(lodo, file.path(LODO_OUT, paste0("astro_composition_lodo_allfits", SFX, ".csv")),
                 row.names = FALSE)
utils::write.csv(res,  file.path(LODO_OUT, paste0("astro_composition_lodo_summary", SFX, ".csv")),
                 row.names = FALSE)
message("\n[lodo-astro] wrote ",
        file.path(LODO_OUT, paste0("astro_composition_lodo_summary", SFX, ".csv")))

## ---- console summary ---------------------------------------------------------
message("\n########## LODO: astrocyte composition ##########")
pr <- data.frame(
  subtype    = substr(res$subtype, 1, 16),
  pct_PWoH   = res$mean_pct_CTR,
  pct_PWH    = res$mean_pct_HIV,
  OR         = res$odds_ratio_full,
  fdr_full   = signif(res$fdr_full, 3),
  lodo_worst = signif(res$lodo_fdr_max, 3),
  OR_range   = paste0(res$lodo_or_min, "-", res$lodo_or_max),
  n_sig      = paste0(res$lodo_n_sig, "/", res$lodo_n_fits),
  robust     = ifelse(res$lodo_all_sig, "YES", "no"),
  excl_2d    = signif(res$fdr_excl_pathology, 3),
  stringsAsFactors = FALSE)
print(pr, row.names = FALSE)

sig <- res[res$fdr_full < FDR_THRESHOLD, ]
if (nrow(sig)) {
  nr <- sig[!sig$lodo_all_sig, ]
  if (nrow(nr)) {
    message("\nNOT robust to every single-donor removal:")
    for (i in seq_len(nrow(nr)))
      message("  ", nr$subtype[i], ": worst FDR ", signif(nr$lodo_fdr_max[i], 3),
              " when dropping donor ", nr$most_influential[i])
  } else {
    message("\nEvery full-cohort finding survives removal of any single donor.")
  }
  if (any(!sig$direction_stable))
    message("\nWARNING: effect direction reversed in some refits - see allfits csv")
}
message("\n[lodo-astro] done")
