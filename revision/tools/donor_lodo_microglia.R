# Leave-one-donor-out (LODO) sensitivity of the microglial state signatures, per region.
#
# Method: reviewer request for leave-one-donor-out analyses or, at minimum,
#   exclusion of the HIVE donor and the donor with multiple reported pathologies,
#   since most viral-positive nuclei come from those two. This covers microglial
#   activation (DAM / IRM / Homeostatic) and the lipid module (LDAM / MG4); the
#   astrocyte program is in revision/tools/donor_lodo_astro.R.
#   microglia_pseudobulk() (revision/tools/microglia_signature_score.R) caches a
#   donor x gene matrix per region (donors with >= 10 microglia), so each LODO
#   iteration is a column subset plus a camera refit; no Seurat object is loaded.
#   Reported per region x signature:
#   - the full-cohort camera FDR (the published-equivalent number)
#   - the LODO range: worst (max) and best (min) FDR over all refits
#   - how many refits still clear FDR < 0.05, and which donor is most
#     influential (the one whose removal weakens the result most)
#   - the targeted two-donor exclusion (1039 HIVE + PML, 1181 multi-pathology)
#   LODO is the stronger claim: no single donor drives a result, not merely
#   those two.
#   microglia_set_tests() falls back to "~ condition" when the full covariate
#   design is rank-deficient, and removing a donor can trigger that. Every fit's
#   formula is recorded and a mismatch with the full fit is flagged
#   (design_varied), so a LODO range that mixes designs is never reported silently.
#
# Inputs:  cached microglia pseudobulk per region (microglia_pseudobulk());
#          with --sets modules, OUT_ROOT/microglia_module_analysis/selected_gene_modules.csv
# Outputs: OUT_ROOT/03-DE_analysis/donor_sensitivity/
#          <REGION>_microglia_lodo_allfits[_modules].csv, microglia_lodo_summary[_modules].csv
#
# Usage: revision/docker/run.sh revision/tools/donor_lodo_microglia.R [--sets signatures|modules]

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

## --sets modules: the paper's own CS-CORE gene modules
##   (microglia_module_analysis/selected_gene_modules.csv; Fig 5d-g) instead of the
##   published signatures. The reviewer's "lipid/hypoxia module" is module 3.
##   Run as a SEPARATE camera call: camera's FDR is adjusted across the sets in one
##   call, so adding the modules to the signature list would change the existing
##   signature FDRs. Outputs get a _modules suffix.
##   Caveat for the text: the modules were built from PWH-vs-PWoH DEGs of these
##   same data, so the full-cohort test is not independent evidence; what the
##   LODO / exclusion refits show is whether the module still separates the groups
##   without a given donor.
.args <- commandArgs(trailingOnly = TRUE)
SETS  <- if ("--sets" %in% .args) .args[match("--sets", .args) + 1] else "signatures"
OUT_SFX <- if (SETS == "modules") "_modules" else ""
if (SETS == "modules") {
  .m <- utils::read.csv(file.path(OUT_ROOT, "microglia_module_analysis", "selected_gene_modules.csv"),
                        stringsAsFactors = FALSE)
  .lab <- c(hiv_down_module = "Module1_homeostatic_down",
            hiv_up_module1  = "Module2_inflammatory_up",
            hiv_up_module2  = "Module3_hypoxia_lipid_up")
  stopifnot(all(unique(.m$module) %in% names(.lab)))
  ## same shape as revision/config/microglia_signatures.R: list(genes =, direction =, source =)
  MICROGLIA_SIGNATURES <- lapply(split(.m$gene, .lab[.m$module]), function(g)
    list(genes = unique(g), direction = NA_character_, source = "paper CS-CORE module"))
  message("[lodo] gene sets: paper CS-CORE modules (",
          paste(sprintf("%s %d", names(MICROGLIA_SIGNATURES),
                        vapply(MICROGLIA_SIGNATURES, function(s) length(s$genes), integer(1))), collapse = ", "), ")")
}
source(file.path(SCORCH_REPO, "revision/tools/microglia_signature_score.R"))

suppressPackageStartupMessages({library(edgeR); library(limma)})

LODO_OUT <- file.path(OUT_ROOT, "03-DE_analysis", "donor_sensitivity")
dir.create(LODO_OUT, recursive = TRUE, showWarnings = FALSE)

## The two donors the reviewer names. 1039 = the only HIV encephalitis case
## (brain_pathology_v2 "HIV Encephalitis, Progressive multifocal
## leukoencephalopathy"); 1181 = "Other infections, Other non-infectious
## pathology". Both are PWH, so excluding them only shrinks the HIV arm.
PATHOLOGY_DONORS <- c("1039", "1181")

FDR_THRESHOLD <- 0.05

## ---- one fit on a donor subset ----------------------------------------------
#' @param pb    the cached pseudobulk list (counts, meta, region)
#' @param drop  donor ids to exclude
#' @return camera table with the design formula used, or NULL if the fit failed
.fit_subset <- function(pb, drop = character(0)) {
  keep <- !(pb$meta$donor %in% drop)
  if (sum(keep) < 6) return(NULL)          # too few donors to say anything
  sub <- list(counts = pb$counts[, keep, drop = FALSE],
              meta   = pb$meta[keep, , drop = FALSE],
              region = pb$region)
  out <- tryCatch(
    suppressMessages(microglia_set_tests(sub)),
    error = function(e) NULL)
  if (is.null(out)) return(NULL)
  data.frame(signature = out$camera$signature,
             direction = out$camera$Direction,
             pvalue    = out$camera$PValue,
             fdr       = out$camera$FDR,
             design    = out$design,
             n_donors  = out$n_donors,
             stringsAsFactors = FALSE)
}

## ---- LODO for one region -----------------------------------------------------
lodo_region <- function(region) {
  message("\n========== ", region, " ==========")
  pb <- microglia_pseudobulk(region)        # cached; no Seurat load
  donors <- pb$meta$donor
  message("[lodo] ", region, ": ", length(donors), " donors with >= 10 microglia")

  full <- .fit_subset(pb)
  if (is.null(full)) { message("[lodo] full fit failed for ", region); return(NULL) }
  full_design <- unique(full$design)

  ## every leave-one-out refit
  fits <- lapply(donors, function(d) {
    f <- .fit_subset(pb, drop = d)
    if (!is.null(f)) f$dropped <- d
    f
  })
  names(fits) <- donors
  ok <- !vapply(fits, is.null, logical(1))
  if (any(!ok))
    message("[lodo] ", sum(!ok), " refit(s) failed: ",
            paste(donors[!ok], collapse = ", "))
  lodo <- do.call(rbind, fits[ok])

  ## the targeted two-donor exclusion
  present <- intersect(PATHOLOGY_DONORS, donors)
  excl <- if (length(present)) .fit_subset(pb, drop = present) else NULL
  if (length(present) < length(PATHOLOGY_DONORS))
    message("[lodo] note: of ", paste(PATHOLOGY_DONORS, collapse = "/"),
            " only ", paste(present, collapse = ", "), " is in ", region,
            " microglia (others had < 10 cells or are absent)")

  ## collapse to one row per signature
  res <- do.call(rbind, lapply(unique(full$signature), function(s) {
    fu <- full[full$signature == s, ]
    lo <- lodo[lodo$signature == s, ]
    if (!nrow(lo)) return(NULL)
    worst_i <- which.max(lo$fdr)
    ex <- if (!is.null(excl)) excl[excl$signature == s, ] else NULL
    data.frame(
      region            = region,
      signature         = s,
      direction         = fu$direction,
      n_donors_full     = fu$n_donors,
      fdr_full          = fu$fdr,
      lodo_fdr_min      = min(lo$fdr),
      lodo_fdr_max      = max(lo$fdr),
      lodo_fdr_median   = stats::median(lo$fdr),
      lodo_n_fits       = nrow(lo),
      lodo_n_sig        = sum(lo$fdr < FDR_THRESHOLD),
      lodo_all_sig      = all(lo$fdr < FDR_THRESHOLD),
      most_influential  = lo$dropped[worst_i],
      direction_stable  = length(unique(lo$direction)) == 1 &&
                          identical(unique(lo$direction), fu$direction),
      fdr_excl_pathology = if (!is.null(ex) && nrow(ex)) ex$fdr else NA_real_,
      excluded_donors   = paste(present, collapse = ","),
      design_full       = full_design,
      design_varied     = length(unique(c(full_design, lo$design))) > 1,
      stringsAsFactors  = FALSE)
  }))
  res <- res[order(res$fdr_full), ]

  ## per-fit detail, for anyone who wants to see the whole distribution
  utils::write.csv(lodo, file.path(LODO_OUT, paste0(region, "_microglia_lodo_allfits", OUT_SFX, ".csv")),
                   row.names = FALSE)
  res
}

## ---- run ---------------------------------------------------------------------
all_res <- do.call(rbind, lapply(c("INS", "PFC", "VST"), function(r)
  tryCatch(lodo_region(r), error = function(e) {
    message("[lodo] ", r, " failed: ", conditionMessage(e)); NULL
  })))

utils::write.csv(all_res, file.path(LODO_OUT, paste0("microglia_lodo_summary", OUT_SFX, ".csv")),
                 row.names = FALSE)
message("\n[lodo] wrote ", file.path(LODO_OUT, paste0("microglia_lodo_summary", OUT_SFX, ".csv")))

## ---- console summary ---------------------------------------------------------
message("\n########## LODO: signatures significant in the full cohort ##########")
sig <- all_res[all_res$fdr_full < FDR_THRESHOLD, ]
if (!nrow(sig)) {
  message("nothing significant at FDR < ", FDR_THRESHOLD, " in the full cohort")
} else {
  pr <- data.frame(
    region      = sig$region,
    signature   = substr(sig$signature, 1, 26),
    dir         = sig$direction,
    fdr_full    = signif(sig$fdr_full, 3),
    lodo_worst  = signif(sig$lodo_fdr_max, 3),
    n_sig       = paste0(sig$lodo_n_sig, "/", sig$lodo_n_fits),
    robust      = ifelse(sig$lodo_all_sig, "YES", "no"),
    excl_2donor = signif(sig$fdr_excl_pathology, 3),
    stringsAsFactors = FALSE)
  print(pr, row.names = FALSE)

  nr <- sig[!sig$lodo_all_sig, ]
  if (nrow(nr)) {
    message("\nNOT robust to every single-donor removal:")
    for (i in seq_len(nrow(nr)))
      message("  ", nr$region[i], " ", nr$signature[i],
              ": worst FDR ", signif(nr$lodo_fdr_max[i], 3),
              " when dropping donor ", nr$most_influential[i])
  }
  if (any(!sig$direction_stable))
    message("\nWARNING: direction changed in some refits - inspect the allfits csv")
  if (any(sig$design_varied))
    message("\nWARNING: the covariate design differed between some fits ",
            "(rank-deficiency fallback). LODO range is not strictly comparable ",
            "to the full fit for those rows - see design_varied.")
}
message("\n[lodo] done")
