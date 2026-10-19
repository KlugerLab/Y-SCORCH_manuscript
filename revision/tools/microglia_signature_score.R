# ==============================================================================
# Donor-level scoring and testing of microglial state signatures (functions).
#
# Every statistic here has n = donors: scoring per nucleus and testing across
# nuclei would be pseudoreplication (40,000 nuclei are not 31 people).
#
# Method: a donor's microglia are summed into one pseudobulk profile before
# scoring. This recovers genes that are sparse in any single nucleus (mitigating
# the Thrupp 2020 depletion of activation genes such as APOE, CST3, SPP1, CD74),
# makes the 4-gene sets usable, and puts the scores on the same footing as the
# pseudobulk DE. Counts are uncorrected counts within CellBender-defined
# microglia, as in the DE: CellBender defines the cell populations, not the
# expression magnitudes. Donors need >= 10 microglia.
#   score   per-donor mean z (each gene z-scored across donors on logCPM,
#           averaged over the set); one value per person.
#   camera  competitive - is the set changed relative to all other genes?
#           Corrects for inter-gene correlation, which matters because
#           DAM_stage2 / DAM_Haney / MGnD share most of their genes.
#   fry     self-contained - is the set changed at all? A set significant by fry
#           but not camera moved with the whole transcriptome, not specifically.
#   Design per region: ~ condition + Age + Sex + Race + PMI (n = 31/31/33 donors).
#   Cross-region: ~ condition + region + covariates + (1 | donor) on the scores;
#   a donor contributes up to 3 regions, so donor is a random effect.
#   Sensitivity: every set is also scored with AMBIENT_SUSPECT genes removed
#   (e.g. IRM_core loses IFITM3 and ISG15); an interferon result that does not
#   survive this is not reported.
#
# Inputs: regional objects (via .ambient_region_object() in revision/tools/ambient_fraction.R),
#         SAMPLE_META, gene sets in revision/config/microglia_signatures.R
# Output: OUT_ROOT/03-DE_analysis/microglia_signatures/ - <R>_microglia_pseudobulk.rds,
#         <R>_{detection,scores,camera,fry}.csv, *_dropAmbient.csv, cross_region_mixedmodel.csv
#
# Usage (sourced; driven by revision/tools/run_microglia_signatures.R):
#   source("revision/config/paths.R")
#   source("revision/tools/microglia_signature_score.R")
#   pb  <- microglia_pseudobulk("INS")
#   res <- microglia_score_region(pb)
#   all <- microglia_score_all()        # 3 regions + cross-region pooled
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(limma); library(edgeR)
  library(dplyr); library(ggplot2)
})

## Sibling files this depends on. Sourced here rather than left to the caller:
## .ambient_region_object() encodes the per-region object paths (VST reads from
## a different stage directory than INS/PFC), and getting that wrong would score
## a cell set that does not correspond to the DE - silently.
.mgsig_dir <- file.path(SCORCH_REPO, "revision")
if (!exists("MICROGLIA_SIGNATURES"))
  source(file.path(.mgsig_dir, "config", "microglia_signatures.R"))
if (!exists(".ambient_region_object"))
  source(file.path(.mgsig_dir, "tools", "ambient_fraction.R"))

MICRO_LABEL   <- "Micro"                  # value in consensus_annotation_v3
MICRO_CT_COL  <- "consensus_annotation_v3"
MICRO_SID_COL <- "sample_ids"
MICRO_OUT     <- file.path(OUT_ROOT, "03-DE_analysis", "microglia_signatures")

# ------------------------------------------------------------------------------
#' Microglia pseudobulk for one region, from UNCORRECTED counts.
#'
#' @return list(counts = genes x donors, meta = donor table, region = region)
#' Cached pseudobulk path. Building it costs a 2.4-3 GB object read, and every
#' gene-level plot needs the same matrix, so it is written once and reused.
microglia_pb_path <- function(region)
  file.path(MICRO_OUT, paste0(toupper(region), "_microglia_pseudobulk.rds"))

microglia_pseudobulk <- function(region, obj = NULL, min_cells = 10,
                                 cache = TRUE) {
  region <- toupper(region)
  cp <- microglia_pb_path(region)
  if (cache && is.null(obj) && file.exists(cp)) {
    message("[mgsig] cached pseudobulk: ", cp)
    return(readRDS(cp))
  }
  if (is.null(obj)) {
    p <- .ambient_region_object(region)     # revision/tools/ambient_fraction.R
    if (!file.exists(p)) stop("missing regional object: ", p)
    message("[mgsig] loading ", p)
    obj <- readRDS(p)
  }
  stopifnot(MICRO_CT_COL %in% colnames(obj@meta.data),
            MICRO_SID_COL %in% colnames(obj@meta.data))

  ct <- as.character(obj@meta.data[[MICRO_CT_COL]])
  if (!MICRO_LABEL %in% ct)
    stop("no cells labelled '", MICRO_LABEL, "' in ", MICRO_CT_COL,
         ". Present: ", paste(sort(unique(ct)), collapse = ", "))
  idx <- which(ct == MICRO_LABEL)
  message("[mgsig] ", region, ": ", length(idx), " microglia")

  ## uncorrected counts, as in the DE
  raw <- uncorrected_counts(obj, region)[, idx, drop = FALSE]
  sid <- as.character(obj@meta.data[[MICRO_SID_COL]])[idx]

  n_per <- table(sid)
  keep_d <- names(n_per)[n_per >= min_cells]
  dropped <- setdiff(names(n_per), keep_d)
  if (length(dropped))
    message("[mgsig] dropping ", length(dropped), " donor(s) with < ", min_cells,
            " microglia: ", paste(dropped, collapse = ", "))
  sel <- sid %in% keep_d
  raw <- raw[, sel, drop = FALSE]; sid <- sid[sel]

  f  <- factor(sid, levels = keep_d)
  pb <- as.matrix(raw %*% Matrix::sparse.model.matrix(~ 0 + f))
  colnames(pb) <- levels(f)
  message("[mgsig] pseudobulk: ", nrow(pb), " genes x ", ncol(pb), " donors")

  meta <- data.frame(donor = colnames(pb),
                     n_cells = as.integer(n_per[colnames(pb)]),
                     region = region, stringsAsFactors = FALSE)
  meta <- .microglia_attach_meta(meta)
  out <- list(counts = pb, meta = meta, region = region)
  if (cache) {
    dir.create(dirname(cp), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cp)
    message("[mgsig] cached pseudobulk -> ", cp)
  }
  out
}

#' Donor covariates from the sample sheet.
.microglia_attach_meta <- function(meta) {
  if (!file.exists(SAMPLE_META)) {
    warning("no sample metadata at ", SAMPLE_META, " - condition unavailable")
    return(meta)
  }
  sm <- utils::read.csv(SAMPLE_META, stringsAsFactors = FALSE)
  sm <- sm[!duplicated(sm$sample_ids), ]
  rownames(sm) <- as.character(sm$sample_ids)
  miss <- setdiff(meta$donor, rownames(sm))
  if (length(miss)) stop("donors absent from the sample sheet: ",
                         paste(miss, collapse = ", "))
  meta$condition <- factor(as.character(pwh_relabel(sm[meta$donor, "hiv_jay_v2"])),
                           levels = c("PWoH", "PWH"))
  meta$Age <- as.numeric(sm[meta$donor, "age_value_v2"])
  meta$Sex <- factor(make.names(sm[meta$donor, "sex_v2"]))
  meta$PMI <- as.numeric(sm[meta$donor, "postmortem_interval_v2"])
  r <- make.names(sm[meta$donor, "race_v2"]); r[grepl("Native", r)] <- "Other"
  meta$Race <- factor(r)
  meta
}

# ------------------------------------------------------------------------------
#' Which signature genes are actually measurable here, and how compromised.
#'
#' Printed with every result. A set whose genes are largely undetected has not
#' been "tested and found negative" - it has not been tested.
microglia_detection_report <- function(pb, sigs = MICROGLIA_SIGNATURES,
                                       ambient_tbl = NULL) {
  present <- rownames(pb$counts)
  det <- present[rowSums(pb$counts > 0) >= max(2, 0.2 * ncol(pb$counts))]
  do.call(rbind, lapply(names(sigs), function(n) {
    g <- sigs[[n]]$genes
    amb <- NA_real_
    if (!is.null(ambient_tbl)) {
      a <- ambient_tbl$ambient_frac[ambient_tbl$gene %in% g &
                                    ambient_tbl$cluster_id == MICRO_LABEL]
      amb <- if (length(a)) max(a, na.rm = TRUE) else NA_real_
    }
    data.frame(signature = n, n_genes = length(g),
               n_in_matrix = sum(g %in% present),
               n_detected  = sum(g %in% det),
               detect_rate = round(sum(g %in% det) / length(g), 3),
               missing = paste(setdiff(g, det), collapse = "/"),
               thrupp = paste(intersect(g, THRUPP_DEPLETED), collapse = "/"),
               max_ambient = round(amb, 3),
               stringsAsFactors = FALSE)
  }))
}

# ------------------------------------------------------------------------------
#' Per-donor score: genes z-scored across donors, averaged over the set.
#'
#' Uses logCPM so donors with more microglia are not systematically higher.
#' Genes absent from the matrix are dropped and counted, never zero-filled - a
#' zero-filled gene would drag the mean toward the middle and look like a
#' moderate score rather than a missing measurement.
microglia_scores <- function(pb, sigs = MICROGLIA_SIGNATURES, drop_ambient = FALSE) {
  y <- DGEList(counts = pb$counts)
  y <- calcNormFactors(y)
  lc <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  z  <- t(scale(t(lc)))                     # per gene, across donors
  z  <- z[is.finite(rowSums(z)), , drop = FALSE]

  sc <- sapply(names(sigs), function(n) {
    g <- sigs[[n]]$genes
    if (drop_ambient) g <- setdiff(g, AMBIENT_SUSPECT)
    g <- intersect(g, rownames(z))
    if (length(g) < 2) return(rep(NA_real_, ncol(z)))
    colMeans(z[g, , drop = FALSE])
  })
  sc <- as.data.frame(sc)
  sc$donor <- colnames(z)
  merge(pb$meta, sc, by = "donor", sort = FALSE)
}

# ------------------------------------------------------------------------------
#' camera + fry for one region, on the same design as the score.
microglia_set_tests <- function(pb, sigs = MICROGLIA_SIGNATURES,
                                drop_ambient = FALSE) {
  meta <- pb$meta
  if (nlevels(droplevels(meta$condition)) < 2)
    stop("only one condition present in ", pb$region)

  form <- "~ condition"
  if (all(!is.na(meta$Age)) && all(!is.na(meta$PMI)) &&
      nlevels(droplevels(meta$Sex)) > 1 && nlevels(droplevels(meta$Race)) > 1 &&
      nrow(meta) > 10)
    form <- "~ condition + Age + Sex + Race + PMI"
  design <- model.matrix(stats::as.formula(form), data = meta)
  if (qr(design)$rank < ncol(design)) {
    message("[mgsig] design rank-deficient; falling back to ~ condition")
    design <- model.matrix(~ condition, data = meta)
    form <- "~ condition"
  }
  message("[mgsig] ", pb$region, " design: ", form, "  (n = ", nrow(meta), ")")

  y <- DGEList(counts = pb$counts)
  k <- filterByExpr(y, design)
  y <- y[k, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  v <- voom(y, design)

  idx <- lapply(sigs, function(s) {
    g <- s$genes
    if (drop_ambient) g <- setdiff(g, AMBIENT_SUSPECT)
    which(rownames(v) %in% g)
  })
  usable <- vapply(idx, length, integer(1)) >= 2
  if (!any(usable)) stop("no signature has >= 2 detected genes")
  if (any(!usable))
    message("[mgsig] not testable (< 2 genes detected): ",
            paste(names(idx)[!usable], collapse = ", "))
  idx <- idx[usable]

  cf <- grep("^conditionPWH$", colnames(design), value = TRUE)
  stopifnot(length(cf) == 1)

  cam <- camera(v, idx, design, contrast = cf)
  cam$signature <- rownames(cam); cam$test <- "camera"
  fr  <- fry(v, idx, design, contrast = cf)
  fr$signature  <- rownames(fr);  fr$test  <- "fry"

  list(camera = cam, fry = fr, design = form, n_donors = nrow(meta),
       n_genes_tested = vapply(idx, length, integer(1)))
}

# ------------------------------------------------------------------------------
#' Everything for one region: scores, tests, sensitivity, detection report.
microglia_score_region <- function(pb, sigs = MICROGLIA_SIGNATURES,
                                   ambient_tbl = NULL, outdir = MICRO_OUT) {
  det <- microglia_detection_report(pb, sigs, ambient_tbl)
  message("[mgsig] detection report:"); print(det[, c(1,2,4,5,7,8)], row.names = FALSE)

  sc_full <- microglia_scores(pb, sigs, drop_ambient = FALSE)
  sc_drop <- microglia_scores(pb, sigs, drop_ambient = TRUE)
  t_full  <- microglia_set_tests(pb, sigs, drop_ambient = FALSE)
  t_drop  <- microglia_set_tests(pb, sigs, drop_ambient = TRUE)

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    r <- pb$region
    utils::write.csv(det, file.path(outdir, paste0(r, "_detection.csv")), row.names = FALSE)
    utils::write.csv(sc_full, file.path(outdir, paste0(r, "_scores.csv")), row.names = FALSE)
    utils::write.csv(sc_drop, file.path(outdir, paste0(r, "_scores_dropAmbient.csv")), row.names = FALSE)
    utils::write.csv(t_full$camera, file.path(outdir, paste0(r, "_camera.csv")), row.names = FALSE)
    utils::write.csv(t_full$fry,    file.path(outdir, paste0(r, "_fry.csv")), row.names = FALSE)
    utils::write.csv(t_drop$camera, file.path(outdir, paste0(r, "_camera_dropAmbient.csv")), row.names = FALSE)
  }
  list(region = pb$region, detection = det, scores = sc_full,
       scores_drop = sc_drop, tests = t_full, tests_drop = t_drop)
}

# ------------------------------------------------------------------------------
#' Cross-region pooled test on the per-donor scores.
#'
#' A donor appears once per region, so donor is a RANDOM EFFECT. Pooling the
#' donor-regions as independent would treat 33 people as ~95 and re-introduce
#' the very inflation this whole approach exists to avoid.
microglia_cross_region <- function(score_list, sigs = MICROGLIA_SIGNATURES,
                                   outdir = MICRO_OUT) {
  if (!requireNamespace("lmerTest", quietly = TRUE))
    stop("lmerTest needed for the cross-region mixed model")
  all_sc <- dplyr::bind_rows(lapply(score_list, function(x) x$scores))
  message("[mgsig] cross-region: ", nrow(all_sc), " donor-region rows, ",
          length(unique(all_sc$donor)), " unique donors, ",
          length(unique(all_sc$region)), " regions")

  res <- do.call(rbind, lapply(names(sigs), function(n) {
    if (!n %in% colnames(all_sc) || all(is.na(all_sc[[n]])))
      return(NULL)
    d <- all_sc[!is.na(all_sc[[n]]), ]
    d$y <- d[[n]]
    fit <- try(lmerTest::lmer(y ~ condition + region + Age + Sex + Race + PMI +
                              (1 | donor), data = d), silent = TRUE)
    if (inherits(fit, "try-error"))
      fit <- try(lmerTest::lmer(y ~ condition + region + (1 | donor), data = d),
                 silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    s <- summary(fit)$coefficients
    row <- grep("^conditionPWH", rownames(s), value = TRUE)
    if (!length(row)) return(NULL)
    data.frame(signature = n, estimate = s[row, "Estimate"],
               se = s[row, "Std. Error"], t = s[row, "t value"],
               p = s[row, "Pr(>|t|)"],
               n_obs = nrow(d), n_donors = length(unique(d$donor)),
               stringsAsFactors = FALSE)
  }))
  if (is.null(res)) { message("[mgsig] no signature estimable cross-region"); return(NULL) }
  res$FDR <- stats::p.adjust(res$p, "BH")
  res <- res[order(res$p), ]
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, file.path(outdir, "cross_region_mixedmodel.csv"), row.names = FALSE)
  }
  res
}

# ------------------------------------------------------------------------------
#' Per-region primary + cross-region pooled, end to end.
microglia_score_all <- function(regions = c("INS", "PFC", "VST"),
                                sigs = MICROGLIA_SIGNATURES) {
  message("[mgsig] signature summary:")
  print(microglia_signature_summary(sigs), row.names = FALSE)
  message("\n[mgsig] Jaccard overlap - concordant scores among overlapping sets ",
          "are ONE result, not several:")
  print(microglia_signature_overlap(sigs))

  out <- list()
  for (r in regions) {
    message("\n########## ", r, " ##########")
    amb <- try(utils::read.csv(ambient_table_path(r), stringsAsFactors = FALSE),
               silent = TRUE)
    if (inherits(amb, "try-error")) amb <- NULL
    pb <- microglia_pseudobulk(r)
    out[[r]] <- microglia_score_region(pb, sigs, ambient_tbl = amb)
  }
  message("\n########## CROSS-REGION ##########")
  cr <- microglia_cross_region(out, sigs)
  if (!is.null(cr)) print(cr, row.names = FALSE)
  invisible(list(per_region = out, cross_region = cr))
}
