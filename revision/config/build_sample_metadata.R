#!/usr/bin/env Rscript
# ==============================================================================
# build_sample_metadata.R: builds the sample sheet (sample_metadata_v3.csv)
# that the published scripts read.
#
# Method: concatenates the CellBender pipeline's per-region manifests, which
# carry the same 76 columns as the published sample sheet, including every
# field the scripts use:
#     sample_name, cellranger_path, region, oud_jay_v2, hiv_jay_v2,
#     sample_ids, age_value_v2, sex_v2, race_v2, postmortem_interval_v2,
#     subject_source_v2
# Boolean columns are normalised to logicals. The result is validated: the
# scripts' own filter
#     !is.na(cellranger_path) & oud_jay_v2 == FALSE & region == <REGION>
# must reproduce the published sample counts and match the sample folders on
# disk exactly, otherwise the script stops:
#     INS 53 -> 31 (19 PWH / 12 PWoH)
#     PFC 53 -> 31 (19 PWH / 12 PWoH)
#     VST 56 -> 33 (20 PWH / 13 PWoH)
#
# Inputs:  PIPELINE_DIR/{ins,pfc,vst}/*_ambientRNA*.csv manifests.
# Outputs: SAMPLE_META (OUT_ROOT/meta/sample_metadata_v3.csv). Not distributed:
#          contains subject-level HIV status, toxicology and cause of death.
#
# Usage: revision/docker/run.sh revision/config/build_sample_metadata.R
# Run once before any stage-01 script.
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

MANIFESTS <- c(
  INS = file.path(PIPELINE_DIR, "ins", "ins_ambientRNA.csv"),
  PFC = file.path(PIPELINE_DIR, "pfc", "pfc_ambientRNA_modified.csv"),
  VST = file.path(PIPELINE_DIR, "vst", "vst_ambientRNA.csv")
)

# Manifests disagree on boolean capitalisation: INS/VST write True/False,
# PFC writes TRUE/FALSE. R's as.logical() accepts both, but normalise
# explicitly so a downstream `== FALSE` can never silently match nothing.
BOOL_COLS <- c("oud_jay_v2", "hiv_jay_v2", "bad_quality",
               "overdose_death_v2", "transgender_v2")

norm_logical <- function(x) {
  if (is.logical(x)) return(x)
  v <- trimws(as.character(x))
  out <- rep(NA, length(v))
  out[tolower(v) %in% c("true",  "t", "yes")] <- TRUE
  out[tolower(v) %in% c("false", "f", "no")]  <- FALSE
  out
}

parts <- lapply(names(MANIFESTS), function(reg) {
  f <- MANIFESTS[[reg]]
  if (!file.exists(f)) stop("missing manifest: ", f)
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  names(d) <- trimws(names(d))          # PFC's last header cell carries a CR
  cat(sprintf("  %s: %d rows, %d cols  (%s)\n", reg, nrow(d), ncol(d), basename(f)))
  d
})
names(parts) <- names(MANIFESTS)

# Columns must line up before rbind.
ref_cols <- names(parts[[1]])
for (reg in names(parts)) {
  extra   <- setdiff(names(parts[[reg]]), ref_cols)
  missing <- setdiff(ref_cols, names(parts[[reg]]))
  if (length(extra) || length(missing)) {
    stop(reg, " column mismatch.\n  extra: ", paste(extra, collapse = ", "),
         "\n  missing: ", paste(missing, collapse = ", "))
  }
  parts[[reg]] <- parts[[reg]][, ref_cols, drop = FALSE]
}

meta <- do.call(rbind, parts)
rownames(meta) <- NULL

for (b in intersect(BOOL_COLS, names(meta))) meta[[b]] <- norm_logical(meta[[b]])
meta$region <- toupper(trimws(meta$region))

# ---- Validate against the counts the original scripts assert -----------------
EXPECTED <- list(INS = c(total = 53, keep = 31, hiv = 19, ctr = 12),
                 PFC = c(total = 53, keep = 31, hiv = 19, ctr = 12),
                 VST = c(total = 56, keep = 33, hiv = 20, ctr = 13))

cat("\nFilter check (!is.na(cellranger_path) & !oud_jay_v2 & region == R):\n")
ok <- TRUE
for (reg in names(EXPECTED)) {
  sub  <- meta[meta$region == reg, , drop = FALSE]
  # `== FALSE` yields NA for an unparsed boolean, and NA subscripts inject
  # phantom rows. `%in% FALSE` is NA-safe.
  keep <- sub[!is.na(sub$cellranger_path) & sub$cellranger_path != "" &
              sub$oud_jay_v2 %in% FALSE, , drop = FALSE]
  got <- c(total = nrow(sub), keep = nrow(keep),
           hiv = sum(keep$hiv_jay_v2 == TRUE), ctr = sum(keep$hiv_jay_v2 == FALSE))
  exp <- EXPECTED[[reg]]
  match <- identical(as.integer(got), as.integer(exp))
  ok <- ok && match
  cat(sprintf("  %s  total=%d keep=%d HIV=%d CTR=%d   %s\n",
              reg, got["total"], got["keep"], got["hiv"], got["ctr"],
              if (match) "OK" else paste0("MISMATCH, expected ",
                                          paste(exp, collapse = "/"))))

  # The kept samples must equal the sample folders present on disk.
  on_disk <- basename(list.dirs(region_dir(reg), recursive = FALSE))
  if (!setequal(keep$sample_name, on_disk)) {
    ok <- FALSE
    cat("     sample_name vs folders differ — only in manifest: ",
        paste(setdiff(keep$sample_name, on_disk), collapse = ", "),
        " | only on disk: ",
        paste(setdiff(on_disk, keep$sample_name), collapse = ", "), "\n")
  }
}
if (!ok) stop("sample metadata validation FAILED — do not proceed.")

dir.create(dirname(SAMPLE_META), recursive = TRUE, showWarnings = FALSE)
write.csv(meta, SAMPLE_META, row.names = FALSE)
cat("\nWrote ", SAMPLE_META, " (", nrow(meta), " rows x ", ncol(meta), " cols)\n", sep = "")
cat("NOTE: contains subject-level HIV status, toxicology and cause of death.\n")
cat("      It is gitignored by design — do not force-add it.\n")
