# ==============================================================================
# Astrocyte cluster and marker supplemental figures (dot plots, labelled tSNEs,
# cluster sizes and marker tables).
#
# Method
#   Driver for astro_save_plots() and astro_features() in
#   revision/tools/astro_cluster_features.R, grouped by named astrocyte subtype
#   (subcluster_labels) rather than cluster number. Uses the astrocyte object
#   without donor 1040's VST astrocytes (33 of 74,894 cells; see
#   revision/tools/astro_drop_1040_vst.R) and warns if only the unfiltered
#   object is available. The Reactive-PWH marker panel is
#   GFAP/CD44/OSMR/SERPINA3/CHI3L1/TNC, the definition used in the paper,
#   rather than interferon-response genes.
#
# Inputs
#   OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled_no1040vst.RDS
#   (falls back to astro_all_clusters_labeled.RDS; --in overrides)
# Outputs  (OUT_ROOT/06-astrocyte_analysis/figures_cluster<suffix>/; suffix
#          defaults to _no1040vst for the filtered object)
#   astro_save_plots():
#     astro_marker_dotplot_annotated.png   markers x subtype, grouped by panel
#     astro_tsne_labeled.png               tSNE labelled by subtype
#     astro_tsne_labeled_by_condition.png  the same, split by condition
#     astro_cluster_sizes.csv
#   astro_features():
#     astro_marker_dotplot.png             FindAllMarkers-derived markers
#     astro_top_marker_dotplot.png         top N per group
#     plus the marker tables behind them
#
# Usage
#   revision/docker/run.sh revision/tools/run_astro_cluster_figures.R
#       [--in <file>] [--suffix _x] [--group <column>] [--skip-findmarkers]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_cluster_features.R"))

suppressPackageStartupMessages({library(Seurat)})

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
SKIP_FM <- "--skip-findmarkers" %in% args

rerun_setup("06-astrocyte_analysis", subdirs = c("data", "figures"))

## Prefer the 1040-filtered object; announce a fallback rather than sliding into
## it, since the two objects are indistinguishable from their outputs.
.cands <- c(
  if (nzchar(getarg("--in", ""))) file.path("./data", getarg("--in", "")) else NA_character_,
  "./data/astro_all_clusters_labeled_no1040vst.RDS",
  "./data/astro_all_clusters_labeled.RDS"
)
IN <- .cands[!is.na(.cands) & file.exists(.cands)][1]
if (is.na(IN)) stop("no astrocyte object found under ./data")
if (!grepl("no1040", IN))
  warning("[astro-fig] ", basename(IN), " still contains donor 1040's VST cells",
          immediate. = TRUE)

SFX <- getarg("--suffix", if (grepl("no1040", IN)) "_no1040vst" else "")

message("[astro-fig] input : ", IN)
message("[astro-fig] suffix: ", if (nzchar(SFX)) SFX else "(none)")

obj <- readRDS(IN)
md  <- obj@meta.data
message("[astro-fig] ", ncol(obj), " cells")

## Group by subtype, not by cluster number: the 21-level harmony_astro_clusters
## gives axes running 0-20, which do not tell the reader which row is
## Reactive-PWH. subcluster_labels carries the named subtypes. Pass --group to
## override.
.grp_default <- "subcluster_labels"
grp <- getarg("--group", .grp_default)
if (!grp %in% colnames(md)) {
  alt <- intersect(c("subcluster_labels", "harmony_astro_clusters", "seurat_clusters"),
                   colnames(md))
  if (!length(alt)) stop("no grouping column found on the object")
  message("[astro-fig] '", grp, "' not present; falling back to '", alt[1], "'")
  grp <- alt[1]
}
message("[astro-fig] grouping column: ", grp, " (",
        length(unique(as.character(md[[grp]]))), " levels: ",
        paste(sort(unique(as.character(md[[grp]]))), collapse = ", "), ")")

## ---- outputs go to a suffixed directory --------------------------------------
## astro_save_plots()/astro_features() hardcode their filenames, so the suffix
## is carried by the directory instead of rewriting both functions.
OUTDIR <- if (nzchar(SFX))
  file.path(OUT_ROOT, "06-astrocyte_analysis", paste0("figures_cluster", SFX)) else
  file.path(OUT_ROOT, "06-astrocyte_analysis")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
message("[astro-fig] writing to ", OUTDIR)

## ---- 1. dot plot + the two tSNEs ---------------------------------------------
message("\n=================== astro_save_plots ===================")
r1 <- try(astro_save_plots(obj, group = grp, outdir = OUTDIR), silent = FALSE)
if (inherits(r1, "try-error"))
  message("[astro-fig] FAILED: ", conditionMessage(attr(r1, "condition")))

## ---- 2. FindAllMarkers-derived dot plots -------------------------------------
## This is the slow half - FindAllMarkers over ~21 clusters on ~73k cells.
if (SKIP_FM) {
  message("\n[astro-fig] --skip-findmarkers given; not running astro_features()")
} else {
  message("\n=================== astro_features ===================")
  r2 <- try(astro_features(obj, group = grp, outdir = OUTDIR), silent = FALSE)
  if (inherits(r2, "try-error"))
    message("[astro-fig] FAILED: ", conditionMessage(attr(r2, "condition")))
}

message("\n[astro-fig] files now in ", OUTDIR, ":")
print(sort(basename(list.files(OUTDIR, pattern = "\\.(png|csv)$"))))
message("\n[astro-fig] done")
