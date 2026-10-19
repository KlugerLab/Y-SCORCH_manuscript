#!/usr/bin/env Rscript
# ==============================================================================
# Astrocyte figures and tables for the reviewer's questions on the Reactive-PWH
# subtype: per-donor contribution, donor-level abundance, QC/doublet metrics,
# donor-aware marker DE, and stability across clustering resolutions.
#
# Method
#   Batch driver for the functions in astro_composition_donor.R,
#   astro_subtype_de_donor.R and astro_reviewer_figures.R (revision/tools/).
#   The "clean" subset drops the labels given by --drop (default Oligo and
#   Endo-Astro: contamination clusters rather than astrocyte states, as in the
#   astrocyte Rmd) and is rebuilt here from subcluster_labels. The
#   resolution-stability analysis re-clusters at seven resolutions and is slow;
#   --only runs a subset of the analyses (ask1, ask2, ask4, subtype_de, ask3,
#   ask3plot), and ask3plot redraws that figure from its saved table.
#
# Inputs
#   OUT_ROOT/06-astrocyte_analysis/data/<--in>  (default
#   astro_all_clusters_labeled.RDS; use astro_all_clusters_labeled_no1040vst.RDS,
#   which excludes donor 1040's VST astrocytes, for manuscript figures)
# Outputs
#   OUT_ROOT/06-astrocyte_analysis/figures/  figures and CSV tables, file names
#   carrying --suffix
#
# Usage
#   revision/docker/run.sh revision/tools/run_astro_reviewer_figures.R
#       [--in <file>] [--suffix _no1040vst] [--subtype "Reactive-PWH"]
#       [--drop Oligo,Endo-Astro] [--only ask1,ask2,...]
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(ggplot2); library(dplyr); library(tidyr)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
TOOLS <- file.path(SCORCH_REPO, "revision/tools")
for (f in c("astro_pub_theme", "astro_cluster_features", "astro_composition_donor",
            "astro_subtype_de_donor", "astro_reviewer_figures"))
  source(file.path(TOOLS, paste0(f, ".R")))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}
SUBTYPE <- getarg("--subtype", "Reactive-PWH")
## clusters that are contamination rather than astrocyte states; dropped for the
## "clean" analyses exactly as in the Rmd
DROP <- strsplit(getarg("--drop", "Oligo,Endo-Astro"), ",")[[1]]

rerun_setup("06-astrocyte_analysis", subdirs = c("figures", "data"))
## --in selects the input object. The manuscript analyses use
##   astro_all_clusters_labeled_no1040vst.RDS (donor 1040's VST astrocytes
##   excluded; see revision/tools/astro_drop_1040_vst.R); pass it explicitly,
##   since the default is the unfiltered object.
IN <- file.path("./data", getarg("--in", "astro_all_clusters_labeled.RDS"))
if (!file.exists(IN)) stop("missing ", IN, " - save the labelled object from the Rmd first")
message("[run] input object: ", IN)

## --suffix is appended to every output filename, so each figure records which
##   sample set produced it.
##   Use "_no1040vst" for runs on astro_all_clusters_labeled_no1040vst.RDS.
##   ASTRO_OUT_SUFFIX is read by .sfxname() in the three sourced tool files.
ASTRO_OUT_SUFFIX <- getarg("--suffix", "")
if (nzchar(ASTRO_OUT_SUFFIX))
  message("[run] output suffix: ", ASTRO_OUT_SUFFIX)

message("[run] loading ", IN, " ...")
obj <- readRDS(IN)
LAB <- if ("subcluster_labels" %in% colnames(obj@meta.data)) "subcluster_labels" else
  stop("no subcluster_labels column on the object")
message("[run] ", ncol(obj), " cells; labels: ",
        paste(sprintf("%s=%d", names(table(obj@meta.data[[LAB]])),
                      table(obj@meta.data[[LAB]])), collapse = "  "))

keep_lab <- setdiff(unique(as.character(obj@meta.data[[LAB]])), DROP)
message("[run] dropping: ", paste(intersect(DROP, unique(as.character(obj@meta.data[[LAB]]))),
                                  collapse = ", "))
Idents(obj) <- LAB
clean <- subset(obj, idents = keep_lab)
message("[run] clean subset: ", ncol(clean), " cells (", ncol(obj) - ncol(clean), " removed)")
if (!SUBTYPE %in% as.character(clean@meta.data[[LAB]]))
  stop("subtype '", SUBTYPE, "' not present after the drop")

ok <- function(label, expr) {
  message("\n=================== ", label, " ===================")
  r <- try(expr, silent = FALSE)
  if (inherits(r, "try-error"))
    message("[run] FAILED: ", label, " -- ", conditionMessage(attr(r, "condition")))
  invisible(r)
}

## --only ask1,ask2,ask4,subtype_de,ask3 runs a subset, so redrawing
##   one figure does not repeat the resolution re-clustering (ask 3).
ONLY <- strsplit(getarg("--only", ""), ",")[[1]]
want <- function(k) !length(ONLY) || k %in% ONLY
if (length(ONLY)) message("[run] only: ", paste(ONLY, collapse = ", "))

## ask 1 - per-donor contribution, and donor origin of each subtype
if (want("ask1"))
  ok("ask 1: per-donor contribution",
     astro_fig_donor_contribution(clean, label_col = LAB))

## ask 2 - donor-level abundance statistics
if (want("ask2"))
  ok("ask 2: donor-level abundance",
     astro_composition(clean, label_col = LAB))

## ask 4 - QC / doublet metrics (joins scDblFinder back from stage-01)
if (want("ask4"))
  qc <- ok("ask 4: QC and doublet metrics", astro_fig_qc_metrics(clean, label_col = LAB))

## donor-aware marker DE for the target subtype, plus the per-donor heatmap
if (want("subtype_de"))
  ok(paste("donor-aware marker DE:", SUBTYPE),
     astro_subtype_de(clean, subtype = SUBTYPE, label_col = LAB))

## ask 3 redrawn from its saved table - no re-clustering
if (want("ask3plot"))
  ok(paste("ask 3 (replot from saved table):", SUBTYPE), {
    f3 <- file.path(OUT_ROOT, "06-astrocyte_analysis", "figures",
                    .sfxname(paste0("ask3_resolution_stability_",
                                    gsub("[^A-Za-z0-9]+", "_", SUBTYPE), ".csv")))
    if (!file.exists(f3)) stop("no saved ask 3 table: ", f3, " - run --only ask3")
    astro_plot_resolution_stability(utils::read.csv(f3), SUBTYPE)
  })

## ask 3 - the slow one: re-clusters at several resolutions
if (want("ask3"))
  ok(paste("ask 3: resolution stability:", SUBTYPE),
     astro_fig_resolution_stability(clean, subtype = SUBTYPE, label_col = LAB))

message("\n[run] DONE - figures and tables in ",
        file.path(OUT_ROOT, "06-astrocyte_analysis"))
