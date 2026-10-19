#!/usr/bin/env Rscript
# ==============================================================================
# Writes the combined object's cell metadata to CSV.
#
# The astrocyte Rmd's S12a chunk needs `scorch_meta` (the combined object's
# @meta.data) to compute each sample's astrocyte fraction, which requires all
# cell types, not just the astrocyte subset. Dumping the metadata once avoids
# loading the 9 GB object just to read it. The donor-1040-filtered object is
# preferred (VST cells of donor 1040 are excluded from the VST analyses).
#
# INPUT   OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3_no1040.RDS
#         (or SCORCH_COMBINED_OBJ; falls back to seurat_combined_v3.RDS with a warning)
# OUTPUT  OUT_ROOT/01-Regional_data_processing/combined/data/combined_v3_cell_metadata.csv
#
# Usage: revision/docker/run.sh revision/tools/dump_combined_metadata.R
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
rerun_setup("01-Regional_data_processing", subdirs = "combined/data")

## Prefer the 1040-filtered combined object. This CSV is what S12a counts cells
##   from (revision/tools/astro_manuscript_figures.R), and a dump from the
##   unfiltered object silently reintroduces donor 1040's 205 VST cells into a
##   figure labelled as excluding them. The fallback is announced, not silent.
.cands <- c(
  Sys.getenv("SCORCH_COMBINED_OBJ", unset = NA_character_),
  "./combined/data/seurat_combined_v3_no1040.RDS",
  "./combined/data/seurat_combined_v3.RDS"
)
IN <- .cands[!is.na(.cands) & file.exists(.cands)][1]
if (is.na(IN)) stop("no combined object found")
if (!grepl("no1040", IN))
  warning("[dump] ", basename(IN), " still contains donor 1040's VST cells",
          immediate. = TRUE)
OUT <- "./combined/data/combined_v3_cell_metadata.csv"
message("[dump] loading ", IN, " ...")
obj <- readRDS(IN)
md <- obj@meta.data

## Fail loudly if 1040's VST cells survived - this file feeds a figure that
## claims to exclude them, and a wrong row count here is invisible downstream.
.n1040 <- sum(as.character(md$sample_ids) == "1040" &
              as.character(md$region) == "VST", na.rm = TRUE)
message("[dump] donor 1040 VST cells present: ", .n1040)
if (.n1040 > 0 && grepl("no1040", IN))
  stop("filtered object still has ", .n1040, " of donor 1040's VST cells")
keep <- intersect(c("sample_ids","sample_name","condition_new","region",
                    "consensus_annotation_v2","consensus_annotation_v3",
                    "consensus_annotation_v3_relabeled","consensus_cell_type_coarse_l1",
                    "nCount_RNA","nFeature_RNA","percent.mt"), colnames(md))
message("[dump] ", nrow(md), " cells x ", length(keep), " columns")
write.csv(md[, keep], OUT)
message("[dump] wrote ", OUT)
