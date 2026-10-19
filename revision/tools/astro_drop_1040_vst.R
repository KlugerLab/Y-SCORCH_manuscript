# ==============================================================================
# Remove donor 1040's VST astrocytes from the labelled (pooled-region) astrocyte object.
#
# Method
#   After ambient RNA correction, VST keeps only samples with low excitatory-neuron
#   content. Donor 1040 scores 0.381 and is the smallest VST library, so its VST
#   cells are excluded everywhere (see VST_SAMPLE_KEEP_OVERRIDE in
#   revision/config/paths.R).
#   The exclusion is VST-only. The astrocyte object pools all three regions, and
#   1040 contributes astrocytes to each (INS 799, PFC 617, VST 33). Only the 33 VST
#   cells are removed; matching on sample_ids alone would also drop the donor's
#   valid INS and PFC data. A guard stops if the donor disappears entirely.
#   Expected impact is small: 1040 keeps 1,416 of 1,449 astrocytes (97.7%) and
#   remains a donor in the pooled analysis. (In the per-region microglia analysis,
#   by contrast, removing 1040 takes VST from 27 donors to 26.)
#
# Inputs   OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled.RDS
# Outputs  OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled_no1040vst.RDS
#
# Usage: revision/docker/run.sh revision/tools/astro_drop_1040_vst.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages(library(Seurat))

rerun_setup("06-astrocyte_analysis", subdirs = c("data"))

IN  <- "./data/astro_all_clusters_labeled.RDS"
OUT <- "./data/astro_all_clusters_labeled_no1040vst.RDS"
if (!file.exists(IN)) stop("missing ", IN)

DROP_DONOR  <- "1040"
DROP_REGION <- "VST"

message("[astro-drop] loading ", IN, " ...")
obj <- readRDS(IN)
md  <- obj@meta.data
md$sample_ids <- as.character(md$sample_ids)
md$region     <- as.character(md$region)

before <- table(md$region[md$sample_ids == DROP_DONOR])
message("[astro-drop] ", DROP_DONOR, " astrocytes by region before: ",
        paste(sprintf("%s=%d", names(before), before), collapse = "  "))

is_target <- md$sample_ids == DROP_DONOR & md$region == DROP_REGION
message("[astro-drop] removing ", sum(is_target), " ", DROP_REGION,
        " astrocytes from ", DROP_DONOR)

n_before <- ncol(obj)
obj <- subset(obj, cells = colnames(obj)[!is_target])

md2 <- obj@meta.data
md2$sample_ids <- as.character(md2$sample_ids)
after <- table(as.character(md2$region)[md2$sample_ids == DROP_DONOR])
message("[astro-drop] ", DROP_DONOR, " astrocytes by region after : ",
        paste(sprintf("%s=%d", names(after), after), collapse = "  "))

## guard: the donor must survive - this is a VST-only exclusion
if (sum(md2$sample_ids == DROP_DONOR) == 0)
  stop("donor ", DROP_DONOR, " was removed entirely - the exclusion is ",
       DROP_REGION, "-only and this is a bug")

message("[astro-drop] ", n_before, " -> ", ncol(obj), " cells")
message("[astro-drop] donors: ", length(unique(md2$sample_ids)))
cat("\ncells by region:\n"); print(table(as.character(md2$region)))

saveRDS(obj, OUT)
message("\n[astro-drop] wrote ", OUT)
message("[astro-drop] run the figures with:")
message("  revision/docker/run.sh revision/tools/run_astro_reviewer_figures.R --in ", basename(OUT))
