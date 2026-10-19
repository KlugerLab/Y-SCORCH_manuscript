# ==============================================================================
# Supplementary violin panels, PWH vs PWoH split by region: Figure S6
# (endothelial cytokine receptors) and, with --panel S10, Figure S10 (microglial
# LDAM markers).
#
# Same call as the published figure
# (03-DE_analysis/validation_scripts/01-validate_de.R:29; S10 uses :38):
#   VlnPlot(features = IL4R, IL6R, IL18R1, BCL6, idents = "Endo",
#           group.by = condition, split.by = region, pt.size = 0.2, alpha = 0.3,
#           cols = region_cols, ncol = 4) & xlab("Condition")
# on the revision's combined object (donor 1040's VST cells excluded), with two
# deliberate differences:
#   * conditions are labelled PWoH / PWH, as in every other revision figure
#   * region colours keep the published S6 mapping (INS orange, PFC blue, VST
#     pink). NB the tSNE/UMAP figures use the OTHER mapping from
#     03-integration_tsne_refined.R (PFC orange, INS blue) - the two published
#     scripts disagree; --region-cols tsne switches to that one.
#
# INPUT   OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3_no1040_msn.RDS
# OUTPUT  OUT_ROOT/03-DE_analysis/figures/s6/   (s10/ for --panel S10)
#           figS6_endo_cytokine_receptors.png/.pdf
#           figS6_endo_cytokine_receptors_cells.csv   nuclei per region x condition behind the violins
#
# Usage: revision/docker/run.sh revision/tools/figS_violin_panels.R
#          [--panel S6|S10] [--region-cols published|tsne]
# ==============================================================================

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
COLS <- match.arg(getarg("--region-cols", "published"), c("published", "tsne"))
## S6 = endothelial cytokine receptors, S10 = microglial LDAM markers
## (01-validate_de.R:29 and :38 - same call, different idents / features / points)
PANEL <- match.arg(getarg("--panel", "S6"), c("S6", "S10"))
if (PANEL == "S6") {
  IDENT <- "Endo"; FEATURES <- c("IL4R", "IL6R", "IL18R1", "BCL6")
  PT <- 0.2; ALPHA <- 0.3; STEM <- "figS6_endo_cytokine_receptors"
} else {
  IDENT <- "Micro"; FEATURES <- c("ACSL1", "NAMPT", "DPYD", "CD163")
  PT <- 0.05; ALPHA <- 0.2; STEM <- "figS10_micro_ldam_markers"
}
OUT <- file.path(OUT_ROOT, "03-DE_analysis", "figures", tolower(PANEL))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

region_cols <- if (COLS == "published") {
  c(INS = "#E69F00", PFC = "#56B4E9", VST = "#CC79A7")
} else {
  c(INS = "#56B4E9", PFC = "#E69F00", VST = "#CC79A7")
}

f <- file.path(OUT_ROOT, "01-Regional_data_processing", "combined", "data",
               "seurat_combined_v3_no1040_msn.RDS")
message("[figS] loading ", basename(f))
obj <- readRDS(f)
Idents(obj) <- "consensus_annotation_v3"
obj$condition_new <- factor(as.character(pwh_relabel(obj$condition_new)), levels = PWH_LEVELS)
obj$region <- factor(as.character(obj$region), levels = c("INS", "PFC", "VST"))

n <- obj@meta.data %>% dplyr::filter(consensus_annotation_v3 == IDENT) %>%
  dplyr::count(region, condition_new)
utils::write.csv(n, file.path(OUT, paste0(STEM, "_cells.csv")), row.names = FALSE)
message("[figS] ", IDENT, " nuclei: ",
        paste(sprintf("%s %s=%d", n$region, n$condition_new, n$n), collapse = ", "))
miss <- setdiff(FEATURES, rownames(obj))
if (length(miss)) stop("features not in the object: ", paste(miss, collapse = ", "))

p <- VlnPlot(obj, features = FEATURES, idents = IDENT, group.by = "condition_new",
             split.by = "region", pt.size = PT, alpha = ALPHA, cols = region_cols, ncol = 4) &
  xlab("Condition")
png(file.path(OUT, paste0(STEM, ".png")), width = 3000, height = 1000, res = 300)
print(p); dev.off()
grDevices::cairo_pdf(file.path(OUT, paste0(STEM, ".pdf")), width = 10, height = 3.33)
print(p); dev.off()
message("[figS] done -> ", OUT)
