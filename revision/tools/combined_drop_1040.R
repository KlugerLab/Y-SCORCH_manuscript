# Removes VST sample 1040 from the combined object and regenerates the cell-counts figure.
#
# Method: after ambient RNA correction, VST keeps only samples with low
#   excitatory-neuron content (dissection contamination). 1040 is 38.1%
#   excitatory, the second-highest in VST, and is also the lowest-yield library
#   in the region (378 nuclei regionally, 205 in the combined object; next lowest
#   is 1,285), so it fails on contamination and on quality independently. The
#   published exclusion list did not remove 1040; the four samples it did remove
#   (1039, 1068, HCCRS, HCTXJ) are already absent from the combined object.
#   Only the VST library is removed; the donor's INS and PFC libraries stay.
#
#   The object is subset rather than rebuilt: rebuilding means rerunning
#   integration, refinement and FIt-SNE on ~642k cells, and 1040 contributes 205
#   of 641,923 cells (0.03%), so its influence on the Harmony integration and
#   clustering is negligible. Per-sample figures change materially, and
#   subsetting fixes those exactly. The embedding was therefore computed with
#   these cells present and they were removed afterwards (to be stated in the
#   methods). If the combined object is rebuilt, apply the exclusion at the VST
#   input instead.
#
#   1181 (multi-pathology) is a sensitivity exclusion, not a QC rule (see
#   revision/tools/donor_lodo_astro.R and the noHIVE_noPath DE runs), so main
#   figures keep it. --drop-1181 also produces the sensitivity variant.
#
# Input:   OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3.RDS
# Outputs: OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3_<SFX>.RDS,
#          n_counts_<SFX>.csv, n_counts_by_sample_<SFX>.csv;
#          combined/figures/n_counts_new_<SFX>.png  (SFX = no1040 or no1040_no1181)
#
# Usage: revision/docker/run.sh revision/tools/combined_drop_1040.R [--drop-1181]

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(cowplot)
})

args      <- commandArgs(trailingOnly = TRUE)
ALSO_1181 <- "--drop-1181" %in% args

rerun_setup("01-Regional_data_processing",
            subdirs = c("combined/data", "combined/figures"))

IN  <- "./combined/data/seurat_combined_v3.RDS"
if (!file.exists(IN)) stop("missing ", IN)

DROP <- c("1040")
SFX  <- "no1040"
if (ALSO_1181) { DROP <- c(DROP, "1181"); SFX <- "no1040_no1181" }

## The exclusion is VST-only. Donors contribute a library to each region, and
## sample_ids repeats across regions - 1040 has 205 cells in VST but ~13,770 in
## INS and PFC. The rule is about VST dissection contamination, so only the VST
## library is removed; the same donor's INS and PFC libraries are unaffected and
## must stay. Matching on sample_ids alone silently deletes all three.
DROP_REGION <- "VST"

message("[drop] loading ", IN, " ...")
obj <- readRDS(IN)
md  <- obj@meta.data
md$sample_ids <- as.character(md$sample_ids)
md$region     <- as.character(md$region)

n_before <- ncol(obj)
is_target <- md$sample_ids %in% DROP & md$region == DROP_REGION
present   <- intersect(DROP, unique(md$sample_ids[md$region == DROP_REGION]))
if (!length(present)) {
  message("[drop] none of ", paste(DROP, collapse = ", "), " present in ",
          DROP_REGION, " - nothing to do")
} else {
  for (d in present) {
    n_vst   <- sum(md$sample_ids == d & md$region == DROP_REGION)
    n_other <- sum(md$sample_ids == d & md$region != DROP_REGION)
    message("[drop] ", d, ": removing ", n_vst, " ", DROP_REGION, " cells; ",
            "KEEPING ", n_other, " cells from other regions")
  }
}

keep_cells <- colnames(obj)[!is_target]
obj <- subset(obj, cells = keep_cells)

message("[drop] ", n_before, " -> ", ncol(obj), " cells (removed ",
        n_before - ncol(obj), ")")

.md <- obj@meta.data
for (rg in c("VST", "INS", "PFC")) {
  sel <- as.character(.md$region) == rg
  message("[drop] ", rg, ": ", length(unique(as.character(.md$sample_ids)[sel])),
          " samples, ", sum(sel), " cells")
}
## guard: the donor must survive in the regions the rule does not touch
for (d in DROP) {
  n_left <- sum(as.character(.md$sample_ids) == d)
  if (n_left == 0)
    warning("donor ", d, " was removed from EVERY region - the exclusion is ",
            "meant to be ", DROP_REGION, "-only")
  else
    message("[drop] donor ", d, " retains ", n_left, " cells outside ", DROP_REGION)
}

OUT <- sprintf("./combined/data/seurat_combined_v3_%s.RDS", SFX)
saveRDS(obj, OUT)
message("[drop] wrote ", OUT)

## ---- regenerate the cell-counts figure ---------------------------------------
## Same aesthetics as the published cell-counts figure in
## 01-Regional_data_processing/combined/03-integration_tsne_refined.R, so the
## only difference is the removed sample. Labels stay
## CTR/HIV to match; pwh_relabel() is available if the figure is ever converted
## to person-first language.
metadata_all <- obj@meta.data

n_counts <- metadata_all %>%
  group_by(region, condition_new) %>%
  summarise(n = n(), .groups = "drop")

n_counts$condition_new <- factor(n_counts$condition_new, levels = c("HIV", "CTR"))
n_counts$region        <- factor(n_counts$region, levels = c("VST", "INS", "PFC"))

cond_col <- c("#0072B2", "darkorange1")
names(cond_col) <- c("CTR", "HIV")

p1 <- ggplot(n_counts, aes(x = region, y = n, fill = condition_new)) +
  scale_fill_manual(values = cond_col) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = n), hjust = -0.05,
            position = position_dodge(width = 1), size = 20) +
  ylab("# Cells") + xlab("") +
  theme_cowplot(font_size = 30) +
  ylim(c(0, 180000)) +
  theme(axis.text = element_text(size = 60),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        title = element_text(size = 60),
        legend.text = element_text(size = 10)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + NoLegend() + coord_flip()

FIG <- sprintf("./combined/figures/n_counts_new_%s.png", SFX)
png(file = FIG, width = 4500, height = 4500, res = 300)
print(p1)
dev.off()
message("[drop] wrote ", FIG)

## the numbers behind the bars, so the figure can be checked without reopening it
utils::write.csv(n_counts,
  sprintf("./combined/data/n_counts_%s.csv", SFX), row.names = FALSE)

n_count_sample <- metadata_all %>%
  group_by(sample_ids, region) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(region, n)
utils::write.csv(n_count_sample,
  sprintf("./combined/data/n_counts_by_sample_%s.csv", SFX), row.names = FALSE)

message("\n########## CELL COUNTS BY REGION x CONDITION ##########")
print(as.data.frame(n_counts), row.names = FALSE)
message("\nsmallest samples per region (check nothing else is marginal):")
print(as.data.frame(n_count_sample %>% group_by(region) %>% slice_head(n = 3)),
      row.names = FALSE)
message("\n[drop] done")
