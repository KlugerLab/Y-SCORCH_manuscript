# ==============================================================================
# Outlier-donor sensitivity of the MultiNicheNet cell-cell communication analysis.
#
# Method: MultiNicheNet settings follow Junchen's analysis
# (revision/external/multiNicheNet.R, not modified): same object, the same six
# cell types, per-region runs, covariates = NA, min_cells 10, min_sample_prop
# 0.5, fraction_cutoff 0.05, logFC_threshold 0.5, p_val_threshold 0.05
# (unadjusted), top_n_target 250, NicheNet v2 human prior (the same two Zenodo
# files, cached under OUT_ROOT/reference/nichenet/ so they are downloaded once).
#
# Two runs per region, both in this environment (RLIBS_CCC):
#   full           all donors - the reference, same package versions
#   noHIVE_noPath  without donors 1039 and 1181
# Both runs use the same package versions so a difference reflects the donors,
# not version drift. The full run is also compared with the original
# full-cohort result (top-50 overlap) so any version drift is reported.
#
# Comparison (group_prioritization_tbl, interactions prioritised as up in PWH):
#   overlap of the top 50 (and top 20 Micro<->Astro, the directions shown in
#   the original figure) between full and exclusion; Spearman of
#   prioritization_score over interactions scored in both; where each
#   full-cohort top-50 interaction ranks without the two donors.
#
# Inputs: OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3_no1040_msn.RDS;
#         MULTINICHENET_REF_DIR (original full-cohort results)
# Output: OUT_ROOT/08-communication/MultiNicheNet_sensitivity/
#           multinichenet_<R>_<run>.RDS, multinichenet_sensitivity_summary.csv,
#           multinichenet_top_interactions_<R>.csv
#
# Usage: revision/docker/run.sh revision/tools/multinichenet_sensitivity.R [--regions INS,PFC,VST]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
.libPaths(c(RLIBS_CCC, .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(dplyr)
  library(nichenetr); library(multinichenetr)
})
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
REGIONS <- strsplit(getarg("--regions", "INS,PFC,VST"), ",")[[1]]
PATHOLOGY_DONORS <- c("1039", "1181")
OUT <- file.path(OUT_ROOT, "08-communication", "MultiNicheNet_sensitivity")
REF <- file.path(OUT_ROOT, "reference", "nichenet")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(REF, recursive = TRUE, showWarnings = FALSE)
JYC <- MULTINICHENET_REF_DIR   # full-cohort results of revision/external/multiNicheNet.R

## ---- data, exactly as the original script ----------------------------------------------------
obj <- readRDS(file.path(OUT_ROOT, "01-Regional_data_processing", "combined", "data",
                         "seurat_combined_v3_no1040_msn.RDS"))
cell_types_oi <- c("L2/3 IT", "L4 IT", "Pvalb", "Micro", "Astro", "Endo")
sub <- subset(obj, cells = colnames(obj)[obj@meta.data$consensus_annotation_v3 %in% cell_types_oi])
rm(obj); invisible(gc())
sub <- JoinLayers(sub)
md <- sub@meta.data
sub$celltype_valid  <- make.names(md$consensus_annotation_v3)
sub$sample_valid    <- make.names(md$sample_name)
sub$condition_valid <- make.names(as.character(md$condition_new))
sub$region_valid    <- make.names(as.character(md$region))
sub$donor_id        <- as.character(md$sample_ids)
sce_full <- Seurat::as.SingleCellExperiment(sub, assay = "RNA")
rm(sub); invisible(gc())
sce_full <- alias_to_symbol_SCE(sce_full, "human") %>% makenames_SCE()

## ---- NicheNet v2 prior (the original two files), cached ---------------------------------------
fetch <- function(url, f) {
  p <- file.path(REF, f)
  if (!file.exists(p)) { options(timeout = 600); utils::download.file(url, p, mode = "wb") }
  readRDS(p)
}
lr_network <- fetch("https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds",
                    "lr_network_human_allInfo_30112033.rds") %>%
  mutate(ligand = convert_alias_to_symbols(ligand, organism = "human"),
         receptor = convert_alias_to_symbols(receptor, organism = "human")) %>%
  mutate(ligand = make.names(ligand), receptor = make.names(receptor)) %>%
  distinct(ligand, receptor)
ltm <- fetch("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds",
             "ligand_target_matrix_nsga2r_final.rds")
colnames(ltm) <- colnames(ltm) %>% convert_alias_to_symbols(organism = "human") %>% make.names()
rownames(ltm) <- rownames(ltm) %>% convert_alias_to_symbols(organism = "human") %>% make.names()
lr_network <- lr_network %>% filter(ligand %in% colnames(ltm))
ltm <- ltm[, unique(lr_network$ligand)]

run_mnn <- function(sce) multi_nichenet_analysis(
  sce = sce, celltype_id = "celltype_valid", sample_id = "sample_valid", group_id = "condition_valid",
  batches = NA, covariates = NA, lr_network = lr_network, ligand_target_matrix = ltm,
  contrasts_oi = c("'HIV-CTR','CTR-HIV'"),
  contrast_tbl = tibble::tibble(contrast = c("HIV-CTR", "CTR-HIV"), group = c("HIV", "CTR")),
  senders_oi = make.names(cell_types_oi), receivers_oi = make.names(cell_types_oi),
  min_cells = 10, fraction_cutoff = 0.05, min_sample_prop = 0.50, scenario = "regular",
  ligand_activity_down = FALSE, logFC_threshold = 0.50, p_val_threshold = 0.05, p_val_adj = FALSE,
  empirical_pval = FALSE, top_n_target = 250, n.cores = 4, verbose = FALSE)

prio_hiv <- function(o) o$prioritization_tables$group_prioritization_tbl %>%
  filter(group == "HIV") %>% arrange(desc(prioritization_score)) %>%
  distinct(id, .keep_all = TRUE) %>% mutate(rank = row_number())

summ <- list()
for (R in REGIONS) {
  message("\n===== ", R)
  sce <- sce_full[, sce_full$region_valid == R & sce_full$condition_valid %in% c("HIV", "CTR")]
  runs <- list(full = sce, noHIVE_noPath = sce[, !sce$donor_id %in% PATHOLOGY_DONORS])
  message(sprintf("[mnn] %s: %d donors full, %d without %s", R, length(unique(sce$donor_id)),
                  length(unique(runs$noHIVE_noPath$donor_id)),
                  paste(intersect(PATHOLOGY_DONORS, unique(sce$donor_id)), collapse = "+")))
  res <- list()
  for (k in names(runs)) {
    f <- file.path(OUT, sprintf("multinichenet_%s_%s.RDS", R, k))
    res[[k]] <- if (file.exists(f)) readRDS(f) else { o <- run_mnn(runs[[k]]); saveRDS(o, f); o }
  }
  a <- prio_hiv(res$full); b <- prio_hiv(res$noHIVE_noPath)
  j <- inner_join(a %>% select(id, sender, receiver, ligand, receptor, score_full = prioritization_score, rank_full = rank),
                  b %>% select(id, score_excl = prioritization_score, rank_excl = rank), by = "id")
  top <- function(x, n, s = NULL, r = NULL) {
    if (!is.null(s)) x <- x[x$sender == s & x$receiver == r, ]
    head(x$id, n)
  }
  ov <- function(n, s = NULL, r = NULL) { A <- top(a, n, s, r); length(intersect(A, top(b, n, s, r))) / max(1, length(A)) }
  jyc <- file.path(JYC, sprintf("multinichenet_%s.RDS", R))
  ov_jyc <- if (file.exists(jyc)) { J <- prio_hiv(readRDS(jyc)); length(intersect(head(J$id, 50), head(a$id, 50))) / 50 } else NA
  summ[[R]] <- data.frame(
    region = R, n_interactions_full = nrow(a), n_interactions_excl = nrow(b),
    top50_overlap = ov(50), top20_micro_to_astro = ov(20, "Micro", "Astro"),
    top20_astro_to_micro = ov(20, "Astro", "Micro"),
    spearman_score = suppressWarnings(cor(j$score_full, j$score_excl, method = "spearman")),
    top50_full_median_rank_excl = stats::median(j$rank_excl[j$rank_full <= 50]),
    top50_overlap_vs_junchen_full = ov_jyc, stringsAsFactors = FALSE)
  utils::write.csv(j %>% arrange(rank_full) %>% head(100),
                   file.path(OUT, sprintf("multinichenet_top_interactions_%s.csv", R)), row.names = FALSE)
  print(summ[[R]], row.names = FALSE, digits = 3)
}
s <- bind_rows(summ)
utils::write.csv(s, file.path(OUT, "multinichenet_sensitivity_summary.csv"), row.names = FALSE)
options(width = 200); print(s, row.names = FALSE, digits = 3)
message("[mnn] done")
