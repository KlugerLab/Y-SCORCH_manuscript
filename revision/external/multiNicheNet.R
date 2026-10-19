# ==============================================================================
# Donor-aware cell-cell communication with MultiNicheNet (analysis by Junchen).
# Pseudobulks each cell type per donor and tests PWH vs PWoH with edgeR through
# muscat, per region, so interactions are compared across donors rather than across
# pooled cells. Replaces the per-cell-pair NICHES analysis of the original submission.
# Outputs: OUT_ROOT/08-communication/MultiNicheNet/output/
# Outlier-donor sensitivity: revision/tools/multinichenet_sensitivity.R
# Usage: revision/docker/run.sh revision/external/multiNicheNet.R
#        (needs the packages installed by revision/tools/install_ccc_rlibs.R)
# ==============================================================================
## ============================================================
## Donor-aware cell-cell communication with MultiNicheNet
## Replaces the NICHES-based (per-cell-pair, pseudoreplicated) pipeline.
## MultiNicheNet pseudobulks per DONOR per cell type, then runs
## edgeR-based DE at the sample level (via muscat), so HIV vs CTR
## differences are tested across donors, not across pooled cells.
##
## Source vignette this is adapted from (real repo, verified):
## https://github.com/saeyslab/multinichenetr/blob/main/vignettes/pairwise_analysis_MISC.Rmd
## ============================================================

## ---- 0. Install (run once) ---------------------------------
## IMPORTANT: install from the official saeyslab repo, not third-party
## mirrors (e.g. some r-universe listings under other usernames are NOT
## the real package/maintainer).
# install.packages("devtools")
# devtools::install_github("saeyslab/nichenetr")
# devtools::install_github("saeyslab/multinichenetr")

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
MNN_OUT <- file.path(OUT_ROOT, "08-communication", "MultiNicheNet", "output")
dir.create(MNN_OUT, recursive = TRUE, showWarnings = FALSE)
rerun_setup("",
            subdirs = c("01-Regional_data_processing/combined/data",
                        "08-communication/MultiNicheNet/figures_new",
                        "08-communication/MultiNicheNet/output"))

## ---- 1. Load and prep the same Seurat subset you already built ----
seurat.combined <- readRDS("./01-Regional_data_processing/combined/data/seurat_combined_v3_no1040_msn.RDS")

cell_types_oi <- c("L2/3 IT", "L4 IT", "Pvalb", "Micro", "Astro", "Endo")

seurat.combined.sub <- subset(
  seurat.combined,
  cells = colnames(seurat.combined)[seurat.combined@meta.data$consensus_annotation_v3 %in% cell_types_oi]
)
seurat.combined.sub <- JoinLayers(seurat.combined.sub)

## MultiNicheNet requires syntactically valid names for group/sample/
## celltype columns AND values (e.g. "L2/3 IT" is not valid -> make.names()).
seurat.combined.sub@meta.data$celltype_valid <- make.names(seurat.combined.sub@meta.data$consensus_annotation_v3)
seurat.combined.sub@meta.data$sample_valid   <- make.names(seurat.combined.sub@meta.data$sample_name)
seurat.combined.sub@meta.data$condition_valid <- make.names(as.character(seurat.combined.sub@meta.data$condition_new))
seurat.combined.sub@meta.data$region_valid    <- make.names(as.character(seurat.combined.sub@meta.data$region))

## Sanity check: each sample must map to exactly one condition
stopifnot(
  seurat.combined.sub@meta.data %>%
    dplyr::distinct(sample_valid, condition_valid) %>%
    dplyr::count(sample_valid) %>%
    dplyr::pull(n) %>%
    all(. == 1)
)

## ---- 2. Convert to SingleCellExperiment (raw counts required) ----
sce_full <- Seurat::as.SingleCellExperiment(seurat.combined.sub, assay = "RNA")
sce_full <- alias_to_symbol_SCE(sce_full, "human") %>% makenames_SCE()

## ---- 3. NicheNet v2 prior LR network + ligand-target matrix (human) ----
## These are downloaded once from Zenodo; cache locally if you'll rerun often.
options(timeout = 300)
lr_network_all <- readRDS(url(
  "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
)) %>%
  mutate(
    ligand = convert_alias_to_symbols(ligand, organism = "human"),
    receptor = convert_alias_to_symbols(receptor, organism = "human")
  ) %>%
  mutate(ligand = make.names(ligand), receptor = make.names(receptor))

lr_network <- lr_network_all %>% distinct(ligand, receptor)

ligand_target_matrix <- readRDS(url(
  "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
))
colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = "human") %>% make.names()
rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = "human") %>% make.names()

lr_network <- lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix <- ligand_target_matrix[, lr_network$ligand %>% unique()]

## ---- 4. Fixed analysis settings ----
sample_id   <- "sample_valid"
group_id    <- "condition_valid"     # "CTR" / "HIV"
celltype_id <- "celltype_valid"

covariates <- NA
batches    <- NA   # see note at bottom re: adding region as a covariate instead of looping

senders_oi   <- make.names(cell_types_oi)
receivers_oi <- make.names(cell_types_oi)

contrasts_oi <- c("'HIV-CTR','CTR-HIV'")
contrast_tbl <- tibble(
  contrast = c("HIV-CTR", "CTR-HIV"),
  group    = c("HIV", "CTR")
)

## Parameters (same spirit as your NICHES min.cells.per.ident = 10 filter)
min_cells        <- 10
min_sample_prop  <- 0.50
fraction_cutoff  <- 0.05
empirical_pval   <- FALSE
logFC_threshold  <- 0.50
p_val_threshold  <- 0.05
p_val_adj        <- FALSE
top_n_target     <- 250
scenario         <- "regular"
ligand_activity_down <- FALSE
n.cores <- 4

## ---- 5. Run MultiNicheNet separately per region ----
## (Mirrors your original per-vec/per-region NICHES plots. If you'd rather
## pool donors across regions for more power, instead of looping you can
## run once on sce_full with covariates = "region_valid" -- but only do
## this if every region has samples in BOTH conditions, which yours does.)

regions <- unique(seurat.combined.sub@meta.data$region_valid)
multinichenet_output_list <- list()

for (region_id in regions) {
  
  message("Running MultiNicheNet for region: ", region_id)
  
  sce <- sce_full[, SummarizedExperiment::colData(sce_full)$region_valid == region_id]
  sce <- sce[, SummarizedExperiment::colData(sce)[, celltype_id] %in% c(senders_oi, receivers_oi)]
  sce <- sce[, SummarizedExperiment::colData(sce)[, group_id] %in% c("HIV", "CTR")]
  
  out <- tryCatch({
    multi_nichenet_analysis(
      sce = sce,
      celltype_id = celltype_id, sample_id = sample_id, group_id = group_id,
      batches = batches, covariates = covariates,
      lr_network = lr_network, ligand_target_matrix = ligand_target_matrix,
      contrasts_oi = contrasts_oi, contrast_tbl = contrast_tbl,
      senders_oi = senders_oi, receivers_oi = receivers_oi,
      min_cells = min_cells,
      fraction_cutoff = fraction_cutoff,
      min_sample_prop = min_sample_prop,
      scenario = scenario,
      ligand_activity_down = ligand_activity_down,
      logFC_threshold = logFC_threshold,
      p_val_threshold = p_val_threshold,
      p_val_adj = p_val_adj,
      empirical_pval = empirical_pval,
      top_n_target = top_n_target,
      n.cores = n.cores,
      verbose = TRUE
    )
  }, error = function(e) {
    message("  -> skipped ", region_id, ": ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(out)) {
    multinichenet_output_list[[region_id]] <- out
    saveRDS(out, file = file.path(MNN_OUT, paste0("multinichenet_", region_id, ".RDS")))
  }
}

## ---- 6. Top prioritized interactions + plots, per region per condition ----
for (region_id in names(multinichenet_output_list)) {
  
  mn_out <- multinichenet_output_list[[region_id]]
  
  for (group_oi in c("HIV", "CTR")) {
    
    top_lr <- get_top_n_lr_pairs(
      mn_out$prioritization_tables,
      top_n = 50,
      groups_oi = group_oi
    )
    
    if (nrow(top_lr) == 0) next
    
    plot_oi <- make_sample_lr_prod_activity_plots(
      mn_out$prioritization_tables,
      top_lr
    )
    
    ggsave(
      filename = paste0(MNN_OUT, "/",
                        region_id, "_", group_oi, "_top50.png"),
      plot = plot_oi, width = 16, height = 13, dpi = 300, limitsize = FALSE
    )
  }
}




## ---- 6. Top 50 prioritized interactions for specific sender-receiver
##         directions of interest, per region ----
## Cell type names here must match the make.names()-valid labels used above:
## "Micro", "Astro", "L4.IT" ("L4 IT" -> make.names() -> "L4.IT").
directions_oi <- list(
  list(sender = "Micro", receiver = "Astro"),
  list(sender = "Astro", receiver = "Micro")
)

for (region_id in names(multinichenet_output_list)) {
  
  mn_out <- multinichenet_output_list[[region_id]]
  
  for (dir_oi in directions_oi) {
    
    ## rank_per_group = FALSE -> top 20 for this sender->receiver direction
    ## across both HIV and CTR jointly (matches the "overview" ranking used
    ## for the circos plots in the vignette); each interaction still carries
    ## which group it's specific to, and the plot below shows both groups.
    top_lr <- get_top_n_lr_pairs(
      mn_out$prioritization_tables,
      top_n = 20,
      groups_oi = "HIV",     # restrict to interactions prioritized as up in HIV
      senders_oi = dir_oi$sender,
      receivers_oi = dir_oi$receiver
    )
    
    if (nrow(top_lr) == 0) {
      message("  -> no interactions found for ", region_id, ": ",
              dir_oi$sender, " -> ", dir_oi$receiver)
      next
    }
    
    ## make_sample_lr_prod_activity_plots() returns a 4-panel patchwork
    ## (LR product dotplot, scaled ligand activity, ligand activity,
    ## cell-type specificity). make_sample_lr_prod_plots() draws only the
    ## first panel -- the LR expression product dotplot -- on its own.
    plot_oi <- make_sample_lr_prod_plots(
      mn_out$prioritization_tables,
      top_lr
    ) &
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 15),
            strip.text.x = element_text(size = 15, face = "bold"),
            strip.text.y = element_text(size = 15, face = "bold", angle = 0),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15))
    
    
    ggsave(
      filename = paste0(MNN_OUT, "/",
                        region_id, "_", dir_oi$sender, "_to_", dir_oi$receiver,
                        "_top20.png"),
      plot = plot_oi, width = 12, height = 8, dpi = 300, limitsize = FALSE
    )
  }
}
