# ==============================================================================
# Supplementary DEG tables (DEG_{INS,PFC,VST,avr_all_regions}.csv) regenerated
# from the CellBender pipeline.
#
# Method
#   regional      edgeR pseudobulk DE, one sample per donor; all tested
#                 gene x cell-type rows, as in the published tables.
#   cross-region  the donor-aware primary model (limma-voom +
#                 duplicateCorrelation, block = donor), with the edgeR
#                 cross-region logFC, adjusted p, test and DEG call beside it
#                 (columns *_edgeR).
#   The published tables had no effect size; logFC and a DEG flag are added.
#
# Columns  gene, cell_type, logFC, logCPM, p_value, p_adj (global BH across genes
#          and cell types), pct_PWoH / pct_PWH (mean fraction of nuclei expressing
#          the gene across PWoH / PWH samples), region, DEG = the paper's rule
#          (expressed in > 5% of nuclei in either group, |logFC| > 0.25,
#          p_adj < 0.05).
#
# Inputs   OUT_ROOT/03-DE_analysis/<R>/results/<R>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#          OUT_ROOT/03-DE_analysis/cross_region/data/all_HIV-CTR_res_donorAware_withFreq.RDS
#          OUT_ROOT/03-DE_analysis/cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS
# Outputs  OUT_ROOT/tables/supplementary_tables/DEG_*.csv and DEG_tables_summary.csv
#          (revision/tools/supp_qc_table.py writes the .xlsx copies)
#
# Usage
#   revision/docker/run.sh revision/tools/supp_deg_tables.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(dplyr); library(data.table) })
DE  <- file.path(OUT_ROOT, "03-DE_analysis")
OUT <- file.path(OUT_ROOT, "tables", "supplementary_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

tidy <- function(x) x %>% ungroup() %>% as.data.frame() %>%
  transmute(gene, cell_type = as.character(cluster_id), logFC,
            logCPM = if ("logCPM" %in% names(.)) logCPM else NA_real_,
            p_value = p_val, p_adj = p_adj.glb,
            pct_PWoH = avg_CTR_frq, pct_PWH = avg_HIV_frq,
            DEG = (avg_CTR_frq > 0.05 | avg_HIV_frq > 0.05) & abs(logFC) > 0.25 & p_adj.glb < 0.05)

rep <- list()
for (R in c("INS", "PFC", "VST")) {
  x <- readRDS(file.path(DE, R, "results", sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq.RDS", R)))
  d <- tidy(x) %>% mutate(region = R) %>% arrange(cell_type, gene)
  fwrite(d, file.path(OUT, sprintf("DEG_%s.csv", R)))
  rep[[R]] <- data.frame(table = paste0("DEG_", R), rows = nrow(d), cell_types = n_distinct(d$cell_type),
                         DEGs = sum(d$DEG), up = sum(d$DEG & d$logFC > 0), down = sum(d$DEG & d$logFC < 0))
}

da <- readRDS(file.path(DE, "cross_region", "data", "all_HIV-CTR_res_donorAware_withFreq.RDS"))
ed <- as.data.frame(readRDS(file.path(DE, "cross_region", "data", "all_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))) %>%
  transmute(gene, cell_type = as.character(cluster_id), logFC_edgeR = logFC, p_adj_edgeR = p_adj.glb,
            test_edgeR = if ("test" %in% names(.)) test else "QL")
d <- tidy(da) %>% select(-logCPM) %>%
  left_join(ed, by = c("gene", "cell_type")) %>%
  mutate(DEG_edgeR = (pct_PWoH > 0.05 | pct_PWH > 0.05) & abs(logFC_edgeR) > 0.25 & p_adj_edgeR < 0.05) %>%
  arrange(cell_type, gene)
fwrite(d, file.path(OUT, "DEG_avr_all_regions.csv"))
rep$xr <- data.frame(table = "DEG_avr_all_regions (donor-aware)", rows = nrow(d), cell_types = n_distinct(d$cell_type),
                     DEGs = sum(d$DEG), up = sum(d$DEG & d$logFC > 0), down = sum(d$DEG & d$logFC < 0))
s <- bind_rows(rep)
fwrite(s, file.path(OUT, "DEG_tables_summary.csv"))
print(s, row.names = FALSE)
message("[supp] DEG tables written to ", OUT)
