# Writes the donor-aware cross-region DE in the edgeR results-table layout, so every cross-region figure can be drawn from it.
#
# Method: figures read the table chosen by SCORCH_XR_DE (xr_de_file() in
#   revision/config/paths.R); SCORCH_XR_DE=donor_aware selects this one.
#   Statistics (gene, cluster_id, logFC, p_val, p_adj.loc, p_adj.glb) come from
#   the limma-voom + duplicateCorrelation(block = donor) fit of
#   revision/tools/cross_region_donor_aware.R (same design and contrast as the
#   edgeR model). The detection frequencies (avg_CTR_frq, avg_HIV_frq, *.frq) and
#   logCPM are joined from the edgeR table by gene x cell type: they describe
#   expression, not the test, and are identical under either model because the
#   pseudobulk counts are the same, so the > 5% expression filter selects the same
#   genes. logCPM is used only for display ordering.
#   Donor-aware rows with no edgeR counterpart (cell types pbDS dropped for rank
#   deficiency, e.g. PC, Lamp5, Lamp5 Lhx6, or genes edgeR's filter removed) have
#   no frequencies and are dropped and reported by cell type, so the figures cover
#   the same cell types as the edgeR versions.
#
# Inputs:  OUT_ROOT/03-DE_analysis/cross_region/donor_aware/cross_region_donor_aware_DE.csv
#          OUT_ROOT/03-DE_analysis/cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS
# Output:  OUT_ROOT/03-DE_analysis/cross_region/data/all_HIV-CTR_res_donorAware_withFreq.RDS
#
# Usage: revision/docker/run.sh revision/tools/cross_region_donor_aware_table.R

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages(library(dplyr))
D <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region")

da <- utils::read.csv(file.path(D, "donor_aware", "cross_region_donor_aware_DE.csv"), stringsAsFactors = FALSE)
ed <- as.data.frame(readRDS(file.path(D, "data", "all_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))) %>%
  mutate(cluster_id = as.character(cluster_id))
expr_cols <- c("logCPM", grep("frq$", colnames(ed), value = TRUE))
message("[xr-da] donor-aware rows: ", nrow(da), " | edgeR rows: ", nrow(ed),
        " | expression columns carried over: ", length(expr_cols))

out <- da %>% select(gene, cluster_id, logFC, p_val, p_adj.loc, p_adj.glb) %>%
  inner_join(ed %>% select(gene, cluster_id, all_of(expr_cols)), by = c("gene", "cluster_id")) %>%
  mutate(test = "limma_voom_duplicateCorrelation", contrast = unique(ed$contrast)[1])

drop <- da %>% anti_join(ed, by = c("gene", "cluster_id")) %>% count(cluster_id, name = "rows_without_edgeR_counterpart")
if (nrow(drop)) { message("[xr-da] dropped (no expression columns):"); print(as.data.frame(drop), row.names = FALSE) }
miss <- ed %>% anti_join(da, by = c("gene", "cluster_id")) %>% count(cluster_id, name = "edgeR_rows_not_in_donor_aware")
if (nrow(miss)) { message("[xr-da] edgeR rows absent from the donor-aware fit:"); print(as.data.frame(miss), row.names = FALSE) }

## p_adj.glb is kept from the donor-aware fit (adjusted over its full table, which
## includes the dropped rows), so it matches the values reported elsewhere.
stopifnot(!anyNA(out$avg_CTR_frq), !anyNA(out$logFC))
saveRDS(out, file.path(D, "data", "all_HIV-CTR_res_donorAware_withFreq.RDS"))

rule <- function(x) (x$avg_CTR_frq > 0.05 | x$avg_HIV_frq > 0.05) & abs(x$logFC) > 0.25 & x$p_adj.glb < 0.05
s <- bind_rows(
  ed %>% mutate(model = "edgeR", deg = rule(.)) %>% group_by(model, cluster_id) %>% summarise(DEGs = sum(deg), up = sum(deg & logFC > 0), .groups = "drop"),
  out %>% mutate(model = "donor_aware", deg = rule(.)) %>% group_by(model, cluster_id) %>% summarise(DEGs = sum(deg), up = sum(deg & logFC > 0), .groups = "drop")) %>%
  tidyr::pivot_wider(names_from = model, values_from = c(DEGs, up))
print(as.data.frame(s), row.names = FALSE)
message("[xr-da] wrote ", file.path(D, "data", "all_HIV-CTR_res_donorAware_withFreq.RDS"), " (", nrow(out), " rows)")
