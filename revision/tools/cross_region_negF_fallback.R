# Builds the cross-region edgeR table with the negative-F fallback, from the saved pseudobulk and design.
#
# Method: in the cross-region design, edgeR 4's default glmQLFit returns
#   negative F statistics for some genes (e.g. SERPINA3 in astrocytes). Those
#   rows are refit; all others keep the default QL result. The 90-sample
#   pseudobulk (the ~1.5 h part of 03-DE_analysis/07-cross_region_edgeR.R) is
#   not rebuilt:
#   1. pbDS on the saved pseudo_bulk.rds / mm.rds with edgeR 4 default QL,
#      checked to reproduce the saved default table (*_edgeR4defaultQL) exactly
#      at the 3 significant digits resDS stores
#   2. pbds_negF_fallback() (revision/config/paths.R): only rows with F < 0 are
#      refit, with legacy QL, else LRT; the `test` column records which
#   3. statistics swapped into the saved default table, so the frequency columns
#      (which need the 9 GB object to recompute) are carried over unchanged
#   07-cross_region_edgeR.R applies the same fallback itself on a full run.
#
# Inputs:  OUT_ROOT/03-DE_analysis/cross_region/data/{pseudo_bulk.rds, mm.rds, mm_edgeR4defaultQL.rds,
#          all_HIV-CTR_res_edgeR_filtered_withFreq_edgeR4defaultQL.RDS, ..._legacyQL.RDS}
# Outputs: OUT_ROOT/03-DE_analysis/cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#          (the file 07 writes; never overwritten if present), negF_fallback_rows.csv
#
# Usage: revision/docker/run.sh revision/tools/cross_region_negF_fallback.R

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(muscat); library(limma); library(dplyr)
})

D        <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "data")
DEF_FILE <- file.path(D, "all_HIV-CTR_res_edgeR_filtered_withFreq_edgeR4defaultQL.RDS")
LEG_FILE <- file.path(D, "all_HIV-CTR_res_edgeR_filtered_withFreq_legacyQL.RDS")
OUT_FILE <- file.path(D, "all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
CONTRAST <- "(INS.HIV-INS.CTR)/3+(PFC.HIV-PFC.CTR)/3+(VST.HIV-VST.CTR)/3"
STAT     <- c("logFC", "logCPM", "F", "p_val", "p_adj.loc", "p_adj.glb")

if (file.exists(OUT_FILE)) stop("primary table already exists, not overwriting: ", OUT_FILE)

pb <- readRDS(file.path(D, "pseudo_bulk.rds"))
mm <- readRDS(file.path(D, "mm.rds"))
stopifnot(isTRUE(all.equal(mm, readRDS(file.path(D, "mm_edgeR4defaultQL.rds")))))
contrast <- limma::makeContrasts(contrasts = CONTRAST, levels = mm)
message("[negF] ", ncol(pb), " samples, ", length(assayNames(pb)), " cell types")

edger_ql_legacy(FALSE)   # make sure muscat's glmQLFit is the edgeR default
res <- pbDS(pb, design = mm, contrast = contrast, filter = "both", min_cells = 10,
            method = "edgeR", verbose = FALSE)
flat <- function(r) dplyr::bind_rows(lapply(r$table, dplyr::bind_rows)) %>%
  mutate(cluster_id = as.character(cluster_id))

## ---- 1. reproduces the saved default run? -------------------------------------
def_saved <- as.data.frame(readRDS(DEF_FILE)) %>% mutate(cluster_id = as.character(cluster_id))
chk <- inner_join(def_saved %>% select(gene, cluster_id, all_of(STAT)),
                  flat(res) %>% select(gene, cluster_id, all_of(STAT)) %>%
                    mutate(across(all_of(STAT), ~ signif(.x, 3))),
                  by = c("gene", "cluster_id"), suffix = c(".saved", ".new"))
mism <- vapply(STAT, function(s) sum(abs(chk[[paste0(s, ".saved")]] - chk[[paste0(s, ".new")]]) >
                                      1e-12 * pmax(1, abs(chk[[paste0(s, ".saved")]]))), numeric(1))
message("[negF] default refit vs saved default table: ", nrow(chk), " of ", nrow(def_saved),
        " rows matched; mismatches per column: ", paste(names(mism), mism, collapse = ", "))
if (nrow(chk) != nrow(def_saved) || any(mism > 0))
  stop("pbDS on the saved pseudobulk does not reproduce the saved default table")

## ---- 2. fallback --------------------------------------------------------------
hyb <- flat(pbds_negF_fallback(res, design = mm, contrast = contrast)) %>%
  select(gene, cluster_id, all_of(STAT), test) %>%
  mutate(across(all_of(STAT), ~ signif(.x, 3)))

## ---- 3. swap statistics into the saved table ----------------------------------
out <- def_saved %>% select(-all_of(STAT)) %>% left_join(hyb, by = c("gene", "cluster_id"))
out <- out[, c(names(def_saved), "test")]
stopifnot(nrow(out) == nrow(def_saved), !anyNA(out$p_val), !anyNA(out$test))
saveRDS(out, OUT_FILE)
message("[negF] wrote ", OUT_FILE)

## ---- report -------------------------------------------------------------------
fb <- out %>% filter(test != "QL") %>%
  left_join(def_saved %>% select(gene, cluster_id, F_default = F, p_default = p_val),
            by = c("gene", "cluster_id")) %>%
  select(cluster_id, gene, test, logFC, F, p_val, p_adj.glb, F_default, p_default,
         avg_CTR_frq, avg_HIV_frq) %>% arrange(cluster_id, p_val)
utils::write.csv(fb, file.path(D, "negF_fallback_rows.csv"), row.names = FALSE)
cat("\n=== tests used ===\n"); print(table(out$test))
cat("\n=== fallback rows ===\n"); print(fb, row.names = FALSE)

leg <- as.data.frame(readRDS(LEG_FILE)) %>% mutate(cluster_id = as.character(cluster_id))
n_sig <- function(x) x %>% group_by(cluster_id) %>% summarise(n = sum(p_adj.glb < 0.05), .groups = "drop")
cmp <- n_sig(def_saved) %>% rename(default_all = n) %>%
  full_join(n_sig(leg) %>% rename(legacy_all = n), by = "cluster_id") %>%
  full_join(n_sig(out) %>% rename(fallback_primary = n), by = "cluster_id")
cat("\n=== significant (p_adj.glb < 0.05) ===\n"); print(as.data.frame(cmp), row.names = FALSE)
cat("totals: default", sum(cmp$default_all), " legacy", sum(cmp$legacy_all),
    " fallback", sum(cmp$fallback_primary), "\n")
j <- inner_join(def_saved %>% select(gene, cluster_id, pd = p_adj.glb),
                out %>% select(gene, cluster_id, pf = p_adj.glb, test), by = c("gene", "cluster_id"))
cat("calls changed vs default: gained", sum(j$pf < 0.05 & j$pd >= 0.05),
    "(of which fallback rows", sum(j$pf < 0.05 & j$pd >= 0.05 & j$test != "QL"), ")",
    " lost", sum(j$pf >= 0.05 & j$pd < 0.05), "\n")
cat("\n=== reactive-astrocyte panel ===\n")
print(out %>% filter(cluster_id == "Astro", gene %in% c("SERPINA3", "OSMR", "CHI3L1", "CD44", "GFAP", "C3")) %>%
        select(gene, logFC, F, p_val, p_adj.glb, test), row.names = FALSE)
message("[negF] done")
