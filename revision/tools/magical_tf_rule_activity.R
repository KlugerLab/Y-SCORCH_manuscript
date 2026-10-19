# ==============================================================================
# Which transcription factors are drawn in the astrocyte MAGICAL network.
#
# A TF passes when, in the region's MAGICAL astrocyte circuits:
#   1. it is expressed in >= 5% of PWH astrocyte nuclei (cross-region table, avg_HIV_frq),
#   2. it binds peaks linked to >= 5 genes up-regulated in that region's astrocytes
#      (regional DE, logFC > 0 and FDR < 0.05), and
#   3. its motif activity is higher in PWH in the donor-aware chromVAR test across regions
#      (limma + duplicateCorrelation, pooled centred model; difference > 0, FDR < 0.05).
#      A TF without a motif of its own uses a composite motif naming it (e.g. Smad2::Smad3).
# Motif activity, not the TF's own RNA, is the evidence of increased activity: TF
# activity is often regulated post-transcriptionally, and the motif test is donor-aware.
# The TF's own cross-region RNA change is reported but not required.
#
# Of the TFs that pass, the network figure draws a representative subset chosen to cover
# the TF families (STAT, AP-1, NF-kB, C/EBP, BACH): --fig-tfs, default STAT1, STAT3,
# CEBPD, FOSL2, BACH1, REL. All passing TFs are in the output table.
#
# Inputs   OUT_ROOT/magical/analyzed_output/<R>_new_Astro/expanded_circuits.csv
#          OUT_ROOT/03-DE_analysis/<R>/results/<R>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#          OUT_ROOT/03-DE_analysis/cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#          OUT_ROOT/04-ATAC_DE/motif_chromvar/crossregion/Astro_chromvar_crossregion.csv
# Output   OUT_ROOT/magical/analyzed_output/<R>_new_Astro/tf_selection_rule_activity.csv
#
# Usage: revision/docker/run.sh revision/tools/magical_tf_rule_activity.R
#          [--region PFC] [--fig-tfs STAT1,STAT3,CEBPD,FOSL2,BACH1,REL] [--out <file name>]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages(library(dplyr))
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
R        <- getarg("--region", "PFC")
CT       <- "Astro"
FIG_TFS  <- strsplit(getarg("--fig-tfs", "STAT1,STAT3,CEBPD,FOSL2,BACH1,REL"), ",")[[1]]
OUT_NAME <- getarg("--out", "tf_selection_rule_activity.csv")
MIN_EXPR <- 0.05; MIN_UP <- 5; FDR <- 0.05
DE <- file.path(OUT_ROOT, "03-DE_analysis")
folder <- file.path(MAGICAL_ANALYZED, magical_tag(R, CT))

## up-regulated targets per TF, from the region's own astrocyte DE
ec  <- utils::read.csv(file.path(folder, "expanded_circuits.csv"), stringsAsFactors = FALSE)
reg <- as.data.frame(readRDS(file.path(DE, R, "results", sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq.RDS", R))))
reg <- reg[as.character(reg$cluster_id) == CT, ]
up_genes <- unique(reg$gene[reg$logFC > 0 & reg$p_adj.glb < 0.05])
t <- ec %>% group_by(TF = TFs) %>%
  summarise(up_targets = length(intersect(unique(Gene_symbol), up_genes)), .groups = "drop") %>%
  arrange(desc(up_targets)) %>% mutate(rank_up_targets = seq_len(dplyr::n()))

## the TF's own expression and cross-region RNA change
xr <- as.data.frame(readRDS(file.path(DE, "cross_region", "data", "all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")))
xr <- xr[as.character(xr$cluster_id) == CT, ]
own <- data.frame(TF = xr$gene, expr_frq_PWH = round(xr$avg_HIV_frq, 4),
                  TF_RNA_up_in_PWH = xr$logFC > 0 & xr$p_adj.glb < 0.05, stringsAsFactors = FALSE)

## motif activity, donor-aware across regions; one motif per TF name (the first listed)
cv <- utils::read.csv(file.path(OUT_ROOT, "04-ATAC_DE", "motif_chromvar", "crossregion",
                                sprintf("%s_chromvar_crossregion.csv", CT)), stringsAsFactors = FALSE)
cv <- cv[!duplicated(toupper(cv$name)), ]
mot <- data.frame(TF = toupper(cv$name), chromvar_motif = cv$name,
                  chromvar_diff_pooled_centered = round(cv$diff_pooled, 3),
                  chromvar_FDR_centered = cv$FDR, chromvar_FDR_plain = cv$FDR_plain,
                  chromvar_FDR_weighted = cv$FDR_weighted, stringsAsFactors = FALSE)

## a TF with no motif of its own takes a composite motif that names it (Smad2::Smad3 -> SMAD2)
comp <- do.call(rbind, lapply(seq_len(nrow(mot)), function(i) {
  parts <- toupper(strsplit(mot$chromvar_motif[i], "::", fixed = TRUE)[[1]])
  if (length(parts) < 2) return(NULL)
  cbind(data.frame(TF = parts, stringsAsFactors = FALSE), mot[i, -1])
}))
comp <- comp[!duplicated(comp$TF) & !comp$TF %in% mot$TF, ]   # first listed composite wins
mot  <- rbind(mot, comp)

res <- t %>% select(TF, rank_up_targets, up_targets) %>%
  left_join(own, by = "TF") %>% left_join(mot, by = "TF") %>%
  ## a TF absent from the expression table is not expressed
  mutate(expr_frq_PWH = ifelse(is.na(expr_frq_PWH), 0, expr_frq_PWH),
         TF_RNA_up_in_PWH = ifelse(is.na(TF_RNA_up_in_PWH), FALSE, TF_RNA_up_in_PWH)) %>%
  mutate(meets_rule_activity = !is.na(expr_frq_PWH) & expr_frq_PWH >= MIN_EXPR & up_targets >= MIN_UP &
           !is.na(chromvar_FDR_centered) & chromvar_diff_pooled_centered > 0 & chromvar_FDR_centered < FDR,
         in_rule_activity_figure = TF %in% FIG_TFS)
bad <- setdiff(FIG_TFS, res$TF[res$meets_rule_activity])
if (length(bad)) stop("figure TFs that do not pass the rule: ", paste(bad, collapse = ", "))
res$meets_rule_activity     <- ifelse(res$meets_rule_activity, "True", "False")
res$in_rule_activity_figure <- ifelse(res$in_rule_activity_figure, "True", "False")
utils::write.csv(res, file.path(folder, OUT_NAME), row.names = FALSE)
message(sprintf("[tf-rule] %s %s: %d TFs pass; drawn: %s", R, CT, sum(res$meets_rule_activity == "True"),
                paste(FIG_TFS, collapse = ", ")))
