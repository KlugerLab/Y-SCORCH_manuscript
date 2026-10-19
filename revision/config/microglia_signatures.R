# ==============================================================================
# Microglial state signatures for donor-level enrichment testing.
#
# Curated gene sets from the literature, with provenance and caveats attached to
# each set; the reporting code reads them and prints them alongside every result.
# Sourced by the microglia signature tools (e.g. revision/tools/microglia_signature_score.R).
#
# Two caveats that are reported with every result:
#
# 1. snRNA-seq depletion (Thrupp 2020). Nuclei are depleted for ~1% of genes
#    that are disproportionately microglial activation genes - APOE, CST3, SPP1,
#    CD74 among them, ~18% of previously defined disease-associated genes. The
#    DAM/LDAM/MGnD sets below are therefore expected to underperform in nuclei.
#    `THRUPP_DEPLETED` lists the named ones so every report can say which of a
#    set's genes are affected. This is an argument for set-level statistics over
#    single-gene calls, and for pseudobulk over per-nucleus scoring: summing
#    thousands of nuclei per donor recovers genes that are sparse in any one.
#
# 2. Ambient RNA, measured in this dataset. In INS microglia:
#       IFITM3 0.33   LGALS3 0.30   ISG15 0.27   CST3 0.23   FTH1 0.18
#       CD63   0.17   APOE   0.02   SPP1  0.01   TREM2/P2RY12/CD74 ~0.00
#    IFITM3 and ISG15 are two of the highest-expressed members of IRM_core, so
#    the interferon signature is the most ambient-exposed set here. APOE and
#    SPP1 come in near zero, so the DAM sets are cleaner on this axis than
#    expected - their problem is Thrupp depletion, not soup.
#    `AMBIENT_SUSPECT` drives the drop-ambient sensitivity analysis.
#
# SPECIES. Sets marked species_origin "mouse" were derived in mouse but the
# symbols here are already human (HGNC). This file does not re-map orthologs;
# the reporting code instead reports the detection rate of each set in the data,
# which is the number that matters. Symbols that are not 1:1 orthologs would
# show up as undetected - check the detection report before concluding a set is
# genuinely absent.
#
# SET SIZE. Four sets have <= 4 genes (LDAM_Haney2024, MG4_lipid_processing,
# HIV_microglia_innate, IRM_high_threshold). They are computed but flagged
# `underpowered`, and should not be headline results.
#
# OVERLAP. DAM_stage2, DAM_Haney2024 and MGnD_Krasemann2017 share most of their
# genes. Three concordant scores there are one result, not three. The reporting
# code prints the Jaccard matrix so this cannot be overlooked.
# ==============================================================================

## Named in Thrupp 2020 as depleted in nuclei relative to cells.
THRUPP_DEPLETED <- c("APOE", "CST3", "SPP1", "CD74")

## Ambient fraction >= 0.15 in INS microglia (revision/tools/ambient_fraction.R output).
## The sensitivity analysis re-scores every set with these removed.
AMBIENT_SUSPECT <- c("IFITM3", "LGALS3", "ISG15", "CST3", "FTH1", "CD63")

MICROGLIA_SIGNATURES <- list(

  # ---- LDAM (lipid-droplet-accumulating microglia) --------------------------
  LDAM_Haney2024 = list(
    genes = c("ACSL1", "NAMPT", "DPYD", "CD163"),
    direction = "up", species_origin = "human", source = "haney2024",
    note = paste("Marker panel used to annotate the LDAM subcluster in human AD",
                 "snRNA-seq. ACSL1 is the defining gene. Only 4 genes -",
                 "underpowered as a set statistic; treat as supporting.")),

  LDAM_Marschallinger2020_up = list(
    genes = c("PLIN3", "ACLY",
              "CD63", "ATP6V1A", "ATP6V1C1", "ATP6V1G1", "TUBA1A",
              "RAB5B", "RAB7A", "CD22",
              "CAT", "KL", "PPP1CB", "RAP1B"),
    direction = "up", species_origin = "mouse", source = "marschallinger2020",
    note = paste("Main-text genes only. The full state is 692 DEGs",
                 "(Supplementary Table 1) - pull that list for a properly",
                 "powered test.")),

  LDAM_Marschallinger2020_down = list(
    genes = c("ADRB1", "ADRB2", "AXL", "CD74", "CLEC7A", "CYBB"),
    direction = "down", species_origin = "mouse", source = "marschallinger2020",
    note = paste("Directionality is the discriminator vs DAM: AXL/CD74/CLEC7A/",
                 "CYBB are UP in DAM but DOWN in LDAM. Marschallinger report",
                 "APOE and TREM2 as NOT regulated in LDAM. A population scoring",
                 "high on APOE/TREM2 is more likely DAM-like.")),

  MG4_lipid_processing_Sun2023 = list(
    genes = c("PPARG", "APOE", "ABCA1", "TREM2"),
    direction = "up", species_origin = "human", source = "sun2023",
    note = paste("Main-text exemplars only. MG4 is APOE/TREM2-HIGH, so it is",
                 "NOT the same state as Marschallinger LDAM - test them as",
                 "separate hypotheses, never as synonyms.")),

  # ---- Homeostatic / DAM ---------------------------------------------------
  Homeostatic = list(
    genes = c("P2RY12", "P2RY13", "CX3CR1", "TMEM119", "TGFBR1", "HEXB",
              "CSF1R", "GPR34", "SELPLG"),
    direction = "up", species_origin = "mouse/human",
    source = "kerenshaul2017,haney2024"),

  DAM_stage1_KerenShaul2017 = list(
    genes = c("TYROBP", "APOE", "B2M", "FTH1", "CTSB", "CTSD"),
    direction = "up", species_origin = "mouse", source = "kerenshaul2017",
    note = "TREM2-independent; accompanied by loss of homeostatic genes."),

  DAM_stage2_KerenShaul2017 = list(
    genes = c("TREM2", "CST7", "LPL", "CD9", "ITGAX", "AXL", "CLEC7A",
              "LILRB4", "CTSL", "CSF1", "CD68", "SPP1"),
    direction = "up", species_origin = "mouse", source = "kerenshaul2017",
    note = "TREM2-dependent; lysosomal, phagocytic, lipid metabolism."),

  DAM_Haney2024 = list(
    genes = c("CD9", "ITGAX", "CLEC7A", "CD63", "SPP1", "LPL", "TREM2", "APOE"),
    direction = "up", species_origin = "human", source = "haney2024",
    note = paste("Human-symbol DAM panel from the same annotation scheme as",
                 "LDAM_Haney2024 - the cleanest like-for-like comparator.")),

  MGnD_Krasemann2017 = list(
    genes = c("APOE", "CLEC7A", "ITGAX", "AXL", "SPP1", "LGALS3", "CCL2",
              "GPNMB", "CD63", "TREM2"),
    direction = "up", species_origin = "mouse", source = "krasemann2017",
    note = "APOE-TREM2 driven; heavily overlapping with DAM stage 2."),

  # ---- IRM (interferon-responsive microglia) --------------------------------
  IRM_core = list(
    genes = c("IFITM3", "IFIT1", "IFIT2", "IFIT3", "IRF7", "STAT1", "ISG15",
              "MX1", "OASL", "OAS1", "RSAD2", "BST2", "USP18", "IFI44L", "IFI6"),
    direction = "up", species_origin = "mouse/human",
    source = "salafrigerio2019,olah2020,lopezatalaya2025",
    note = paste("Essentially the type-I IFN program. ALSO test",
                 "HALLMARK_INTERFERON_ALPHA/GAMMA_RESPONSE - better powered and",
                 "more credible than a 15-gene hand-curated list. NOTE: IFITM3",
                 "(0.33) and ISG15 (0.27) are the most ambient-exposed genes in",
                 "this dataset's microglia; the drop-ambient sensitivity run is",
                 "load-bearing for this set specifically.")),

  IRM_high_threshold = list(
    genes = c("CXCL10", "CD38"),
    direction = "up", species_origin = "human", source = "irm_crispri2026",
    note = "2 genes. Separates weak from strong IFN response; not a set test."),

  # ---- HIV-associated ------------------------------------------------------
  # There is no consensus named HIV-microglia signature. This is a placeholder
  # set (exploratory = TRUE) and is not reported on its own four genes; it is
  # to be populated from the Plaza-Jennings 2022 and Trim 2025 supplementary DEG
  # tables before reporting anything from it.
  HIV_microglia_innate_TRIM2025 = list(
    genes = c("NLRP3", "TLR2", "TLR4", "CD163"),
    direction = "up", species_origin = "human", source = "trim2025",
    exploratory = TRUE,
    note = paste("PLACEHOLDER - main-text exemplars only. Replace with the full",
                 "supplementary DEG list before this appears in any figure."))
)

MICROGLIA_CITATIONS <- c(
  marschallinger2020 = "Marschallinger J, et al. Nat Neurosci. 2020;23(2):194-208. doi:10.1038/s41593-019-0566-1",
  haney2024          = "Haney MS, et al. Nature. 2024;628(8006):154-161. doi:10.1038/s41586-024-07185-7",
  kerenshaul2017     = "Keren-Shaul H, et al. Cell. 2017;169(7):1276-1290.e17. doi:10.1016/j.cell.2017.05.018",
  krasemann2017      = "Krasemann S, et al. Immunity. 2017;47(3):566-581.e9. doi:10.1016/j.immuni.2017.08.008",
  salafrigerio2019   = "Sala Frigerio C, et al. Cell Rep. 2019;27(4):1293-1306.e6. doi:10.1016/j.celrep.2019.03.099",
  olah2020           = "Olah M, et al. Nat Commun. 2020;11:6129. doi:10.1038/s41467-020-19737-2",
  sun2023            = "Sun N, et al. Cell. 2023;186(20):4386-4403.e29. doi:10.1016/j.cell.2023.08.037",
  plazajennings2022  = "Plaza-Jennings AL, et al. Mol Cell. 2022;82(24):4647-4663.e8. doi:10.1016/j.molcel.2022.11.016",
  farhadian2018      = "Farhadian SF, et al. JCI Insight. 2018;3(18):e121718. doi:10.1172/jci.insight.121718",
  trim2025           = "PMID 42357873. [VERIFY full citation before submission]",
  thrupp2020         = "Thrupp N, et al. Cell Rep. 2020;32(13):108189. doi:10.1016/j.celrep.2020.108189",
  lopezatalaya2025   = "Lopez-Atalaya JP, et al. FEBS J. 2025. doi:10.1111/febs.70126",
  irm_crispri2026    = "npj Dementia. 2026. doi:10.1038/s44400-026-00125-4"
)

#' Genes and citations for one signature.
microglia_signature <- function(name) {
  s <- MICROGLIA_SIGNATURES[[name]]
  if (is.null(s)) stop("unknown signature '", name, "'")
  list(genes = s$genes,
       citations = unname(MICROGLIA_CITATIONS[strsplit(s$source, ",")[[1]]]))
}

#' Jaccard overlap between every pair of signatures.
#'
#' Print this before interpreting. DAM_stage2 / DAM_Haney2024 /
#' MGnD_Krasemann2017 overlap heavily and their scores will move together;
#' reporting them as three independent confirmations would be wrong.
microglia_signature_overlap <- function(sigs = MICROGLIA_SIGNATURES) {
  n <- names(sigs)
  m <- outer(n, n, Vectorize(function(a, b) {
    ga <- sigs[[a]]$genes; gb <- sigs[[b]]$genes
    length(intersect(ga, gb)) / length(union(ga, gb))
  }))
  dimnames(m) <- list(n, n)
  round(m, 3)
}

#' One row per signature: size, provenance, and the two caveat counts.
microglia_signature_summary <- function(sigs = MICROGLIA_SIGNATURES) {
  do.call(rbind, lapply(names(sigs), function(n) {
    s <- sigs[[n]]
    data.frame(
      signature      = n,
      n_genes        = length(s$genes),
      direction      = s$direction,
      species_origin = s$species_origin,
      n_thrupp       = sum(s$genes %in% THRUPP_DEPLETED),
      n_ambient      = sum(s$genes %in% AMBIENT_SUSPECT),
      underpowered   = length(s$genes) <= 4,
      exploratory    = isTRUE(s$exploratory),
      stringsAsFactors = FALSE)
  }))
}
