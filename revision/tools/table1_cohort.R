# ==============================================================================
# Table 1: summarised description of the donor cohort, PWoH vs PWH (reviewer
# request for a summary table in addition to the full per-donor data file).
#
# Replaces the grouped table in the published supplement
# (patient_level_information_and_statistics.xlsx) with one covering what the
# analyses depend on; the full per-donor sheet (detailed_sample_metadata.xlsx)
# remains the data file.
#
# Method
#   Cohort: rows of SAMPLE_META with a cellranger_path and oud_jay_v2 == FALSE,
#   33 donors (20 PWH / 13 PWoH). Samples analysed: snRNA = every regional
#   library minus the VST exclusions; snATAC = libraries in
#   atac_nuclei_by_sample.csv.
#   Statistics: continuous variables as median (IQR) [range] with Mann-Whitney U
#   (wilcox.test, normal approximation with ties); categorical as n (%) with
#   Fisher's exact test (r x 2). Tests only for variables recorded in both
#   groups. Brain bank is completely confounded with condition, so it has no
#   p-value.
#   Corrections to the published grouped table, checked against the source:
#   "% HIV RNA < 50 copies/ml = 0.33" was the < 200 count (6/18); < 50 is 4/18.
#   HIV-duration bins were labelled < 5 / 6-14 / > 15; the source bins are
#   0-5 / 6-10 / 11-15 / >= 16. CD4 units are cells/uL, not cells/ml. CD4 and
#   nadir are available for 19/20 PWH, not 18/20.
#
# Inputs   OUT_ROOT/meta/sample_metadata_v3.csv (SAMPLE_META; subject-level, not
#          distributed)
#          OUT_ROOT/01-Regional_data_processing/VST/data/04-VST-tsne_marker_plotting/
#            sample_exclusion_decision.csv
#          OUT_ROOT/02-ATAC_processing/atac_nuclei_by_sample.csv
# Outputs  OUT_ROOT/tables/table1_cohort/table1_cohort.{csv,tex,md}
#          (revision/tools/table1_to_xlsx.py then writes the .xlsx)
#
# Usage
#   revision/docker/run.sh revision/tools/table1_cohort.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
OUT <- file.path(OUT_ROOT, "tables", "table1_cohort")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

m <- utils::read.csv(SAMPLE_META, stringsAsFactors = FALSE, colClasses = "character")
m <- m[!is.na(m$cellranger_path) & nzchar(m$cellranger_path) & tolower(m$oud_jay_v2) == "false", ]
d <- m[!duplicated(m$sample_ids), ]
d$grp <- ifelse(toupper(d$hiv_jay_v2) == "TRUE", "PWH", "PWoH")
G <- c("PWoH", "PWH")
N <- table(factor(d$grp, levels = G))
stopifnot(N[["PWH"]] == 20, N[["PWoH"]] == 13)
na_txt <- function(x) is.na(x) | !nzchar(x) | x %in% c("NA", "not_applicable")

rows <- list()
add <- function(section, variable, vals, p = NA, note = "")
  rows[[length(rows) + 1]] <<- data.frame(section = section, variable = variable,
    PWoH = vals[["PWoH"]], PWH = vals[["PWH"]], p = p, note = note, stringsAsFactors = FALSE)
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.2f", p)))
num_summ <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (!length(x)) return("—")
  q <- stats::quantile(x, c(.25, .5, .75))
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f) [%.", digits, "f–%.", digits, "f]"),
          q[2], q[1], q[3], min(x), max(x))
}
cnt <- function(k, n) if (n == 0) "—" else sprintf("%d (%.0f%%)", k, 100 * k / n)
by_grp <- function(f) setNames(lapply(G, function(g) f(d[d$grp == g, , drop = FALSE])), G)
num_row <- function(section, variable, col, digits = 1, test = TRUE, note = "") {
  x <- suppressWarnings(as.numeric(d[[col]]))
  v <- by_grp(function(s) num_summ(suppressWarnings(as.numeric(s[[col]])), digits))
  avail <- tapply(is.finite(x), d$grp, sum)[G]
  part <- avail > 0 & avail < N[G]                      # a group with nothing recorded shows an em-dash instead
  if (any(part)) note <- trimws(paste(note, sprintf("available for %s", paste(sprintf("%d/%d %s", avail, N[G], G)[part], collapse = ", "))))
  p <- if (test) suppressWarnings(stats::wilcox.test(x[d$grp == "PWH"], x[d$grp == "PWoH"])$p.value) else NA
  add(section, variable, v, p, note)
}
## One header row plus one row per level. A group with nothing recorded shows an
## em-dash throughout, so PWH-only variables do not read "0 (0%)" for PWoH.
cat_rows <- function(section, variable, lev, levs = NULL, test = TRUE, note = "") {
  f <- factor(lev, levels = if (is.null(levs)) sort(unique(stats::na.omit(lev))) else levs)
  g <- factor(d$grp, levels = G)
  tb <- table(f, g)
  rec <- tapply(!is.na(f), g, sum)
  p <- if (test && nrow(tb) > 1 && all(rec > 0)) stats::fisher.test(tb)$p.value else NA
  add(section, variable, c(PWoH = "", PWH = ""), p, note)
  for (l in levels(f))
    add(section, paste0("   ", l), sapply(G, function(gg) if (rec[[gg]] == 0) "\u2014" else cnt(tb[l, gg], N[[gg]])))
  nr <- N[G] - rec
  if (any(nr[rec > 0] > 0))
    add(section, "   Not recorded", sapply(G, function(gg) if (rec[[gg]] == 0) "\u2014" else cnt(nr[[gg]], N[[gg]])))
}
## a single yes/no characteristic on one row, Fisher's exact
bin_row <- function(section, variable, flag) {
  v <- sapply(G, function(gg) cnt(sum(flag[d$grp == gg], na.rm = TRUE), N[[gg]]))
  tb <- table(factor(flag, levels = c(TRUE, FALSE)), factor(d$grp, levels = G))
  add(section, variable, v, stats::fisher.test(tb)$p.value)
}
only_pwh <- function(v) ifelse(d$grp == "PWH", v, NA)

## ---- donors and samples -------------------------------------------------------------
add("Cohort", "Donors, n", c(PWoH = as.character(N[["PWoH"]]), PWH = as.character(N[["PWH"]])))
excl <- utils::read.csv(file.path(OUT_ROOT, "01-Regional_data_processing", "VST", "data", "04-VST-tsne_marker_plotting",
                                  "sample_exclusion_decision.csv"), stringsAsFactors = FALSE, colClasses = "character")
vst_out <- excl$sample_ids[toupper(excl$excluded) == "TRUE"]
rna <- m[, c("sample_ids", "region")]; rna <- rna[!(rna$region == "VST" & rna$sample_ids %in% vst_out), ]
atac <- utils::read.csv(file.path(OUT_ROOT, "02-ATAC_processing", "atac_nuclei_by_sample.csv"), stringsAsFactors = FALSE,
                        colClasses = c("character", "character", "numeric", "character"))
names(atac)[names(atac) == "sample"] <- "sample_ids"
for (mod in c("snRNA-seq", "snATAC-seq")) {
  s <- if (mod == "snRNA-seq") rna else atac
  add("Samples analysed", paste0(mod, " libraries passing QC"), c(PWoH = "", PWH = ""))
  for (R in c("INS", "PFC", "VST")) {
    ids <- unique(s$sample_ids[s$region == R])
    add("Samples analysed", paste0("   ", R), c(PWoH = cnt(sum(d$sample_ids[d$grp == "PWoH"] %in% ids), N[["PWoH"]]),
                                             PWH = cnt(sum(d$sample_ids[d$grp == "PWH"] %in% ids), N[["PWH"]])))
  }
  k <- table(s$sample_ids)
  all3 <- names(k)[k == 3]
  add("Samples analysed", "   All three regions", c(PWoH = cnt(sum(d$sample_ids[d$grp == "PWoH"] %in% all3), N[["PWoH"]]),
                                                    PWH = cnt(sum(d$sample_ids[d$grp == "PWH"] %in% all3), N[["PWH"]])))
}

## ---- demographics and post-mortem ---------------------------------------------------
num_row("Demographics", "Age at death, years", "age_value_v2")
bin_row("Demographics", "Female sex", d$sex_v2 == "female")
race <- ifelse(grepl("Native", d$race_v2), "American Indian or Alaska Native (alone or with White)",
         ifelse(d$race_v2 == "White", "White", ifelse(grepl("Black", d$race_v2), "Black or African American", NA)))
cat_rows("Demographics", "Race", race, levs = c("Black or African American", "White",
                                                 "American Indian or Alaska Native (alone or with White)"))
bin_row("Demographics", "Hispanic or Latino ethnicity", d$ethnicity_v2 == "Hispanic_or_Latino")
num_row("Post-mortem", "Post-mortem interval, hours", "postmortem_interval_v2")
bank <- c(NSU = "Nova Southeastern University", `NNTC:MHBB` = "NNTC: Manhattan HIV Brain Bank",
          `NNTC:NNAB` = "NNTC: National Neurological AIDS Bank", `NNTC:TNRC` = "NNTC: Texas NeuroAIDS Research Center")
cat_rows("Post-mortem", "Brain bank", unname(bank[d$subject_source_v2]), levs = unname(bank), test = FALSE,
         note = "Condition and brain bank are completely confounded; no test")
cod <- ifelse(na_txt(d$cause_of_death_v2), NA,
        ifelse(grepl("GSW|trauma", d$cause_of_death_v2, ignore.case = TRUE), "Traumatic (gunshot wound, blunt force)",
        ifelse(grepl("cardi|myocard|heart|atheros", d$cause_of_death_v2, ignore.case = TRUE), "Cardiac", "Other")))
cat_rows("Post-mortem", "Cause of death", cod, levs = c("Cardiac", "Traumatic (gunshot wound, blunt force)"), test = FALSE,
         note = "Not recorded for PWH in the source data; one PWoH cardiac death was viral myocarditis")

## ---- HIV clinical (PWH only) --------------------------------------------------------
dur <- c(`0-5_years` = "0–5", `6-10_years` = "6–10", `11-15_years` = "11–15", `16_or_more_years` = "≥16")
cat_rows("HIV clinical (PWH)", "Estimated duration of HIV infection, years", only_pwh(unname(dur[d$years_with_hiv_v2])),
         levs = unname(dur), test = FALSE)
art <- ifelse(na_txt(d$antiretroviral_treatment_v2), NA,
        ifelse(d$antiretroviral_treatment_v2 == "no_art", "None recorded", "Regimen recorded"))
cat_rows("HIV clinical (PWH)", "Antiretroviral therapy at last assessment", only_pwh(art),
         levs = c("Regimen recorded", "None recorded"), test = FALSE)
vl_raw <- d$plasma_viral_load_measurement_v2
vl_num <- suppressWarnings(as.numeric(vl_raw))
vl <- ifelse(grepl("below_det_limit_50", vl_raw) | (!is.na(vl_num) & vl_num < 50), "<50",
       ifelse(grepl("below_det_limit_200", vl_raw) | (!is.na(vl_num) & vl_num < 200), "50–199 (or <200, assay limit)",
       ifelse(grepl("above_det_limit", vl_raw) | (!is.na(vl_num) & vl_num >= 10000), "≥10,000",
       ifelse(!is.na(vl_num), "200–9,999", NA))))
vl <- factor(vl, levels = c("<50", "50–199 (or <200, assay limit)", "200–9,999", "≥10,000"))
cat_rows("HIV clinical (PWH)", "Last plasma HIV RNA, copies/mL", only_pwh(as.character(vl)), levs = levels(vl), test = FALSE)
num_row("HIV clinical (PWH)", "Months from last viral load to death", "months_since_last_plasma_viral_measurement_v2", test = FALSE)
num_row("HIV clinical (PWH)", "Last CD4+ T-cell count, cells/µL", "cd4_counts_v2", digits = 0, test = FALSE)
cd4 <- suppressWarnings(as.numeric(d$cd4_counts_v2))
add("HIV clinical (PWH)", "   Last CD4+ <200 cells/µL", c(PWoH = "—", PWH = cnt(sum(cd4[d$grp == "PWH"] < 200, na.rm = TRUE), sum(is.finite(cd4[d$grp == "PWH"])))),
    note = sprintf("%% of the %d PWH with a count", sum(is.finite(cd4[d$grp == "PWH"]))))
num_row("HIV clinical (PWH)", "Nadir CD4+ T-cell count, cells/µL", "cd4_nadir_v2", digits = 0, test = FALSE)

cog <- d$cognitive_status_v2
cogc <- ifelse(na_txt(cog) | cog == "Test not done", NA,
         ifelse(grepl("HAD|Mild Neurocognitive", cog), "HIV-associated neurocognitive disorder (HAD or MND)",
         ifelse(grepl("Normal", cog), "Neurocognitively normal",
         ifelse(grepl("other cause", cog, ignore.case = TRUE), "Impairment attributed to other causes",
         ifelse(grepl("Unable", cog), "Unable to assign a diagnosis", "Other")))))
cat_rows("HIV clinical (PWH)", "Neurocognitive diagnosis (NNTC)", only_pwh(cogc),
         levs = c("HIV-associated neurocognitive disorder (HAD or MND)", "Neurocognitively normal",
                  "Impairment attributed to other causes", "Unable to assign a diagnosis"), test = FALSE,
         note = "Not recorded includes 'test not done'")

## neuropathology: findings are not mutually exclusive, so one row per finding
bp <- d$brain_pathology_v2
pw <- d$grp == "PWH"
np <- list(
  "HIV encephalitis" = grepl("HIV Encephalitis", bp),
  "Other CNS infection (incl. PML, bacterial)" = grepl("Progressive multifocal|Other infections|Bacterial", bp),
  "Cerebrovascular or hypoxic-ischaemic change" = grepl("infarct|Atherosclerosis|Hemorrhage|Hypoxic|Arteriolar|Thrombus|siderocalcinosis", bp),
  "Alzheimer type II gliosis" = grepl("Alzheimer type 2", bp),
  "Neurofibrillary or senile-plaque pathology" = grepl("Neurofibrillary|Senile plaque", bp),
  "Normal or minimal non-diagnostic abnormalities only" = bp %in% c("Brain Normal", "Minimal non-diagnostic abnormalities"))
add("HIV clinical (PWH)", "Neuropathology (NNTC; findings overlap)", c(PWoH = "", PWH = ""),
    note = "A donor can have several findings; not recorded for PWoH")
for (k in names(np)) add("HIV clinical (PWH)", paste0("   ", k), c(PWoH = "—", PWH = cnt(sum(np[[k]] & pw), N[["PWH"]])))
nrec <- sum(na_txt(bp) & pw)
if (nrec) add("HIV clinical (PWH)", "   Not recorded", c(PWoH = "—", PWH = cnt(nrec, N[["PWH"]])))

## ---- write ---------------------------------------------------------------------------
T1 <- do.call(rbind, rows)
T1$p <- fmt_p(T1$p)
names(T1)[3:4] <- c(sprintf("PWoH (n = %d)", N[["PWoH"]]), sprintf("PWH (n = %d)", N[["PWH"]]))
utils::write.csv(T1, file.path(OUT, "table1_cohort.csv"), row.names = FALSE, fileEncoding = "UTF-8")

foot <- c("Continuous variables: median (interquartile range) [range]; categorical: n (%). P values: Mann-Whitney U",
          "(continuous) or Fisher's exact test (categorical), for variables recorded in both groups. Age, sex, race and",
          "post-mortem interval are covariates in every differential model. Donors with opioid",
          "use disorder were excluded by design. Post-mortem toxicology was recorded for 3 PWH only and is given in the",
          "full data file. NNTC, National NeuroAIDS Tissue Consortium; HAD, HIV-associated dementia; MND, mild",
          "neurocognitive disorder; PML, progressive multifocal leukoencephalopathy.")
md <- c(sprintf("| %s | %s | %s | P | Note |", "Characteristic", names(T1)[3], names(T1)[4]), "|---|---|---|---|---|")
sec <- ""
for (i in seq_len(nrow(T1))) {
  if (T1$section[i] != sec) { sec <- T1$section[i]; md <- c(md, sprintf("| **%s** | | | | |", sec)) }
  md <- c(md, sprintf("| %s | %s | %s | %s | %s |", sub("^   ", "&nbsp;&nbsp;&nbsp;", T1$variable[i]), T1[i, 3], T1[i, 4], T1$p[i], T1$note[i]))
}
writeLines(c(md, "", paste(foot, collapse = " ")), file.path(OUT, "table1_cohort.md"), useBytes = TRUE)

tex_esc <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x); x <- gsub("([%&#_])", "\\\\\\1", x)
  x <- gsub("≥", "$\\\\geq$", x); x <- gsub("–", "--", x); x <- gsub("—", "---", x)
  x <- gsub("µ", "$\\\\mu$", x); x <- gsub("<", "$<$", x); x
}
tex <- c("\\begin{table}[!t]", "\\centering", "\\footnotesize",
         "\\caption{Characteristics of the donors.}", "\\label{tab:cohort}",
         "\\begin{tabular}{p{6.4cm}lll}", "\\toprule",
         sprintf("Characteristic & %s & %s & \\textit{P} \\\\", tex_esc(names(T1)[3]), tex_esc(names(T1)[4])), "\\midrule")
## the Note column becomes lettered footnotes (the .csv/.xlsx keep it as a column)
sec <- ""; notes <- character(0)
for (i in seq_len(nrow(T1))) {
  if (T1$section[i] != sec) { sec <- T1$section[i]; tex <- c(tex, sprintf("\\multicolumn{4}{l}{\\textbf{%s}} \\\\", tex_esc(sec))) }
  v <- T1$variable[i]; ind <- grepl("^   ", v); v <- tex_esc(sub("^   ", "", v))
  if (nzchar(T1$note[i])) { notes <- c(notes, T1$note[i]); v <- sprintf("%s\\textsuperscript{%s}", v, letters[length(notes)]) }
  tex <- c(tex, sprintf("%s%s & %s & %s & %s \\\\", if (ind) "\\quad " else "", v, tex_esc(T1[i, 3]), tex_esc(T1[i, 4]), tex_esc(T1$p[i])))
}
note_txt <- paste(sprintf("\\textsuperscript{%s}%s.", letters[seq_along(notes)], tex_esc(notes)), collapse = " ")
tex <- c(tex, "\\bottomrule", "\\end{tabular}",
         sprintf("\\par\\smallskip\\parbox{\\linewidth}{\\scriptsize %s %s}", tex_esc(paste(foot, collapse = " ")), note_txt),
         "\\end{table}", "% needs \\usepackage{booktabs} in the preamble")
writeLines(tex, file.path(OUT, "table1_cohort.tex"), useBytes = TRUE)
print(T1[, c(2:5)], row.names = FALSE, right = FALSE)
message("[t1] wrote ", OUT)
