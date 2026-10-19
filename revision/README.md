# JCI revision

Code for the revised manuscript. The stage folders (`01-`…`07-`) hold the published
scripts with the revision's changes applied; every change is marked `## REVISION:`.
Everything new for the revision is in this folder.

## What changed from the published analysis

1. **Ambient RNA removal.** All transcriptomic analyses were rerun on CellBender-processed
   libraries (CellBender 0.3, FPR 0.01, 150 epochs). CellBender defines the nuclei used for
   QC, clustering and annotation; differential expression is estimated from the
   uncorrected counts of those nuclei (see `tools/de_concordance_corrected.R` for the
   corrected-vs-uncorrected comparison).
2. **Donor-aware statistics.**
   - Regional DE: unchanged model (muscat pseudobulk, edgeR, `~ HIV + Age + Sex + Race + PMI`),
     one sample per donor per region.
   - Cross-region DE: limma-voom with donor as a blocking factor (`duplicateCorrelation`),
     because most donors contribute more than one region. edgeR is kept as a concordance check.
   - Astrocyte subtype markers: donor-paired pseudobulk (`~ donor + subtype`) instead of
     cell-level Wilcoxon tests; subtype abundance tested per donor.
   - Cell-cell communication: MultiNicheNet (donor-level) instead of NICHES.
3. **VST sample exclusion** by a stated criterion (share of cortical excitatory nuclei, threshold
   in `config/paths.R`) instead of a fixed list.
4. **Outlier-donor sensitivity analyses** (excluding the HIV encephalitis donor and the donor
   with additional pathologies; leave-one-donor-out).

## Setup

1. Container: `docker/Dockerfile` (Bioconductor image plus the R packages, MAGICAL and
   FIt-SNE used here). Build it as `scorch_seurat_v5:latest` or set `SCORCH_IMAGE`.
2. Data locations: copy `config/local_paths.R.example` to `config/local_paths.R` and fill in
   the paths (pipeline outputs, original analysis outputs, output folder). It is not tracked.
3. Subject metadata: `config/build_sample_metadata.R` rebuilds the sample sheet from the
   pipeline manifests into `OUT_ROOT/meta/` (not tracked; contains subject-level data).
4. Extra R packages for MultiNicheNet/GSVA: `tools/install_ccc_rlibs.R` (into `RLIBS_CCC`).

Run any script through the launcher, from the repository root:

    revision/docker/run.sh 03-DE_analysis/01-INS_edgeR.R
    revision/docker/run.sh revision/tools/cross_region_donor_aware.R

It mounts the locations in `local_paths.R` read-only, runs in the background and writes a
log to `logs/`. All outputs go under `OUT_ROOT`.

## Order

1. Stages `01-` to `07-` in order (as in the original repository).
2. After stage 01: `tools/patch_msn_subtypes.R`, `tools/combined_drop_1040.R`,
   `tools/astro_drop_1040_vst.R`, `tools/make_v3_layer_joined.R`, `tools/dump_combined_metadata.R`.
3. Cross-region donor-aware DE: `tools/cross_region_donor_aware.R`, then
   `tools/cross_region_donor_aware_table.R` (figures read it with `SCORCH_XR_DE=donor_aware`).
4. Any figure or table script below.

## Scripts by result

| Result | Scripts (`tools/` unless noted) |
|---|---|
| Ambient-RNA sensitivity (reviewer) | run_ambient_fraction.R, ambient_fraction.R, ambient_vs_lowexpr.R, de_concordance_corrected.R, ambient_concordance_figure.R, compare_combined_to_published.R |
| Cross-region DE, donor-aware | cross_region_donor_aware.R, cross_region_donor_aware_table.R, cross_region_negF_fallback.R, cross_region_paper_rule_compare.R |
| Fig 2 (ATAC, motifs) | atac_motif_chromvar.R, atac_motif_crossregion.R, fig2_motif_heatmap.R, fig2_tf_rna_vs_motif.R, fig2_atac_summary.R, fig2E_scatter_by_class.R |
| Fig 3 (DE pathways, heatmaps) | fig3_pathway_panels.R, fig3_updown_panels.R, fig3_crossregion_pathways.R, cross_region_celltype_heatmap.R |
| Fig 4 (region-pair scatters) | fig4_cross_region_scatter.R, fig4_common.R |
| Microglia | run_microglia_signatures.R, microglia_signature_score.R, microglia_signature_heatmap.R, microglia_signature_atac.R, micro_module_networks.R, micro_module_umap.R, micro_subgroup_umap.R; `external/ldam_scores.R` |
| Astrocytes (Fig 6 and supplement) | astro_reactive_contrast.R, fig1e_heatmap_rebuild.R, astro_gsea_hallmark_reactome.R, astro_subtype_de_donor.R, astro_composition_donor.R, astro_lodo_forest.R, astro_cluster_robustness.R, astro_subtype_condition_de.R, astro_manuscript_figures.R, astro_reviewer_figures.R, run_astro_reviewer_figures.R, run_astro_cluster_figures.R, astro_cluster_features.R, astro_marker_featureplot.R, astro_nuclei_counts.R, figS11_astro_motif_heatmap.R, figS14b_reactive_atac_pathways.R, figS14b_donor_gsea.R, figS14b_pooled_donor.R |
| Gene-regulatory networks (MAGICAL) | magical_frq_variant.R, magical_tf_selection.R, magical_network.R, magical_variant_compare.R |
| Cell-cell communication | `external/multiNicheNet.R` |
| Outlier-donor sensitivity | donor_lodo_microglia.R, donor_lodo_astro.R, sensitivity_negF_refit.R, donor_sensitivity_supp.R, ldam_ssgsea_sensitivity.R, multinichenet_sensitivity.R |
| VST sample exclusion | vst_ext_rank_figure.R, vst_ext_marker_table.R, vst_addback_de.R, run_vst_sensitivity.sh |
| Viral reads per nucleus | viral_threshold_table.R, viral_nuclei_table.R, viral_figure.R |
| Other supplementary figures | figS_violin_panels.R, immune_subset_features.R |
| Supplementary tables | table1_cohort.R + table1_to_xlsx.py, supp_deg_tables.R, supp_qc_table.py, legend_donor_counts.R, donor_representation.R |

`external/` holds analyses contributed by Junchen (MultiNicheNet, LDAM ssGSEA).
