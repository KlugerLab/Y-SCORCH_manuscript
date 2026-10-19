#!/usr/bin/env bash
# ==============================================================================
# VST sample-exclusion sensitivity: regional VST edgeR DE under two sample sets,
# to show how much the VST findings depend on the exclusions (reviewer request).
#
# Method
#   Runs 03-DE_analysis/10-sensitivity_DE.R twice with VST_SAMPLE_EXCLUDE set:
#     published26  the paper's 7-sample removal               -> 16 HIV / 10 CTR
#     criterion29  the stated exclusion criterion applied to  -> 18 HIV / 11 CTR
#                  the corrected data; adds back 1119, HCTME, 7200776574
#   Requires the VST chain to be complete (03 -> 04-01 -> combined -> VST/05).
#   The two results are compared with 03-DE_analysis/11-compare_sensitivity.R
#   (command printed at the end).
#
# Outputs  OUT_ROOT/03-DE_analysis/VST/results/
#            VST_published26_res_edgeR_filtered_withFreq.RDS
#            VST_criterion29_res_edgeR_filtered_withFreq.RDS
#
# Usage
#   revision/tools/run_vst_sensitivity.sh
# ==============================================================================
set -euo pipefail
cd "$(dirname "$0")/../.."

PUBLISHED="1039,1068,1119,7200776574,HCCRS,HCTME,HCTXJ"
CRITERION="1039,1068,HCCRS,HCTXJ"

echo "== published cohort (26 samples) =="
VST_SAMPLE_EXCLUDE="$PUBLISHED" revision/docker/run.sh \
  03-DE_analysis/10-sensitivity_DE.R --region VST --label published26

echo "== criterion cohort (29 samples) =="
VST_SAMPLE_EXCLUDE="$CRITERION" revision/docker/run.sh \
  03-DE_analysis/10-sensitivity_DE.R --region VST --label criterion29

cat <<'MSG'

Both launched in the background. When they finish, compare with:

  revision/docker/run.sh 03-DE_analysis/11-compare_sensitivity.R \
    --a VST/results/VST_published26_res_edgeR_filtered_withFreq.RDS --a-label published26 \
    --b VST/results/VST_criterion29_res_edgeR_filtered_withFreq.RDS --b-label criterion29 \
    --out vst_exclusions
MSG
