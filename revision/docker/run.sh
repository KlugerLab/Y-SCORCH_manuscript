#!/usr/bin/env bash
# ==============================================================================
# Launch one analysis script inside the container, in the background, with a log.
#
#   revision/docker/run.sh 01-Regional_data_processing/INS/01-INS-mt_removal.R [Rscript args...]
#
# Jobs are started with nohup. The log lands in
# logs/<script-basename>.<timestamp>.<pid>.log (not tracked). The exit code alone
# is not enough: read the log and check the expected outputs exist before
# starting the next stage.
# ==============================================================================
set -euo pipefail

IMAGE="${SCORCH_IMAGE:-scorch_seurat_v5:latest}"
## repository root = two levels above this file; data locations from local_paths.R
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOCAL="$REPO/revision/config/local_paths.R"
[[ -f "$LOCAL" ]] || { echo "missing $LOCAL - copy local_paths.R.example" >&2; exit 1; }
lp() { sed -nE "s/^$1[[:space:]]*<-[[:space:]]*\"([^\"]+)\".*/\1/p" "$LOCAL"; }
OUT_ROOT="$(lp OUT_ROOT)"

SCRIPT="${1:?usage: run.sh <script.R> [Rscript args...]}"
shift || true

[[ -f "$REPO/$SCRIPT" || -f "$SCRIPT" ]] || { echo "no such script: $SCRIPT" >&2; exit 1; }

mkdir -p "$REPO/logs" "$OUT_ROOT"
# The $$ suffix keeps two launches of the same script in the same second (e.g.
# one per region) from sharing and interleaving a log.
LOG="$REPO/logs/$(basename "${SCRIPT%.*}").$(date +%Y%m%d-%H%M%S).$$.log"

# R's tempdir() is kept off the container's writable layer, which sits on the
#   host root filesystem. Signac writes per-group BED files to tempdir() before
#   handing them to MACS2; if that filesystem fills, the files are silently
#   truncated mid-line and MACS2 dies with
#     StrandFormatError: ... "chr2 230127066 230127183 AGGAACCAGchr1 10103 ..."
#   (two BED records run together), which looks like corrupt data but is a
#   failed write. TMPDIR points under OUT_ROOT instead.
SCRATCH="$OUT_ROOT/.tmp"
mkdir -p "$SCRATCH"

# Read-only mounts for every data source; only OUT_ROOT and the repo are writable.
#
# HDF5_USE_FILE_LOCKING=FALSE is required for reading .h5 files from an NFS
#   mount without a working lock manager (NFSv3, `hard`, `local_lock=none`).
#   HDF5 >= 1.10 takes a POSIX lock on every open, including read-only, so the
#   first Read10X_h5() blocks and, because the mount is `hard`, retries forever
#   at 0% CPU. Plain reads (cp, readRDS) never request a lock, which is why
#   copying files works while opening HDF5 does not.
#   Disabling the lock is safe here: the data mounts are read-only, and the lock
#   exists solely to serialise writers.
# Forward the environment variables the scripts read (SCORCH_DE_COUNTS =
#   uncorrected|corrected, SCORCH_RACE_CODING, VST_SAMPLE_EXCLUDE,
#   SCORCH_OUT_SUFFIX, SCORCH_XR_DE). Without this,
#   `SCORCH_DE_COUNTS=corrected revision/docker/run.sh ...` would silently run
#   the primary analysis.
ENVFWD=()
for v in SCORCH_DE_COUNTS SCORCH_RACE_CODING SCORCH_OUT_SUFFIX VST_SAMPLE_EXCLUDE SCORCH_XR_DE SCORCH_EDGER_QL; do
  [[ -n "${!v:-}" ]] && ENVFWD+=(-e "$v=${!v}")
done
[[ ${#ENVFWD[@]} -gt 0 ]] && echo "env      : ${ENVFWD[*]}"

## every input location in local_paths.R, read-only (files: their folder)
MOUNTS=()
for v in PIPELINE_DIR ORIG_PIPELINE_DIR ORIG_ANALYSIS_DIR FIXED_META_DIR VIRAL_REF_DIR RAW_DIR \
         MSN_REF_ORIG RLIBS_DIR RLIBS_CCC MULTINICHENET_REF_DIR PUB_SUPP_DIR; do
  p="$(lp $v)"; [[ -z "$p" ]] && continue
  [[ -f "$p" ]] && p="$(dirname "$p")"
  [[ -d "$p" ]] && MOUNTS+=(-v "$p:$p:ro")
done

nohup docker run --rm \
  "${ENVFWD[@]+"${ENVFWD[@]}"}" \
  -v "$REPO:$REPO" \
  -v "$OUT_ROOT:$OUT_ROOT" \
  "${MOUNTS[@]}" \
  -w "$REPO" \
  -e SCORCH_REPO="$REPO" \
  --user "$(id -u):$(id -g)" \
  -e HDF5_USE_FILE_LOCKING=FALSE \
  -e TMPDIR="$SCRATCH" -e TMP="$SCRATCH" -e TEMP="$SCRATCH" \
  "$IMAGE" \
  Rscript "$SCRIPT" "$@" \
  > "$LOG" 2>&1 &

echo "launched : $SCRIPT"
echo "log      : $LOG"
echo "pid      : $!"
echo
echo "  tail -f $LOG"
