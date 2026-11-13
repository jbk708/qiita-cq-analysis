#!/bin/bash
#SBATCH --job-name=ena-shim
#SBATCH -p short
#SBATCH -A knightlab
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH -t 12:00:00
# Usage (from the state dir, for the log): sbatch ena-shim.sbatch <manifest> [strict|proxy]
set -euo pipefail

STATE=${SHIM_STATE_DIR:-/ddn_scratch/$USER/qiita-pilot/ena-shim}
REPO=${PILOT_REPO:-/ddn_scratch/$USER/qiita-pilot/qiita-ena-manifests}
PY=${SHIM_PY:-/ddn_scratch/$USER/qiita-dev/Qiita/qiita-compute-orchestrator/.venv/bin/python}
manifest=${1:?usage: sbatch ena-shim.sbatch <manifest>}

exec "$PY" -I "$REPO/shim/ena_shim.py" --manifest "$manifest" \
    --audit "$STATE/$(basename "$manifest" .manifest.tsv).audit.jsonl" \
    --mode "${2:-strict}" --env-out "$STATE/portal.env" --runs-out "$STATE/portal.runs"
