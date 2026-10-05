#!/bin/bash
# Submit miint_profile.py as one Slurm job. Extra sbatch args via SBATCH_ARGS (e.g. "-p gpu").
# Usage: submit_profile.sh <out_dir> PRJ... [-- miint_profile.py options]
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); source "$P/../config.sh"
: "${QDEV_ACCOUNT:?set QDEV_ACCOUNT to your Slurm account}"
OUT=$1; shift; mkdir -p "$OUT/logs"
PY=$QDEV_ROOT/Qiita/qiita-compute-orchestrator/.venv/bin/python   # DuckDB 1.5.4 + the staged miint build
sbatch --parsable -J miint-sylph -A "$QDEV_ACCOUNT" -p short -c 16 --mem=64G -t 12:00:00 ${SBATCH_ARGS:-} \
  -o "$OUT/logs/%j.log" --export=ALL --wrap "$PY -u $P/miint_profile.py $OUT $*"
