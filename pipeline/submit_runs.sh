#!/bin/bash
# Submit deplete+sketch as a Slurm array over the metagenomic runs of one or more studies.
# Usage: submit_runs.sh <out_dir> <study_accession> [<study_accession> ...]
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); source "$P/../config.sh"
SH=$QDEV_SHARE; OUT=$1; shift
mkdir -p "$OUT/logs"
LIST=$OUT/runs-$(date +%Y%m%d-%H%M%S).txt
awk -F'\t' -v studies=" $* " 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
  index(studies, " "$c["study_accession"]" ") && $c["prep_protocol"] ~ /_metagenomics$/ {print $c["run_accession"]}' \
  $SH/manifest.tsv > "$LIST"
N=$(wc -l < "$LIST"); [ "$N" -gt 0 ] || { echo "no metagenomic runs for: $*"; exit 1; }
echo "$N runs -> $LIST"
sbatch --parsable -J deacon-sylph -p short -c 8 --mem=16G -t 08:00:00 --array=1-"$N"%12 \
  -o "$OUT/logs/%A_%a.log" \
  --wrap "bash $P/deplete_sketch.sh \$(sed -n \${SLURM_ARRAY_TASK_ID}p $LIST) $OUT 8"
