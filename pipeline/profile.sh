#!/bin/bash
# Profile every sketch under <out_dir> against GTDB r232 and add taxonomy.
# Usage: profile.sh <out_dir> [threads]   -> <out_dir>/profile.tsv, <out_dir>/taxprof/*.sylphmpa
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); source "$P/../config.sh"
SH=$QDEV_SHARE; BIN=$QDEV_CONDA/bin
OUT=$1; T=${2:-16}
mapfile -t S < <(find "$OUT" -name '*.sylsp' -not -name '*.tmp')
echo "${#S[@]} sketches"
"$BIN/sylph" profile "$SH/db/gtdb-r232-c200-dbv2.syl2db" "${S[@]}" -t "$T" > "$OUT/profile.tsv"
mkdir -p "$OUT/taxprof"
"$BIN/sylph-tax" --taxonomy-dir "$SH/db/sylph-tax" taxprof "$OUT/profile.tsv" -t GTDB_r232 -o "$OUT/taxprof/"
