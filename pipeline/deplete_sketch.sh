#!/bin/bash
# One run: Qiita parquet -> miint FASTQ (stream) -> deacon host depletion -> sylph sketch.
# Usage: deplete_sketch.sh <run_accession> <out_dir> [threads]
# Writes <out_dir>/<study>/<run>/{<run>_1.fq.gz,<run>_2.fq.gz | <run>.fq.gz, deacon.json, <run>.sylsp}
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); source "$P/../config.sh"
SH=$QDEV_SHARE; BIN=$QDEV_CONDA/bin
IDX=$SH/db/panhuman-1.k31w15.idx
RUN=$1; OUT=$2; T=${3:-8}
row=$(awk -F'\t' -v r="$RUN" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $c["run_accession"]==r{print $c["study_accession"]"\t"$c["prep_protocol"]"\t"$c["layout"]"\t"$c["path"]}' $SH/manifest.tsv)
[ -n "$row" ] || { echo "$RUN not in manifest" >&2; exit 2; }
IFS=$'\t' read -r STUDY PROTO LAYOUT PARQ <<<"$row"
[[ $PROTO == *_metagenomics ]] || { echo "$RUN is $PROTO, not metagenomics; skipping" >&2; exit 0; }
D=$OUT/$STUDY/$RUN; mkdir -p "$D"
[ -f "$D/$RUN.sylsp" ] && { echo "$RUN already done"; exit 0; }
STREAM=("$BIN/python" "$P/stream_fastq.py" "$PARQ" --threads 2)
if [ "$LAYOUT" = PAIRED ]; then
  "${STREAM[@]}" | "$BIN/deacon" filter -d --interleaved -t "$T" "$IDX" - \
      -o "$D/${RUN}_1.fq.gz" -O "$D/${RUN}_2.fq.gz" -s "$D/deacon.json"
  "$BIN/sylph" sketch -1 "$D/${RUN}_1.fq.gz" -2 "$D/${RUN}_2.fq.gz" -S "$RUN" -d "$D" -t "$T"
else
  "${STREAM[@]}" | "$BIN/deacon" filter -d -t "$T" "$IDX" - -o "$D/$RUN.fq.gz" -s "$D/deacon.json"
  "$BIN/sylph" sketch -r "$D/$RUN.fq.gz" -S "$RUN" -d "$D" -t "$T"
fi
mv "$D"/*.sylsp "$D/$RUN.sylsp.tmp" && mv "$D/$RUN.sylsp.tmp" "$D/$RUN.sylsp"
