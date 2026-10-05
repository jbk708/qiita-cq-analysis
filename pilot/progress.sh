#!/bin/bash
# Per-ticket download progress for the barnacle dev stack. Usage: ./progress.sh  (or: watch -n 60 ./progress.sh)
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
B=$QDEV_ROOT
PY=$B/Qiita/qiita-compute-orchestrator/.venv/bin/python
printf "%-7s %-12s %10s %10s\n" ticket runs_done "/ roster" workspace
for d in $B/scratch/ticket/*/; do
  t=$(basename $d); [ -f $d/ena_run_map.parquet ] || continue
  read total done < <($PY - "$d/ena_run_map.parquet" "$B/scratch/staging/reads" <<'PYEOF'
import sys, os, duckdb
ids = [r[0] for r in duckdb.sql(f"select prep_sample_idx from read_parquet('{sys.argv[1]}')").fetchall()]
root = sys.argv[2]
print(len(ids), sum(os.path.exists(f"{root}/{i}/read.parquet") for i in ids))
PYEOF
)
  printf "%-7s %-12s %10s %10s\n" "$t" "$done" "/ $total" "$(du -sh $d 2>/dev/null | cut -f1)"
done
echo "lake: $(du -sh $B/persistent/ducklake 2>/dev/null | cut -f1)   staged reads: $(du -sh $B/scratch/staging/reads 2>/dev/null | cut -f1)"
