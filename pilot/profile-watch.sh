#!/bin/bash
# Rolling sylph profiling: when a study's download tickets have all completed, profile its metagenomic runs that
# are not yet in <out_dir>/profile.tsv (pipeline/submit_profile.sh). One profiling job at a time; a study set is
# submitted at most 3 times.
# Usage: nohup setsid pilot/profile-watch.sh <out_dir> >/dev/null 2>&1 &   (SBATCH_ARGS passes through, e.g. "-p gpu")
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
: "${QDEV_ACCOUNT:?set QDEV_ACCOUNT to your Slurm account}"
OUT=${1:?usage: profile-watch.sh <out_dir>}
B=$QDEV_ROOT; LOG=$B/logs/profile-watch.log; TRIES=$B/logs/profile-watch.tries; source $B/env/ports
PY=$QDEV_CONDA/bin/python
log() { echo "$(date "+%F %T") $*" >> $LOG; }
psql_() {
  local j; j=$(squeue -u $USER -h -n qiita-dev-stack -t R -o %i | head -1)
  [ -n "$j" ] || return 1
  srun --jobid=$j --overlap -n 1 -c 1 --mem=1G apptainer exec $B/postgres_17.sif \
    psql -h localhost -p $PGPORT -d qiita -At -c "$1"
}
touch $TRIES; echo $$ > $B/logs/profile-watch.pid; mkdir -p "$OUT"
log "started (out=$OUT)"
while :; do
  if squeue -u $USER -h -n miint-sylph | grep -q .; then sleep 900; continue; fi
  completed=$(psql_ "
    SELECT DISTINCT i.ena_study_accession FROM qiita.ena_import_batch_item i
    WHERE coalesce(cardinality(i.download_work_ticket_idxs), 0) > 0
      AND NOT EXISTS (SELECT 1 FROM qiita.work_ticket w
                      WHERE w.work_ticket_idx = ANY(i.download_work_ticket_idxs) AND w.state <> 'completed')
    ORDER BY 1") || { sleep 300; continue; }
  [ -n "$completed" ] || { sleep 900; continue; }
  read -r unlisted todo < <($PY - "$QDEV_SHARE/manifest.tsv" "$OUT/profile.tsv" $completed <<'PY'
import os, sys, duckdb
man, prof, studies = sys.argv[1], sys.argv[2], sys.argv[3:]
c = duckdb.connect()
c.execute(f"CREATE TABLE m AS SELECT * FROM read_csv('{man}', delim='\t', header=true)")
listed = {s for (s,) in c.execute("SELECT DISTINCT study_accession FROM m").fetchall()}
done = set()
if os.path.exists(prof):
    done = {r for (r,) in c.execute(f"SELECT DISTINCT run FROM read_csv('{prof}', delim='\t', header=true)").fetchall()}
nohits = os.path.join(os.path.dirname(prof), "no_hits.tsv")
if os.path.exists(nohits):
    done |= {line.split("\t")[1].strip() for line in open(nohits) if line.strip()}
todo = sorted({s for s, r in c.execute("SELECT study_accession, run_accession FROM m WHERE prep_protocol LIKE '%metagenomics' "
                                        "AND study_accession IN (SELECT unnest($s))", {"s": studies}).fetchall() if r not in done})
unlisted = [s for s in studies if s not in listed]
print(",".join(unlisted) or "-", ",".join(todo) or "-")
PY
)
  if [ "$unlisted" != "-" ] && [ -n "$unlisted" ]; then
    if [ "$unlisted" = "${refreshed_for:-}" ]; then
      [ "$unlisted" = "${warned_for:-}" ] || { log "still not in manifest after a refresh: $unlisted; skipping them"; warned_for=$unlisted; }
    else
      log "completed but not in manifest: $unlisted; refreshing share"
      $PILOT_REPO/pilot/share-refresh.sh > $B/logs/share-refresh-profile-watch.log 2>&1 || log "share refresh failed"
      refreshed_for=$unlisted; sleep 60; continue
    fi
  fi
  if [ "$todo" != "-" ] && [ -n "$todo" ]; then
    n=$(grep -c "^$todo$" $TRIES)
    if [ "$n" -ge 3 ]; then
      grep -q "^gave-up $todo$" $TRIES || { echo "gave-up $todo" >> $TRIES; log "gave up on $todo after $n submissions"; }
    else
      echo "$todo" >> $TRIES
      job=$($PILOT_REPO/pipeline/submit_profile.sh "$OUT" ${todo//,/ } 2>&1 | tail -1)
      log "submitted profiling job $job for ${todo//,/ }"
    fi
  fi
  sleep 900
done
