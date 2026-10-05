#!/bin/bash
# Submit study lists one after another: each waits until every item has registered (or failed) and its
# download tickets are done, then refreshes the share. Failures are logged and do not stop the chain.
# Usage: nohup setsid pilot/wave-chain.sh runs/wave_A.txt runs/wave_B.txt runs/wave_C.txt >/dev/null 2>&1 &
#        pilot/wave-chain.sh --check <batch_idx>   (prints: open completed failed)
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
B=$QDEV_ROOT; LOG=$B/logs/wave-chain.log; source $B/env/ports
log() { echo "$(date "+%F %T") $*" >> $LOG; }
psql_() {
  local j; j=$(squeue -u $USER -h -n qiita-dev-stack -t R -o %i | head -1)
  srun --jobid=$j --overlap -n 1 -c 1 --mem=1G apptainer exec $B/postgres_17.sif \
    psql -h localhost -p $PGPORT -d qiita -At -F' ' -c "$1"
}
counts() {
  psql_ "
    WITH i AS (SELECT * FROM qiita.ena_import_batch_item WHERE batch_idx = $1),
         t AS (SELECT w.state FROM qiita.work_ticket w JOIN i ON w.work_ticket_idx = ANY(i.download_work_ticket_idxs))
    SELECT (SELECT count(*) FROM i WHERE state <> 'failed' AND coalesce(cardinality(download_work_ticket_idxs), 0) = 0)
             + (SELECT count(*) FROM t WHERE state NOT IN ('completed','failed','cancelled')),
           (SELECT count(*) FROM t WHERE state = 'completed'),
           (SELECT count(*) FROM t WHERE state IN ('failed','cancelled'))
             + (SELECT count(*) FROM i WHERE state = 'failed')"
}
[ "$1" = --check ] && { counts "$2"; exit; }

echo $$ > $B/logs/wave-chain.pid
for list in "$@"; do
  log "submit $list"
  out=$($QDEV_CLI submit-ena-import --from-file "$list" --no-watch 2>&1)
  batch=$(grep -o '"ena_import_batch_idx": *[0-9]*' <<<"$out" | head -1 | grep -o '[0-9]*$')
  [ -n "$batch" ] || { log "submit failed: $(tail -3 <<<"$out")"; exit 1; }
  log "batch $batch"
  while :; do
    read -r open done failed < <(counts $batch)
    [ "${open:-1}" = 0 ] && break
    sleep 600
  done
  log "batch $batch done: $done tickets completed, $failed failed"
  psql_ "SELECT i.ena_study_accession, coalesce(w.work_ticket_idx::text, '-'), coalesce(w.state, i.state),
           left(coalesce(w.failure_reason, i.failure_reason, ''), 160)
         FROM qiita.ena_import_batch_item i
         LEFT JOIN qiita.work_ticket w ON w.work_ticket_idx = ANY(i.download_work_ticket_idxs)
         WHERE i.batch_idx = $batch AND (i.state = 'failed' OR w.state IN ('failed','cancelled'))" >> $LOG
  $PILOT_REPO/pilot/share-refresh.sh > $B/logs/share-refresh-$batch.log 2>&1 && log "share refreshed" || log "share refresh failed"
done
log "all lists done"
