#!/bin/bash
# Redrive download tickets that failed on a transient ENA/network error, in rounds, for every batch at or
# above <first_batch_idx>. Each ticket is redriven at most MAX_ROUNDS times; other failures are left alone.
# Usage: nohup setsid pilot/auto-redrive.sh <first_batch_idx> >/dev/null 2>&1 &
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
B=$QDEV_ROOT; LOG=$B/logs/auto-redrive.log; ROUNDS=$B/logs/auto-redrive.rounds; source $B/env/ports
FIRST=${1:?usage: auto-redrive.sh <first_batch_idx>}; MAX_ROUNDS=${MAX_ROUNDS:-6}
log() { echo "$(date "+%F %T") $*" >> $LOG; }
psql_() {
  local j; j=$(squeue -u $USER -h -n qiita-dev-stack -t R -o %i | head -1)
  [ -n "$j" ] || return 1
  srun --jobid=$j --overlap -n 1 -c 1 --mem=1G apptainer exec $B/postgres_17.sif \
    psql -h localhost -p $PGPORT -d qiita -At -F' ' -c "$1"
}
touch $ROUNDS; echo $$ > $B/logs/auto-redrive.pid
log "started (batches >= $FIRST, max $MAX_ROUNDS rounds per ticket)"
while :; do
  tickets=$(psql_ "
    SELECT DISTINCT w.work_ticket_idx FROM qiita.ena_import_batch_item i
    JOIN qiita.work_ticket w ON w.work_ticket_idx = ANY(i.download_work_ticket_idxs)
    WHERE i.batch_idx >= $FIRST AND w.state = 'failed'
      AND (w.failure_reason LIKE '%Could not connect%' OR w.failure_reason LIKE '%Could not set lock%'
           OR w.failure_reason LIKE '%data is missing or partial%')
    ORDER BY 1") || { sleep 300; continue; }
  for t in $tickets; do
    n=$(grep -c "^$t$" $ROUNDS)
    if [ "$n" -ge "$MAX_ROUNDS" ]; then
      grep -q "^gave-up $t$" $ROUNDS || { echo "gave-up $t" >> $ROUNDS; log "ticket $t: gave up after $n rounds"; }
      continue
    fi
    echo "$t" >> $ROUNDS
    log "ticket $t: redrive round $((n + 1)) -> $($QDEV_CLI ticket run $t 2>&1 | tr -d ' \n')"
  done
  sleep 900
done
