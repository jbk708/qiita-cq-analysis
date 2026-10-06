#!/bin/bash
# Run the Rapid 16S `amplicon` workflow per pool, one at a time, from a plan file.
# Plan lines (tab-separated, '#' comments): study  pool_idx  run_idx  trim  orient_primer(true|false)  primer
# Usage: nohup setsid pilot/amplicon-run.sh <plan.tsv> [sortmerna_reference_idx] >/dev/null 2>&1 &
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
PLAN=${1:?usage: amplicon-run.sh <plan.tsv> [sortmerna_reference_idx]}; REF=${2:-1}
LOG=$QDEV_ROOT/logs/amplicon-run.log
log() { echo "$(date "+%F %T") $*" >> $LOG; }
state() { $QDEV_CLI ticket status "$1" 2>/dev/null | python3 -c "import json,sys; d=json.load(sys.stdin); print(d.get('state','?'), (d.get('failure_reason') or '')[:200])" 2>/dev/null; }
echo $$ > $QDEV_ROOT/logs/amplicon-run.pid
log "started: $PLAN (sortmerna_reference_idx=$REF)"
while IFS=$'\t' read -r study pool run trim orient primer <&3; do   # fd 3: srun below would swallow stdin
  [ -n "$pool" ] || continue
  scope="{\"kind\":\"sequenced_pool\",\"sequenced_pool_idx\":$pool,\"sequencing_run_idx\":$run}"
  ctx="{\"sortmerna_reference_idx\":$REF,\"trim\":$trim,\"orient_primer\":$orient${primer:+,\"primer\":\"$primer\"}}"
  out=$($QDEV_CLI ticket submit --action-id amplicon --action-version 1.0.0 --scope-target-json "$scope" --context-json "$ctx" 2>&1)
  t=$(grep -o '"work_ticket_idx": *[0-9]*' <<<"$out" | head -1 | grep -o '[0-9]*$')
  [ -n "$t" ] || { log "$study pool $pool: submit failed: $(tail -2 <<<"$out" | tr '\n' ' ')"; continue; }
  log "$study pool $pool: ticket $t (trim $trim, orient_primer $orient)"
  while :; do
    read -r s reason < <(state "$t")
    case "$s" in
      completed) log "$study pool $pool: ticket $t completed"; break;;
      failed|cancelled) log "$study pool $pool: ticket $t $s: $reason"; break;;
      *) sleep 120;;
    esac
  done
done 3< <(grep -v '^#' "$PLAN")
log "plan done"
