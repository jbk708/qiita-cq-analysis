#!/bin/bash
# Wait for each download ticket to finish, then redrive the next. Stops on a failure.
# Usage: nohup setsid pilot/redrive-seq.sh 10 8 9 7 >/dev/null 2>&1 &   (the first idx is assumed already redriven)
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
B=$QDEV_ROOT; LOG=$B/logs/redrive-seq.log
log() { echo "$(date "+%F %T") $*" >> $LOG; }
state() { $QDEV_CLI ticket status "$1" 2>/dev/null | python3 -c "import json,sys; print(json.load(sys.stdin).get(\"state\",\"?\"))" 2>/dev/null; }
echo $$ > $B/logs/redrive-seq.pid
first=1
for t in "$@"; do
  if [ $first = 0 ]; then log "redrive $t"; $QDEV_CLI ticket run "$t" >> $LOG 2>&1; fi
  first=0
  while :; do
    s=$(state "$t"); case "$s" in
      completed) log "ticket $t completed"; break;;
      failed|cancelled) log "ticket $t $s; stopping"; exit 1;;
      *) sleep 120;;
    esac
  done
done
log "all done"
