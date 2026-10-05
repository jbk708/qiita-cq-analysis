#!/bin/bash
# Wait for the ENA Portal API to answer, then submit the 16S retry + tranche 2 once.
# Start: nohup setsid ./portal-watch.sh >/dev/null 2>&1 &    Log: logs/portal-watch.log
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
B=$QDEV_ROOT; LOG=$B/logs/portal-watch.log; DONE=$B/runs/16s_t2.submitted
URL='https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=study_accession%3D%22PRJNA1150116%22&fields=study_accession&format=tsv'
echo "$$" > $B/logs/portal-watch.pid
printf "# 16S retry (BioSampleModel hotfix) + tranche 2\nPRJNA1150116\nPRJNA685389\nPRJNA905714\nPRJNA1197483\n" > $B/runs/16s_t2.txt
log() { echo "$(date '+%F %T') $*" >> $LOG; }
log "watcher started (pid $$)"
while [ ! -f $DONE ]; do
  code=$(curl -s -m 30 -o /tmp/portal-watch.$$ -w '%{http_code}' "$URL")
  if [ "$code" = 200 ] && grep -q PRJNA1150116 /tmp/portal-watch.$$; then
    if squeue -u $USER -h -n qiita-dev-stack -t R | grep -q .; then
      log "portal up (200); submitting"
      if $QDEV_CLI submit-ena-import --from-file $B/runs/16s_t2.txt --no-watch >> $LOG 2>&1; then
        touch $DONE; log "submitted"
      else
        log "submit failed; retrying in 10 min"
      fi
    else
      log "portal up but stack job not running; waiting"
    fi
  else
    log "portal $code"
  fi
  rm -f /tmp/portal-watch.$$
  [ -f $DONE ] || sleep 600
done
log "watcher exiting"
