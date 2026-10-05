#!/bin/bash
# Status of one ENA import batch and its tickets. Usage: watch.sh <batch_idx>
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
B=$QDEV_ROOT; BATCH=${1:?usage: watch.sh <batch_idx>}; source $B/env/ports
J=$(squeue -u $USER -h -n qiita-dev-stack -t R -o %i)
srun --jobid=$J --overlap -n 1 -c 1 --mem=1G apptainer exec $B/postgres_17.sif psql -h localhost -p $PGPORT -d qiita -At -F' | ' \
 -c "select i.ena_study_accession, i.state, i.download_work_ticket_idxs, coalesce(left(i.failure_reason,120),'') from qiita.ena_import_batch_item i where batch_idx=$BATCH order by 1" \
 -c "select work_ticket_idx, state, coalesce(left(failure_reason,160),'') from qiita.work_ticket where work_ticket_idx in (select unnest(download_work_ticket_idxs) from qiita.ena_import_batch_item where batch_idx=$BATCH) order by 1"
du -sh $B/persistent/ducklake/read $B/scratch/staging/reads 2>/dev/null
