#!/bin/bash
# Rebuild the shared read manifest and study index (studies.tsv), and make new per-sample parquet group-readable.
# Run from the login node after each wave: ./share-refresh.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
set -euo pipefail
B=$QDEV_ROOT; SH=$QDEV_SHARE; source $B/env/ports; R=$B/scratch/staging/reads
J=$(squeue -u $USER -h -n qiita-dev-stack -t R -o %i | head -1); [ -n "$J" ] || { echo "stack not running"; exit 1; }
mkdir -p $SH; chmod 750 $SH
chmod g+x $B $B/scratch $B/scratch/staging
chmod -R g+rX $R
srun --jobid=$J --overlap -n 1 -c 1 --mem=1G apptainer exec -B $SH $B/postgres_17.sif \
  psql -h localhost -p $PGPORT -d qiita -At -F $'\t' -o $SH/.runs.tsv -c "
  SELECT st.bioproject_accession, st.ena_study_accession, st.idx, replace(st.title, E'\t', ' '), ss.ena_run_accession, ss.ena_experiment_accession, ps.idx,
         pp.name, sr.platform
  FROM qiita.sequenced_sample ss
  JOIN qiita.prep_sample ps ON ps.idx = ss.prep_sample_idx
  JOIN qiita.prep_protocol pp ON pp.idx = ps.prep_protocol_idx
  JOIN qiita.sequenced_pool sp ON sp.idx = ss.sequenced_pool_idx
  JOIN qiita.sequencing_run sr ON sr.idx = sp.sequencing_run_idx
  JOIN qiita.prep_sample_to_study pst ON pst.prep_sample_idx = ps.idx
  JOIN qiita.study st ON st.idx = pst.study_idx
  WHERE ss.ena_run_accession IS NOT NULL"
srun --jobid=$J --overlap -n 1 -c 1 --mem=1G apptainer exec -B $SH $B/postgres_17.sif \
  psql -h localhost -p $PGPORT -d qiita -At -F $'\t' -o $SH/.studies.tsv -c "
  SELECT st.idx, st.bioproject_accession, st.ena_study_accession, replace(st.title, E'\t', ' '), st.default_tier,
         string_agg(u.email || ':' || sa.access_tier, ',' ORDER BY u.email)
  FROM qiita.study st
  LEFT JOIN qiita.study_access sa ON sa.study_idx = st.idx
  LEFT JOIN qiita.user u ON u.principal_idx = sa.principal_idx
  GROUP BY st.idx ORDER BY st.idx"
srun --jobid=$J --overlap -n 1 -c 4 --mem=8G $B/Qiita/qiita-compute-orchestrator/.venv/bin/python - $SH $R <<'PY'
import sys, os, duckdb
sh, root = sys.argv[1], sys.argv[2]
c = duckdb.connect()
c.execute(f"""CREATE TABLE runs AS SELECT column0 study_accession, column1 secondary_study_accession, column2::BIGINT qiita_study_idx,
  column3 study_title, column4 run_accession, column5 experiment_accession, column6::BIGINT prep_sample_idx,
  column7 prep_protocol, column8 platform
  FROM read_csv('{sh}/.runs.tsv', delim='\t', header=false, all_varchar=true)""")
rows = []
for (ps,) in c.execute("SELECT prep_sample_idx FROM runs").fetchall():
    p = f"{root}/{ps}/read.parquet"
    if os.path.exists(p):
        n, paired = c.execute(f"SELECT count(*), count(sequence2) FROM read_parquet('{p}')").fetchone()
        rows.append((ps, p, n, "PAIRED" if paired else "SINGLE", os.path.getsize(p)))
c.execute("CREATE TABLE files(prep_sample_idx BIGINT, path VARCHAR, read_count BIGINT, layout VARCHAR, parquet_bytes BIGINT)")
c.executemany("INSERT INTO files VALUES (?,?,?,?,?)", rows)
c.execute(f"""COPY (SELECT r.*, f.layout, f.read_count, f.parquet_bytes, f.path FROM runs r JOIN files f USING (prep_sample_idx)
  ORDER BY study_accession, run_accession) TO '{sh}/manifest.tsv' (DELIMITER '\t', HEADER)""")
c.execute(f"""COPY (SELECT s.column0::BIGINT qiita_study_idx, s.column1 study_accession, s.column2 secondary_study_accession,
  s.column3 study_title, s.column4 default_tier, s.column5 access, count(m.run_accession) runs_with_reads
  FROM read_csv('{sh}/.studies.tsv', delim='\t', header=false, all_varchar=true) s
  LEFT JOIN read_csv('{sh}/manifest.tsv', delim='\t', header=true) m ON m.qiita_study_idx = s.column0::BIGINT
  GROUP BY ALL ORDER BY 1) TO '{sh}/studies.tsv' (DELIMITER '\t', HEADER)""")
print(c.sql(f"""SELECT study_accession, prep_protocol, platform, layout, count(*) runs, sum(read_count) reads,
  round(sum(parquet_bytes)/1e9, 2) gb FROM read_csv('{sh}/manifest.tsv') GROUP BY ALL ORDER BY 1, 2"""))
print(c.sql(f"""SELECT prep_protocol, count(DISTINCT study_accession) studies, count(*) runs, sum(read_count) reads,
  round(sum(parquet_bytes)/1e9, 2) gb FROM read_csv('{sh}/manifest.tsv') GROUP BY 1 ORDER BY 1"""))
PY
rm -f $SH/.runs.tsv $SH/.studies.tsv
rsync -a --delete $PILOT_REPO/pipeline/ $SH/pipeline/ \
  && sed -e "s|{{QDEV_SHARE}}|$SH|g" -e "s|{{QDEV_CONDA}}|$QDEV_CONDA|g" $PILOT_REPO/docs/reads-share.md > $SH/README.md
printf 'export QDEV_SHARE=%s\nexport QDEV_CONDA=%s\nexport QDEV_ACCOUNT=%s\n' "$SH" "$QDEV_CONDA" "$QDEV_ACCOUNT" > $SH/config.sh
chmod -R g+rX $SH
