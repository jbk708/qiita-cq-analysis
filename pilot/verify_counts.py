"""Compare stored read counts per run (staging copy and lake) with ENA read_count. Usage: verify_counts.py PRJ..."""
import os
import sys

import duckdb

root = os.environ.get("QDEV_ROOT", f"/ddn_scratch/{os.environ['USER']}/qiita-pilot/qiita-dev")
share = os.environ.get("QDEV_SHARE", f"/ddn_scratch/{os.environ['USER']}/qiita-pilot/qiita-dev-share")
studies = sys.argv[1:] or sys.exit(__doc__)

c = duckdb.connect(config={"allow_unsigned_extensions": "true", "extension_directory": f"{root}/duckdb-ext",
                           "memory_limit": "8GB", "threads": 8})
c.execute("LOAD miint; LOAD httpfs")
c.execute(f"CREATE TABLE m AS SELECT * FROM read_csv('{share}/manifest.tsv', delim='\t', header=true) "
          "WHERE study_accession IN (SELECT unnest($s))", {"s": studies})
ps = ",".join(str(r[0]) for r in c.execute("SELECT prep_sample_idx FROM m").fetchall()) or "NULL"
c.execute(f"CREATE TABLE lake AS SELECT prep_sample_idx, count(*) n FROM read_parquet('{root}/persistent/ducklake/read/*.parquet') "
          f"WHERE prep_sample_idx IN ({ps}) GROUP BY 1")
rows = [(s, r, n) for s in studies for r, n in c.execute(
    "SELECT run_accession, read_count FROM read_ena(?, fields => 'run_accession,read_count')", [s]).fetchall()]
c.execute("CREATE TABLE ena(study VARCHAR, run VARCHAR, ena BIGINT)")
c.executemany("INSERT INTO ena VALUES (?,?,?)", rows)
c.sql("""SELECT e.study, count(*) ena_runs, count(m.run_accession) imported,
  sum((m.read_count = e.ena)::INT) staged_match, sum((l.n = e.ena)::INT) lake_match
  FROM ena e LEFT JOIN m ON m.run_accession = e.run LEFT JOIN lake l USING (prep_sample_idx) GROUP BY 1 ORDER BY 1""").show()
c.sql("""SELECT e.study, e.run, e.ena, m.read_count staged, l.n lake FROM ena e LEFT JOIN m ON m.run_accession = e.run
  LEFT JOIN lake l USING (prep_sample_idx) WHERE coalesce(l.n, -1) <> e.ena OR coalesce(m.read_count, -1) <> e.ena""").show(max_rows=50)
