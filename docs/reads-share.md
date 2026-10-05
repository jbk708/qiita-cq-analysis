# ENA reads from the Qiita dev stack (barnacle)

Reads imported from ENA by a personal Qiita stack, one Parquet file per sequencing run.
Public ENA data. Read-only for the owner's Unix group. The set grows as imports land.

- `manifest.tsv` — one row per run: `study_accession` (PRJ…), `secondary_study_accession`
  (SRP/ERP/DRP), `run_accession`, `experiment_accession`, `prep_sample_idx`, `prep_protocol`,
  `platform`, `layout` (SINGLE/PAIRED), `read_count` (equals ENA's), `parquet_bytes`, `path`.
- Each `path` is a Parquet file with columns `sequence_idx`, `read_id`, `sequence1`, `qual1`,
  `sequence2`, `qual2` (the `*2` columns are NULL for single-end runs).

## pandas

```python
import pandas as pd

m = pd.read_csv("{{QDEV_SHARE}}/manifest.tsv", sep="\t")
run = m[m.run_accession == "SRR8063290"].iloc[0]
reads = pd.read_parquet(run.path, columns=["read_id", "sequence1", "qual1", "sequence2", "qual2"])
```

Needs `pyarrow`. Large runs are multi-GB, so read only the columns you need, or use DuckDB.

## DuckDB (with or without miint)

```python
import duckdb

con = duckdb.connect()
m = con.sql("SELECT * FROM read_csv('{{QDEV_SHARE}}/manifest.tsv')")
paths = [p for (p,) in con.sql("SELECT path FROM m WHERE study_accession = 'PRJNA1083602'").fetchall()]
con.sql(f"SELECT count(*) FROM read_parquet({paths})").show()      # all runs of a study at once
```

Write FASTQ with miint's writer (`{ORIENTATION}` gives `_R1`/`_R2` for paired runs):

```python
con.execute("SET allow_unsigned_extensions = true")
con.execute("INSTALL miint FROM 'https://ftp.microbio.me/pub/miint'; LOAD miint")
con.execute(f"""COPY (SELECT read_id, sequence1, qual1, sequence2, qual2
                      FROM read_parquet('{paths[0]}') ORDER BY sequence_idx)
               TO 'SRR_out_{{ORIENTATION}}.fastq.gz' (FORMAT FASTQ)""")
```

## Notebook and conda env

- Env: `{{QDEV_CONDA}}` (DuckDB 1.5.4 to match miint, pandas, pyarrow,
  JupyterLab, deacon 0.18, sylph 1.0, sylph-tax). Register it once:
  `{{QDEV_CONDA}}/bin/python -m ipykernel install --user --name qiita-reads`
- Notebook: `qiita_reads.ipynb` covers the manifest, reading a run, miint QC, FASTQ export, and loading profiles.
- JupyterLab on the `jupyter` partition: `sbatch jupyter.sbatch`, then follow the tunnel line in `jupyter-<job>.log`.

## Host depletion + profiling (metagenomic runs only)

```
P={{QDEV_SHARE}}/pipeline
$P/submit_runs.sh /ddn_scratch/<you>/profiling/<name> PRJNA493153 [more PRJ...]   # Slurm array, 12 at a time
$P/profile.sh     /ddn_scratch/<you>/profiling/<name>                              # after the array; needs ~24 GB
```

Each run streams from its parquet as FASTQ (miint), goes through `deacon filter -d` against `panhuman-1`
(depleted FASTQ plus `deacon.json`), then `sylph sketch`. `profile.sh` runs one `sylph profile` against
GTDB r232 over every sketch, giving `profile.tsv`, then `sylph-tax`, giving `taxprof/<run>.sylphmpa`.
Databases are under `db/`. Runs are skipped unless their protocol is `*_metagenomics`; a run that already
has a sketch is skipped, so re-submitting only redoes failures.
