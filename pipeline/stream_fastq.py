"""Stream one run's reads from its Qiita parquet as FASTQ on stdout (miint writer).

Reads come out in file order (the parquet is written sorted by sequence_idx), so nothing is
sorted in memory. Paired runs are written interleaved (R1, R2, R1, R2, ...) for `deacon filter --interleaved`.
Usage: stream_fastq.py <read.parquet> [--threads N]
"""
import argparse
import duckdb

MIINT_REPO = "https://ftp.microbio.me/pub/miint"

p = argparse.ArgumentParser()
p.add_argument("parquet")
p.add_argument("--threads", type=int, default=4)
p.add_argument("--memory", default="4GB", help="DuckDB memory_limit; DuckDB sizes to the node, not the Slurm allocation")
a = p.parse_args()

con = duckdb.connect(config={"allow_unsigned_extensions": "true", "threads": a.threads,
                             "memory_limit": a.memory, "preserve_insertion_order": "true"})
con.execute(f"INSTALL miint FROM '{MIINT_REPO}'; LOAD miint")
src = f"read_parquet('{a.parquet}')"
(paired,) = con.execute(f"SELECT count(sequence2) > 0 FROM (SELECT sequence2 FROM {src} LIMIT 1)").fetchone()
opts = "FORMAT FASTQ, INTERLEAVE true" if paired else "FORMAT FASTQ"
con.execute(
    f"COPY (SELECT read_id, sequence1, qual1, sequence2, qual2 FROM {src}) "
    f"TO '/dev/stdout' ({opts})"
)
