"""Print in, out or mixed: whether this ingest_ena_reads job's runs are all, none or some of the shim's runs.

Usage: scope_check.py <portal.runs>. Reads the job's params.json from $QIITA_INPUT_PATH; any error exits non-zero.
"""

from __future__ import annotations

import json
import os
import sys


def classify(job_runs, portal_runs):
    if not job_runs:
        raise ValueError("job has no ENA runs")
    hits = job_runs & portal_runs
    return "in" if hits == job_runs else "out" if not hits else "mixed"


def verdict(portal_runs_path):
    import duckdb

    with open(os.path.join(os.environ["QIITA_INPUT_PATH"], "params.json")) as f:
        run_map = json.load(f)["inputs"]["ena_run_map"]
    job_runs = {r for (r,) in duckdb.connect().execute(
        "SELECT ena_run_accession FROM read_parquet(?)", [run_map]).fetchall()}
    with open(portal_runs_path) as f:
        return classify(job_runs, set(f.read().split()))


if __name__ == "__main__":
    print(verdict(sys.argv[1]))
