"""Fetch one ENA run through ingest_ena_reads' own helper and print the outcome as JSON.

Run with the orchestrator venv python; honors MIINT_ENA_PORTAL_URL_BASE.
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import duckdb
from qiita_compute_orchestrator.jobs.ingest_ena_reads import _classify_ena_fetch_error, _stage_run_reads


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run")
    ap.add_argument("--tmp-dir", type=Path)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory(dir=a.tmp_dir) as tmp:
        try:
            result = _stage_run_reads(a.run, "http", Path(tmp) / "reads.parquet", Path(tmp), 4, 2)
            count, warnings = result if isinstance(result, tuple) else (result, [])
            out = {"run": a.run, "count": count, "warnings": warnings}
        except duckdb.Error as exc:
            failure = _classify_ena_fetch_error(a.run, exc, step_name="ingest_ena_reads")
            out = {"run": a.run, "error": str(exc), "kind": failure.kind.name}
    print(json.dumps(out))


if __name__ == "__main__":
    main()
