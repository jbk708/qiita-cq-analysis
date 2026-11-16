"""Stage runs from NCBI's SRA bucket (miint read_ncbi_sra_sequences) as gzipped FASTQ plus a shim manifest per study.

Usage: stage_sra.py --duckdb <miint duckdb CLI> --ext-dir <dir> --out <dir> --manifest-dir <dir> <study>:<run list> ...
Each study gets <manifest-dir>/<study>.tsv once all its listed runs stage and match ENA's read_count; a study with a
mismatch or a run SRA can't serve gets <study>.failed instead and stays on ENA. A study with either file is skipped.
--partial keeps the runs that stage and lists the rest in <study>.skipped (for a proxy-mode shim, which sends them to ENA).
--expected-counts takes NCBI runinfo CSV (Run, spots, spots_with_mates) for runs ENA never processed (read_count 0):
those runs are checked against NCBI's count instead (spots_with_mates for paired runs, as only full pairs are kept).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import os
import re
import subprocess
import sys
import urllib.parse
import urllib.request

COLUMNS = [
    "run_accession", "study_accession", "sample_accession", "experiment_accession", "library_layout",
    "r1_path", "r1_bytes", "r1_md5", "r2_path", "r2_bytes", "r2_md5",
]
FIELDS = "run_accession,study_accession,sample_accession,experiment_accession,library_layout,read_count"


def ena_runs(study):
    q = urllib.parse.urlencode({"accession": study, "result": "read_run", "fields": FIELDS, "format": "tsv", "limit": 0})
    with urllib.request.urlopen(f"https://www.ebi.ac.uk/ena/portal/api/filereport?{q}", timeout=300) as r:
        return list(csv.DictReader(io.StringIO(r.read().decode()), delimiter="\t"))


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def stage(duckdb, ext_dir, out, run):
    acc, paired = run["run_accession"], run["library_layout"] == "PAIRED"
    keep = "WHERE sequence2 IS NOT NULL" if paired else ""
    target = os.path.join(out, f"{acc}_{{ORIENTATION}}.fastq.gz")
    for m in (1, 2):
        stale = os.path.join(out, f"{acc}_R{m}.fastq.gz")
        if os.path.exists(stale):
            os.remove(stale)
    sql = (f"SET extension_directory='{ext_dir}'; LOAD httpfs;"
           f"COPY (SELECT * FROM read_ncbi_sra_sequences('{acc}') {keep}) TO '{target}' (FORMAT FASTQ);"
           f"SELECT count(*) FROM read_fastx('{os.path.join(out, acc + '_R1.fastq.gz')}');")
    res = subprocess.run([duckdb, "-csv", "-noheader", "-c", sql], capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(res.stderr.strip()[-300:])
    reads = int(re.findall(r"\d+", res.stdout)[-1])
    files = [os.path.join(out, f"{acc}_R{m}.fastq.gz") for m in ((1, 2) if paired else (1,))]
    return reads, files


def count_matches(spots, run):
    """ENA's read_count is spots for most runs but reads (two per spot) for some ENA-submitted paired runs."""
    expected = int(run["read_count"])
    return spots == expected or (run["library_layout"] == "PAIRED" and 2 * spots == expected)


def write_tsv(path, rows):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def stage_study(a, study, wanted):
    runs = [r for r in ena_runs(study) if r["run_accession"] in wanted]
    missing = wanted - {r["run_accession"] for r in runs}
    if missing:
        return f"{len(missing)} listed run(s) not in ENA's run report"
    rows, skipped = [], []
    for run in runs:
        acc = run["run_accession"]
        try:
            spots, files = stage(a.duckdb, a.ext_dir, a.out, run)
            if acc in a.expected:
                want = a.expected[acc][run["library_layout"] == "PAIRED"]
                reason = None if spots == want else f"{spots} spots vs expected {want}"
            else:
                reason = None if count_matches(spots, run) else f"{spots} spots vs ENA read_count {run['read_count']}"
        except RuntimeError as exc:
            reason = str(exc)
        if reason:
            if not a.partial:
                return f"{acc}: {reason}"
            for m in (1, 2):
                path = os.path.join(a.out, f"{acc}_R{m}.fastq.gz")
                if os.path.exists(path):
                    os.remove(path)
            skipped.append(f"{acc}: {reason}")
            print(f"{study} {acc} skipped: {reason}", flush=True)
            continue
        print(f"{study} {acc} spots={spots} ena={run['read_count']} ok", flush=True)
        row = {k: run[k] for k in COLUMNS[:5]}
        for m, f in enumerate(files, 1):
            row.update({f"r{m}_path": f, f"r{m}_bytes": os.path.getsize(f), f"r{m}_md5": md5(f)})
        rows.append(row)
    if not rows:
        return f"no run staged; first: {skipped[0] if skipped else 'none listed'}"
    if skipped:
        with open(os.path.join(a.manifest_dir, f"{study}.skipped"), "w") as f:
            f.write("\n".join(skipped) + "\n")
    write_tsv(os.path.join(a.manifest_dir, f"{study}.tsv"), rows)
    return None


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--duckdb", required=True)
    p.add_argument("--ext-dir", required=True, help="DuckDB extension_directory holding httpfs")
    p.add_argument("--out", required=True)
    p.add_argument("--manifest-dir", required=True)
    p.add_argument("--partial", action="store_true", help="stage what SRA can serve; list the rest in <study>.skipped")
    p.add_argument("--expected-counts", help="NCBI runinfo CSV whose spot counts replace ENA's read_count")
    p.add_argument("studies", nargs="+", help="<study>:<file of its run accessions>")
    a = p.parse_args(argv)
    a.expected = {}
    if a.expected_counts:
        with open(a.expected_counts, newline="") as f:
            # (single-end count, paired count), indexed by whether the run is PAIRED
            a.expected = {r["Run"]: (int(r["spots"]), int(r["spots_with_mates"])) for r in csv.DictReader(f)}
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.manifest_dir, exist_ok=True)
    for spec in a.studies:
        study, runs_file = spec.split(":", 1)
        done = [os.path.join(a.manifest_dir, f"{study}{ext}") for ext in (".tsv", ".failed")]
        if any(os.path.exists(d) for d in done):
            print(f"{study}: already staged or failed, skipping", flush=True)
            continue
        reason = stage_study(a, study, set(open(runs_file).read().split()))
        if reason:
            with open(done[1], "w") as f:
                f.write(reason + "\n")
            print(f"{study}: FAILED, stays on ENA: {reason}", flush=True)
        else:
            print(f"{study}: staged", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
