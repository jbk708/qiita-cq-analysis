"""Build the shim manifest (v1 TSV) for a study by matching ENA submitted_bytes to local per_sample_FASTQ files."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys

COLUMNS = [
    "run_accession", "study_accession", "sample_accession", "experiment_accession", "library_layout",
    "r1_path", "r1_bytes", "r1_md5", "r2_path", "r2_bytes", "r2_md5",
]
REQUIRED = ["run_accession", "sample_accession", "experiment_accession", "run_alias",
            "library_layout", "submitted_bytes", "submitted_md5"]
ALIAS = re.compile(r"^qiita_ppdid_(\d+):")


class ManifestBuildError(Exception):
    pass


def _sizes(root, pdid, cache):
    if pdid not in cache:
        d = os.path.join(root, "per_sample_FASTQ", pdid)
        try:
            entries = [e for e in os.scandir(d) if e.is_file()]
        except OSError as e:
            raise ManifestBuildError(f"cannot list {d}: {e}")
        cache[pdid] = {e.path: e.stat().st_size for e in entries}
    return cache[pdid]


def _match(run, mate, size, sizes):
    hits = [p for p, s in sizes.items() if s == size]
    if not hits:
        raise ManifestBuildError(f"{run}: no file of {size} bytes for R{mate}")
    if len(hits) > 1:
        raise ManifestBuildError(f"{run}: several files of {size} bytes for R{mate}: {sorted(hits)}")
    if f"_R{mate}_" not in os.path.basename(hits[0]):
        raise ManifestBuildError(f"{run}: R{mate} bytes matched {hits[0]}, which has no _R{mate}_")
    return hits[0]


def _build_row(ena, study, root, cache, claimed):
    run = ena["run_accession"]
    m = ALIAS.match(ena["run_alias"])
    if not m:
        raise ManifestBuildError(f"{run}: run_alias {ena['run_alias']!r} is not qiita_ppdid_<N>:...")
    try:
        sizes = _sizes(root, m.group(1), cache)
    except ManifestBuildError as e:
        raise ManifestBuildError(f"{run}: {e}")
    nbytes = ena["submitted_bytes"].split(";")
    md5s = ena["submitted_md5"].split(";")
    want = 2 if ena["library_layout"] == "PAIRED" else 1
    if len(nbytes) != want or len(md5s) != want:
        raise ManifestBuildError(
            f"{run}: {ena['library_layout']} needs {want} submitted file(s), got "
            f"{len(nbytes)} sizes and {len(md5s)} md5s")
    row = dict.fromkeys(COLUMNS, "")
    row.update(run_accession=run, study_accession=study, sample_accession=ena["sample_accession"],
               experiment_accession=ena["experiment_accession"], library_layout=ena["library_layout"])
    for mate, (n, md5) in enumerate(zip(nbytes, md5s), 1):
        path = _match(run, mate, int(n), sizes)
        if path in claimed:
            raise ManifestBuildError(f"{run}: {path} already claimed by {claimed[path]}")
        claimed[path] = run
        row.update({f"r{mate}_path": path, f"r{mate}_bytes": n, f"r{mate}_md5": md5})
    return row


def build_rows(ena_rows, study, qiita_root):
    """One manifest row per ENA run; raises ManifestBuildError listing every run that failed."""
    cache, claimed, rows, errors = {}, {}, [], []
    for ena in ena_rows:
        missing = [c for c in REQUIRED if c not in ena]
        if missing:
            raise ManifestBuildError(f"ENA TSV is missing columns: {missing}")
        try:
            rows.append(_build_row(ena, study, qiita_root, cache, claimed))
        except (ManifestBuildError, ValueError) as e:
            errors.append(str(e))
    if errors:
        raise ManifestBuildError("\n".join(errors))
    return rows


def write_manifest(fh, rows):
    w = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ena-tsv", required=True, help="pre-fetched read_run TSV")
    ap.add_argument("--qiita-root", required=True, help="mount holding per_sample_FASTQ/")
    ap.add_argument("--out", required=True)
    ap.add_argument("--study", required=True)
    a = ap.parse_args(argv)
    with open(a.ena_tsv, newline="") as f:
        ena_rows = list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    try:
        rows = build_rows(ena_rows, a.study, a.qiita_root)
    except ManifestBuildError as e:
        print(f"build_manifest: {e}", file=sys.stderr)
        return 1
    with open(a.out, "w", newline="") as f:
        write_manifest(f, rows)
    print(f"build_manifest: {len(rows)}/{len(ena_rows)} runs -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
