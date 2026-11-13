"""Unit tests for stage_sra.py's per-study outcome, with SRA and ENA stubbed out."""

from __future__ import annotations

import os
import tempfile
import unittest
from unittest import mock

import stage_sra

RUNS = [dict(run_accession=r, study_accession="PRJ1", sample_accession="S" + r, experiment_accession="X" + r,
             library_layout="SINGLE", read_count="10") for r in ("SRR1", "SRR2", "SRR3")]


class TestStageStudy(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = tmp.name
        self.out, self.mdir = os.path.join(self.dir, "staged"), os.path.join(self.dir, "manifests")
        self.runs_file = os.path.join(self.dir, "PRJ1.txt")
        with open(self.runs_file, "w") as f:
            f.write("SRR1 SRR2 SRR3\n")
        self.spots = {"SRR1": 10, "SRR2": 4, "SRR3": 10}
        self.raises = set()
        patches = [mock.patch.object(stage_sra, "ena_runs", return_value=RUNS),
                   mock.patch.object(stage_sra, "stage", side_effect=self.fake_stage)]
        for p in patches:
            p.start()
            self.addCleanup(p.stop)

    def fake_stage(self, duckdb, ext_dir, out, run):
        acc = run["run_accession"]
        path = os.path.join(out, f"{acc}_R1.fastq.gz")
        with open(path, "wb") as f:
            f.write(b"x" * self.spots[acc])
        if acc in self.raises:
            raise RuntimeError("HTTP 404")
        return self.spots[acc], [path]

    def main(self, *flags):
        return stage_sra.main(["--duckdb", "d", "--ext-dir", "e", "--out", self.out, "--manifest-dir", self.mdir,
                               *flags, f"PRJ1:{self.runs_file}"])

    def manifest_runs(self):
        with open(os.path.join(self.mdir, "PRJ1.tsv")) as f:
            return [ln.split("\t")[0] for ln in f.read().splitlines()[1:]]

    def test_default_fails_the_whole_study_on_one_bad_run(self):
        self.main()
        self.assertFalse(os.path.exists(os.path.join(self.mdir, "PRJ1.tsv")))
        with open(os.path.join(self.mdir, "PRJ1.failed")) as f:
            self.assertIn("SRR2: 4 spots vs ENA read_count 10", f.read())

    def test_partial_stages_good_runs_and_lists_skipped(self):
        self.raises.add("SRR3")
        self.main("--partial")
        self.assertEqual(self.manifest_runs(), ["SRR1"])
        with open(os.path.join(self.mdir, "PRJ1.skipped")) as f:
            skipped = f.read().splitlines()
        self.assertEqual([s.split(":")[0] for s in skipped], ["SRR2", "SRR3"])
        self.assertEqual(sorted(os.listdir(self.out)), ["SRR1_R1.fastq.gz"])

    def test_partial_with_no_good_run_fails_the_study(self):
        self.spots["SRR1"] = self.spots["SRR3"] = 1
        self.main("--partial")
        self.assertFalse(os.path.exists(os.path.join(self.mdir, "PRJ1.tsv")))
        self.assertTrue(os.path.exists(os.path.join(self.mdir, "PRJ1.failed")))


if __name__ == "__main__":
    unittest.main()
