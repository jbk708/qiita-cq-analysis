"""Unit tests for build_manifest.py: size matching against a fake per_sample_FASTQ tree."""

from __future__ import annotations

import csv
import io
import os
import tempfile
import unittest

import build_manifest as bm

MD5 = "0123456789abcdef0123456789abcdef"


def ena_row(run, pdid, sizes, layout="PAIRED", **kw):
    row = {
        "run_accession": run,
        "sample_accession": "SAMEA1",
        "experiment_accession": "ERX1",
        "run_alias": f"qiita_ppdid_{pdid}:{run}.s",
        "library_layout": layout,
        "submitted_bytes": ";".join(str(s) for s in sizes),
        "submitted_md5": ";".join([MD5] * len(sizes)),
    }
    row.update(kw)
    return row


class Tree:
    def __init__(self, test):
        self.tmp = tempfile.TemporaryDirectory()
        test.addCleanup(self.tmp.cleanup)
        self.root = self.tmp.name

    def add(self, pdid, name, size):
        d = os.path.join(self.root, "per_sample_FASTQ", str(pdid))
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, name), "wb") as f:
            f.write(b"x" * size)
        return os.path.join(d, name)


class TestBuildRows(unittest.TestCase):
    def setUp(self):
        self.tree = Tree(self)

    def build(self, rows):
        return bm.build_rows(rows, "PRJEB1", self.tree.root)

    def test_unique_size_match_gives_row(self):
        r1 = self.tree.add(7, "a_R1_001.fastq.gz", 10)
        r2 = self.tree.add(7, "a_R2_001.fastq.gz", 11)
        (row,) = self.build([ena_row("ERR1", 7, [10, 11])])
        self.assertEqual(row["run_accession"], "ERR1")
        self.assertEqual(row["study_accession"], "PRJEB1")
        self.assertEqual((row["r1_path"], row["r1_bytes"], row["r1_md5"]), (r1, "10", MD5))
        self.assertEqual((row["r2_path"], row["r2_bytes"], row["r2_md5"]), (r2, "11", MD5))

    def test_single_layout_leaves_r2_empty(self):
        r1 = self.tree.add(7, "a_R1_001.fastq.gz", 10)
        (row,) = self.build([ena_row("ERR1", 7, [10], layout="SINGLE")])
        self.assertEqual(row["r1_path"], r1)
        self.assertEqual((row["r2_path"], row["r2_bytes"], row["r2_md5"]), ("", "", ""))

    def test_header_is_v1_columns(self):
        self.assertEqual(
            bm.COLUMNS,
            "run_accession study_accession sample_accession experiment_accession library_layout "
            "r1_path r1_bytes r1_md5 r2_path r2_bytes r2_md5".split(),
        )
        buf = io.StringIO()
        bm.write_manifest(buf, [])
        self.assertEqual(buf.getvalue().rstrip("\n").split("\t"), bm.COLUMNS)

    def test_zero_candidates_names_run(self):
        self.tree.add(7, "a_R1_001.fastq.gz", 10)
        with self.assertRaisesRegex(bm.ManifestBuildError, "ERR1.*no file"):
            self.build([ena_row("ERR1", 7, [10, 11])])

    def test_size_collision_names_run(self):
        self.tree.add(7, "a_R1_001.fastq.gz", 10)
        self.tree.add(7, "b_R1_001.fastq.gz", 10)
        self.tree.add(7, "a_R2_001.fastq.gz", 11)
        with self.assertRaisesRegex(bm.ManifestBuildError, "ERR1.*several"):
            self.build([ena_row("ERR1", 7, [10, 11])])

    def test_file_claimed_by_two_runs(self):
        self.tree.add(7, "a_R1_001.fastq.gz", 10)
        self.tree.add(7, "a_R2_001.fastq.gz", 11)
        with self.assertRaisesRegex(bm.ManifestBuildError, "ERR2.*already claimed by ERR1"):
            self.build([ena_row("ERR1", 7, [10, 11]), ena_row("ERR2", 7, [10, 11])])

    def test_r1_bytes_landing_on_r2_file(self):
        self.tree.add(7, "a_R2_001.fastq.gz", 10)
        self.tree.add(7, "a_R1_001.fastq.gz", 11)
        with self.assertRaisesRegex(bm.ManifestBuildError, "ERR1.*_R1_"):
            self.build([ena_row("ERR1", 7, [10, 11])])

    def test_non_qiita_alias(self):
        with self.assertRaisesRegex(bm.ManifestBuildError, "ERR1.*run_alias"):
            self.build([ena_row("ERR1", 7, [10, 11], run_alias="something_else")])

    def test_all_errors_reported(self):
        with self.assertRaises(bm.ManifestBuildError) as cm:
            self.build([ena_row("ERR1", 7, [10, 11]), ena_row("ERR2", 7, [12, 13])])
        self.assertIn("ERR1", str(cm.exception))
        self.assertIn("ERR2", str(cm.exception))

    def test_layout_and_size_count_disagree(self):
        self.tree.add(7, "a_R1_001.fastq.gz", 10)
        with self.assertRaisesRegex(bm.ManifestBuildError, "ERR1.*PAIRED"):
            self.build([ena_row("ERR1", 7, [10], layout="PAIRED")])

    def test_missing_ena_column(self):
        row = ena_row("ERR1", 7, [10, 11])
        del row["submitted_md5"]
        with self.assertRaisesRegex(bm.ManifestBuildError, "submitted_md5"):
            self.build([row])


class TestCli(unittest.TestCase):
    def test_main_writes_tsv_and_fails_without_writing(self):
        tree = Tree(self)
        tree.add(7, "a_R1_001.fastq.gz", 10)
        tree.add(7, "a_R2_001.fastq.gz", 11)
        fields = ["run_accession", "sample_accession", "experiment_accession", "run_alias",
                  "library_layout", "submitted_bytes", "submitted_md5"]
        ena = os.path.join(tree.root, "ena.tsv")
        out = os.path.join(tree.root, "m.tsv")

        def write_ena(rows):
            with open(ena, "w", newline="") as f:
                w = csv.DictWriter(f, fields, delimiter="\t")
                w.writeheader()
                w.writerows(rows)

        args = ["--ena-tsv", ena, "--qiita-root", tree.root, "--out", out, "--study", "PRJEB1"]
        write_ena([ena_row("ERR1", 7, [10, 11])])
        self.assertEqual(bm.main(args), 0)
        with open(out) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual([r["run_accession"] for r in rows], ["ERR1"])

        os.remove(out)
        write_ena([ena_row("ERR1", 7, [10, 99])])
        self.assertNotEqual(bm.main(args), 0)
        self.assertFalse(os.path.exists(out))


if __name__ == "__main__":
    unittest.main()
