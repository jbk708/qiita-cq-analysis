"""Tests for scope_check.py and native-python.sh (the per-job gate on the shim env var)."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import unittest

import scope_check as sc

try:
    import duckdb
except ImportError:
    duckdb = None

HERE = os.path.dirname(os.path.abspath(__file__))
WRAPPER = os.path.join(HERE, "native-python.sh")
URL = "http://node1:5000/portal/api"

STUB = """#!/bin/bash
if [[ $1 == -I ]]; then echo "$STUB_VERDICT"; [[ $STUB_VERDICT != crash ]]; exit; fi
{ printf '%s\\n' "$@"; echo "URL=${MIINT_ENA_PORTAL_URL_BASE-unset}"; } > "$STUB_OUT"
"""


class TestClassify(unittest.TestCase):
    def test_in_out_mixed(self):
        self.assertEqual(sc.classify({"A", "B"}, {"A", "B", "C"}), "in")
        self.assertEqual(sc.classify({"X"}, {"A"}), "out")
        self.assertEqual(sc.classify({"A", "X"}, {"A"}), "mixed")

    def test_empty_job_is_an_error(self):
        with self.assertRaises(ValueError):
            sc.classify(set(), {"A"})


@unittest.skipIf(duckdb is None, "duckdb not installed")
class TestMain(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = tmp.name
        self.portal_runs = os.path.join(self.dir, "portal.runs")
        with open(self.portal_runs, "w") as f:
            f.write("ERR1\nERR2\n")

    def params(self, runs):
        pq = os.path.join(self.dir, "ena_run_map.parquet")
        con = duckdb.connect()
        con.execute("CREATE TABLE t(prep_sample_idx BIGINT, ena_run_accession VARCHAR)")
        con.executemany("INSERT INTO t VALUES (?, ?)", list(enumerate(runs, 1)))
        con.execute(f"COPY t TO '{pq}' (FORMAT parquet)")
        inp = os.path.join(self.dir, "input")
        os.makedirs(inp, exist_ok=True)
        with open(os.path.join(inp, "params.json"), "w") as f:
            json.dump({"step_name": "ingest_ena_reads", "inputs": {"ena_run_map": pq}}, f)
        return inp

    def test_verdicts(self):
        for runs, want in ((["ERR1", "ERR2"], "in"), (["SRR9"], "out"), (["ERR1", "SRR9"], "mixed")):
            with self.subTest(runs=runs):
                os.environ["QIITA_INPUT_PATH"] = self.params(runs)
                self.addCleanup(os.environ.pop, "QIITA_INPUT_PATH", None)
                self.assertEqual(sc.verdict(self.portal_runs), want)

    def test_missing_params_raises(self):
        os.environ["QIITA_INPUT_PATH"] = self.dir
        self.addCleanup(os.environ.pop, "QIITA_INPUT_PATH", None)
        with self.assertRaises(OSError):
            sc.verdict(self.portal_runs)


class TestWrapper(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = tmp.name
        self.stub = os.path.join(self.dir, "stub-python")
        with open(self.stub, "w") as f:
            f.write(STUB)
        os.chmod(self.stub, 0o755)
        self.portal_env = os.path.join(self.dir, "portal.env")
        self.out = os.path.join(self.dir, "out.txt")

    def run_wrapper(self, *argv, verdict="in", env_file=URL, runs=True):
        if env_file is not None:
            with open(self.portal_env, "w") as f:
                f.write(f"MIINT_ENA_PORTAL_URL_BASE={env_file}\n")
        if runs:
            with open(os.path.join(self.dir, "portal.runs"), "w") as f:
                f.write("ERR1\n")
        env = {"PATH": os.environ["PATH"], "SHIM_PY": self.stub, "SHIM_STATE_DIR": self.dir,
               "STUB_VERDICT": verdict, "STUB_OUT": self.out}
        p = subprocess.run(["/bin/bash", WRAPPER, *argv], env=env, capture_output=True, text=True)
        if not os.path.exists(self.out):
            return p, None
        with open(self.out) as f:
            return p, f.read().split("\n")

    JOB = ("-m", "qiita_compute_orchestrator.jobs", "--job", "ingest_ena_reads")

    def test_no_portal_env_is_pure_passthrough(self):
        p, got = self.run_wrapper(*self.JOB, env_file=None)
        self.assertEqual(p.returncode, 0)
        self.assertEqual(got[:-2], list(self.JOB))
        self.assertEqual(got[-2], "URL=unset")

    def test_other_job_does_not_export(self):
        p, got = self.run_wrapper("-m", "x", "--job", "ingest_reads")
        self.assertEqual(p.returncode, 0)
        self.assertEqual(got[-2], "URL=unset")
        self.assertNotIn("ena-shim", p.stderr)

    def test_in_scope_exports_and_logs(self):
        p, got = self.run_wrapper(*self.JOB, verdict="in")
        self.assertEqual(p.returncode, 0)
        self.assertEqual(got[:-2], list(self.JOB))
        self.assertEqual(got[-2], f"URL={URL}")
        self.assertIn(f"ena-shim: portal={URL}", p.stderr)

    def test_out_of_scope_passes_through(self):
        p, got = self.run_wrapper(*self.JOB, verdict="out")
        self.assertEqual(p.returncode, 0)
        self.assertEqual(got[-2], "URL=unset")

    def test_mixed_scope_exits_1(self):
        p, got = self.run_wrapper(*self.JOB, verdict="mixed")
        self.assertEqual((p.returncode, got), (1, None))

    def test_crashed_scope_check_exits_1(self):
        p, got = self.run_wrapper(*self.JOB, verdict="crash")
        self.assertEqual((p.returncode, got), (1, None))

    def test_non_http_value_exits_1(self):
        p, got = self.run_wrapper(*self.JOB, env_file="https://node1:5000/portal/api")
        self.assertEqual((p.returncode, got), (1, None))
        p, got = self.run_wrapper(*self.JOB, env_file="")
        self.assertEqual((p.returncode, got), (1, None))


if __name__ == "__main__":
    unittest.main()
