"""End-to-end: real miint + httpfs fetch local fixtures through the shim.

Skipped unless SHIM_MIINT_PYTHON (orchestrator venv python) and MIINT_EXTENSION_DIRECTORY are set.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest

import build_manifest as bm

HERE = os.path.dirname(os.path.abspath(__file__))
PY = os.environ.get("SHIM_MIINT_PYTHON")
EXT = os.environ.get("MIINT_EXTENSION_DIRECTORY")
SMALL_PAIRS = 200
BIG_READS = 70000


def fastq_text(n, mate, rng, length=150):
    recs = []
    for i in range(n):
        seq = bytes(rng.choices(b"ACGT", k=length))
        qual = bytes(rng.choices(range(35, 74), k=length))
        recs.append(b"@M1:1:FC:1:1:%d:1 %d:N:0:1\n%s\n+\n%s\n" % (i, mate, seq, qual))
    return b"".join(recs)


def write_gz(path, data, level=1):
    with gzip.open(path, "wb", compresslevel=level) as f:
        f.write(data)


def flip_first_base(data):
    i = data.index(b"\n") + 1
    return data[:i] + (b"C" if data[i:i + 1] == b"A" else b"A") + data[i + 1:]


def md5(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


@unittest.skipUnless(PY and EXT, "needs SHIM_MIINT_PYTHON and MIINT_EXTENSION_DIRECTORY")
class TestMiintThroughShim(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base = os.environ.get("SHIM_TEST_DIR")
        if base:
            os.makedirs(base, exist_ok=True)
        cls.tmp = tempfile.TemporaryDirectory(dir=base)
        d, rng = cls.tmp.name, random.Random(7)
        files, bad, ref = {}, {}, {}
        for name, n, mate in (("small_R1", SMALL_PAIRS, 1), ("small_R2", SMALL_PAIRS, 2), ("big_R1", BIG_READS, 1)):
            files[name] = os.path.join(d, f"{name}_001.fastq.gz")
            write_gz(files[name], fastq_text(n, mate, rng))
        # Stored (level 0) gzip keeps a one-base change the same size and the CRC valid, so only md5 can catch it.
        for mate in (1, 2):
            text = fastq_text(SMALL_PAIRS, mate, rng)
            ref[mate], bad[mate] = (os.path.join(d, f"{k}_R{mate}_001.fastq.gz") for k in ("ref", "bad"))
            write_gz(ref[mate], text, 0)
            write_gz(bad[mate], flip_first_base(text), 0)
            assert os.path.getsize(ref[mate]) == os.path.getsize(bad[mate])
        rows = []
        for run, layout, r1, r2, md5s in (
            ("ERR9000001", "PAIRED", files["small_R1"], files["small_R2"], None),
            ("ERR9000002", "SINGLE", files["big_R1"], None, None),
            ("ERR9000003", "PAIRED", bad[1], bad[2], (md5(ref[1]), md5(ref[2]))),
        ):
            r = dict(dict.fromkeys(bm.COLUMNS, ""), run_accession=run, study_accession="PRJEB0", sample_accession="SAMEA0",
                     experiment_accession="ERX0", library_layout=layout, r1_path=r1, r1_bytes=os.path.getsize(r1),
                     r1_md5=md5s[0] if md5s else md5(r1))
            if r2:
                r.update(r2_path=r2, r2_bytes=os.path.getsize(r2), r2_md5=md5s[1] if md5s else md5(r2))
            rows.append(r)
        cls.manifest = os.path.join(d, "t.manifest.tsv")
        with open(cls.manifest, "w", newline="") as f:
            bm.write_manifest(f, rows)
        cls.audit = os.path.join(d, "t.audit.jsonl")
        cls.shim = subprocess.Popen(
            [sys.executable, os.path.join(HERE, "ena_shim.py"), "--manifest", cls.manifest, "--mount-root", d,
             "--audit", cls.audit, "--bind", "127.0.0.1"], stdout=subprocess.PIPE, text=True)
        port = cls.shim.stdout.readline().rsplit(" ", 1)[-1].strip()
        if not port.isdigit():
            cls.shim.kill()
            raise RuntimeError("shim did not start")
        cls.url = f"http://127.0.0.1:{port}/portal/api"

    @classmethod
    def tearDownClass(cls):
        cls.shim.terminate()
        cls.shim.wait(10)
        cls.shim.stdout.close()
        cls.tmp.cleanup()

    def probe(self, run):
        env = dict(os.environ, MIINT_ENA_PORTAL_URL_BASE=self.url)
        p = subprocess.run([PY, "-I", os.path.join(HERE, "miint_probe.py"), run, "--tmp-dir", self.tmp.name],
                           env=env, capture_output=True, text=True, timeout=600)
        self.assertEqual(p.returncode, 0, p.stderr)
        return json.loads(p.stdout.strip().splitlines()[-1])

    def audit_lines(self):
        with open(self.audit) as f:
            return [json.loads(line) for line in f]

    def test_happy_paths_use_range_requests_and_no_full_get(self):
        start = len(self.audit_lines())
        small, big = self.probe("ERR9000001"), self.probe("ERR9000002")
        self.assertEqual(small["count"], SMALL_PAIRS, small)
        self.assertEqual(big["count"], BIG_READS, big)
        lines = self.audit_lines()[start:]
        self.assertTrue(any(ln["path"].startswith("/portal/api/search") for ln in lines))
        self.assertTrue(any(ln["method"] == "HEAD" and ln["path"].startswith("/f/") for ln in lines))
        self.assertTrue(any(ln["status"] == 206 for ln in lines))
        self.assertEqual([ln["path"] for ln in lines if ln.get("full_get")], [])

    def test_same_size_corruption_fails_on_md5(self):
        out = self.probe("ERR9000003")
        self.assertIn("md5", out.get("error", "").lower(), out)
        self.assertEqual(out["kind"], "EXTERNAL_FETCH_TRANSIENT")

    def test_unmapped_run_is_http_422(self):
        out = self.probe("ERR9999999")
        self.assertIn("HTTP 422", out.get("error", ""), out)

    # Qiita's classifier matches "md5" in the error text, which echoes the search URL's fastq_md5 field.
    @unittest.expectedFailure
    def test_unmapped_run_is_bad_input(self):
        self.assertEqual(self.probe("ERR9999999")["kind"], "BAD_INPUT")


if __name__ == "__main__":
    unittest.main()
