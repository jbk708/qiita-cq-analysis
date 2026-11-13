"""Unit tests for ena_shim.py: manifest loader, Portal search, Range file server, audit log."""

from __future__ import annotations

import csv
import http.client
import http.server
import json
import os
import re
import tempfile
import threading
import unittest
import urllib.parse

import ena_shim as shim

COLUMNS = ("run_accession study_accession sample_accession experiment_accession library_layout "
           "r1_path r1_bytes r1_md5 r2_path r2_bytes r2_md5").split()
MD5 = "0123456789abcdef0123456789abcdef"
FIELDS = ("run_accession,sample_accession,experiment_accession,study_accession,fastq_ftp,fastq_aspera,"
          "fastq_bytes,fastq_md5,library_layout,submitted_ftp,submitted_aspera,submitted_bytes,submitted_format")


def write_manifest(path, rows, columns=COLUMNS):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


class Env(unittest.TestCase):
    """Two paired runs and one single run backed by real files in a temp dir."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = tmp.name
        self.files = {}
        for name, data in {
            "e1_R1_001.fastq.gz": bytes(range(256)) * 4,
            "e1_R2_001.fastq.gz": b"B" * 300,
            "e2_R1_001.fastq.gz": b"C" * 50,
            "e2_R2_001.fastq.gz": b"D" * 60,
            "s1_R1_001.fastq.gz": b"E" * 70,
        }.items():
            p = os.path.join(self.dir, name)
            with open(p, "wb") as f:
                f.write(data)
            self.files[name] = (p, data)
        self.manifest = os.path.join(self.dir, "m.manifest.tsv")
        self.rows = [
            self.row("ERR1", "e1", paired=True),
            self.row("ERR2", "e2", paired=True),
            self.row("ERR3", "s1", paired=False),
        ]
        write_manifest(self.manifest, self.rows)

    def row(self, run, stem, paired):
        r1 = self.files[f"{stem}_R1_001.fastq.gz"]
        row = dict(run_accession=run, study_accession="PRJEB1", sample_accession="SAME" + run,
                   experiment_accession="ERX" + run, library_layout="PAIRED" if paired else "SINGLE",
                   r1_path=r1[0], r1_bytes=len(r1[1]), r1_md5=MD5, r2_path="", r2_bytes="", r2_md5="")
        if paired:
            r2 = self.files[f"{stem}_R2_001.fastq.gz"]
            row.update(r2_path=r2[0], r2_bytes=len(r2[1]), r2_md5=MD5)
        return row


class TestLoader(Env):
    def load(self, rows, columns=COLUMNS):
        write_manifest(self.manifest, rows, columns)
        return shim.load_manifest(self.manifest)

    def refuses(self, rows, pattern, columns=COLUMNS):
        with self.assertRaisesRegex(shim.ManifestError, pattern):
            self.load(rows, columns)

    def test_loads_valid(self):
        runs = shim.load_manifest(self.manifest)
        self.assertEqual(sorted(runs), ["ERR1", "ERR2", "ERR3"])
        self.assertEqual(len(runs["ERR1"].files), 2)
        self.assertEqual(len(runs["ERR3"].files), 1)

    def test_extra_columns_accepted(self):
        rows = [dict(r, match="size") for r in self.rows]
        self.assertEqual(sorted(self.load(rows, COLUMNS + ["match"])), ["ERR1", "ERR2", "ERR3"])

    def test_missing_column(self):
        self.refuses(self.rows, "r2_md5", [c for c in COLUMNS if c != "r2_md5"])

    def test_duplicate_run(self):
        self.refuses(self.rows + [self.rows[0]], "duplicate run ERR1")

    def test_duplicate_path(self):
        self.refuses([self.rows[0], dict(self.rows[1], r1_path=self.rows[0]["r1_path"],
                                         r1_bytes=self.rows[0]["r1_bytes"])], "duplicate path")

    def test_relative_path(self):
        self.refuses([dict(self.rows[0], r1_path="rel/e1.fastq.gz")], "ERR1.*not absolute")

    def test_md5_must_be_32_hex(self):
        self.refuses([dict(self.rows[0], r1_md5="")], "ERR1.*md5")
        self.refuses([dict(self.rows[0], r2_md5="xyz")], "ERR1.*md5")

    def test_layout_disagrees_with_r2(self):
        self.refuses([dict(self.rows[0], library_layout="SINGLE")], "ERR1.*layout")
        self.refuses([dict(self.rows[2], library_layout="PAIRED")], "ERR3.*layout")

    def test_basename_must_be_gz(self):
        p = os.path.join(self.dir, "plain.fastq")
        with open(p, "wb") as f:
            f.write(b"x")
        self.refuses([dict(self.rows[2], r1_path=p, r1_bytes=1)], "ERR3.*\\.gz")

    def test_missing_file(self):
        os.remove(self.files["e1_R2_001.fastq.gz"][0])
        self.refuses(self.rows, "ERR1.*cannot stat")

    def test_unreadable_file(self):
        if os.geteuid() == 0:
            self.skipTest("root reads everything")
        os.chmod(self.files["e1_R2_001.fastq.gz"][0], 0)
        self.refuses(self.rows, "ERR1.*not readable")

    def test_size_mismatch(self):
        self.refuses([dict(self.rows[0], r1_bytes=1)], "ERR1.*size")


class Served(Env):
    def setUp(self):
        super().setUp()
        self.audit = os.path.join(self.dir, "t.audit.jsonl")
        self.server = shim.make_server(shim.load_manifest(self.manifest), self.audit, self.dir)
        self.port = self.server.server_address[1]
        t = threading.Thread(target=self.server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True)
        t.start()
        self.addCleanup(self.stop)

    def stop(self):
        self.server.shutdown()
        self.server.server_close()

    def req(self, method, path, headers=None, conn=None):
        c = conn or http.client.HTTPConnection("127.0.0.1", self.port, timeout=10)
        c.request(method, path, headers=headers or {})
        r = c.getresponse()
        body = r.read()
        if conn is None:
            c.close()
        return r, body

    def audit_lines(self):
        with open(self.audit) as f:
            return [json.loads(line) for line in f]


class TestSearch(Served):
    def search(self, query, fields=FIELDS, **extra):
        params = dict(result="read_run", query=query, fields=fields, limit="0", format="tsv")
        params.update(extra)
        params = {k: v for k, v in params.items() if v is not None}
        return self.req("GET", "/portal/api/search?" + urllib.parse.urlencode(params, quote_via=urllib.parse.quote))

    def table(self, body):
        lines = body.decode().rstrip("\n").split("\n")
        header = lines[0].split("\t")
        return header, [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]

    def test_equals_form(self):
        r, body = self.search('run_accession="ERR1"')
        self.assertEqual(r.status, 200)
        header, rows = self.table(body)
        self.assertEqual(header, FIELDS.split(","))
        (row,) = rows
        urls = row["fastq_ftp"].split(";")
        self.assertEqual(len(urls), 2)
        self.assertEqual(urls[0], f"http://127.0.0.1:{self.port}/f/ERR1/1/e1_R1_001.fastq.gz")
        self.assertTrue(urls[1].endswith("/f/ERR1/2/e1_R2_001.fastq.gz"))
        self.assertEqual(row["fastq_bytes"], "1024;300")
        self.assertEqual(row["fastq_md5"], f"{MD5};{MD5}")
        self.assertEqual(row["library_layout"], "PAIRED")
        self.assertEqual(row["study_accession"], "PRJEB1")
        for blank in ("fastq_aspera", "submitted_ftp", "submitted_aspera", "submitted_bytes", "submitted_format"):
            self.assertEqual(row[blank], "")

    def test_in_form_and_single_layout(self):
        r, body = self.search('run_accession IN ("ERR3","ERR1")')
        _, rows = self.table(body)
        self.assertEqual([x["run_accession"] for x in rows], ["ERR3", "ERR1"])
        self.assertEqual(rows[0]["fastq_ftp"].count(";"), 0)
        self.assertEqual(rows[0]["fastq_bytes"], "70")

    def test_header_follows_requested_order(self):
        r, body = self.search('run_accession="ERR2"', fields="fastq_md5,run_accession")
        header, rows = self.table(body)
        self.assertEqual(header, ["fastq_md5", "run_accession"])
        self.assertEqual(rows[0]["run_accession"], "ERR2")

    def test_unmapped_run_is_422_naming_run(self):
        r, body = self.search('run_accession="ERR9"')
        self.assertEqual(r.status, 422)
        self.assertIn(b"ERR9", body)

    def test_mixed_mapped_and_unmapped_is_422(self):
        r, body = self.search('run_accession IN ("ERR1","ERR9")')
        self.assertEqual(r.status, 422)
        self.assertIn(b"ERR9", body)

    def test_strictness(self):
        for kw in (
            dict(result="read_experiment"),
            dict(format="json"),
            dict(fields="run_accession,read_count"),
            dict(fields=None),
            dict(unknown="1"),
        ):
            with self.subTest(kw=kw):
                self.assertEqual(self.search('run_accession="ERR1"', **kw)[0].status, 422)
        for q in ('study_accession="PRJEB1"', 'run_accession="ERR1" OR 1=1', "run_accession=ERR1", ""):
            with self.subTest(q=q):
                self.assertEqual(self.search(q)[0].status, 422)

    def test_more_than_50_accessions_is_422(self):
        q = "run_accession IN (" + ",".join(f'"ERR{i}"' for i in range(51)) + ")"
        r, body = self.search(q)
        self.assertEqual(r.status, 422)
        self.assertIn(b"50", body)

    def test_unknown_path_is_404(self):
        self.assertEqual(self.req("GET", "/portal/api/other")[0].status, 404)


class FakeEna(http.server.BaseHTTPRequestHandler):
    """Answers /search from `rows` ({run: {field: value}}); `fail` makes every request a 500."""

    rows, fail, queries = {}, False, []

    def log_message(self, *args):
        pass

    def do_GET(self):
        qs = urllib.parse.parse_qs(urllib.parse.urlsplit(self.path).query)
        FakeEna.queries.append(qs)
        if FakeEna.fail:
            body, status = b"boom", 500
        else:
            fields = qs["fields"][0].split(",")
            accs = re.findall(r'"([A-Za-z0-9]+)"', qs["query"][0])
            lines = ["\t".join(fields)] + ["\t".join(FakeEna.rows[a].get(f, "") for f in fields)
                                           for a in accs if a in FakeEna.rows]
            body, status = ("\n".join(lines) + "\n").encode(), 200
        self.send_response(status)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


class TestProxy(TestSearch):
    def setUp(self):
        FakeEna.rows = {"SRR7": {"run_accession": "SRR7", "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR7_1.fastq.gz",
                                 "fastq_bytes": "11", "fastq_md5": MD5, "library_layout": "SINGLE"}}
        FakeEna.fail, FakeEna.queries = False, []
        self.ena = http.server.ThreadingHTTPServer(("127.0.0.1", 0), FakeEna)
        threading.Thread(target=self.ena.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True).start()
        self.addCleanup(self.ena.server_close)
        self.addCleanup(self.ena.shutdown)
        super().setUp()
        self.stop()
        self.server = shim.make_server(shim.load_manifest(self.manifest), self.audit, self.dir,
                                       ena_portal=f"http://127.0.0.1:{self.ena.server_address[1]}")
        self.port = self.server.server_address[1]
        threading.Thread(target=self.server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True).start()

    def test_unmapped_run_is_422_naming_run(self):
        """Proxy mode: a run neither the shim nor ENA has is a 502 naming it, not a silently short table."""
        r, body = self.search('run_accession="ERR9"')
        self.assertEqual(r.status, 502)
        self.assertIn(b"ERR9", body)

    def test_mixed_mapped_and_unmapped_is_422(self):
        r, body = self.search('run_accession IN ("ERR1","SRR7","ERR9")')
        self.assertEqual(r.status, 502)
        self.assertIn(b"ERR9", body)

    def test_unmapped_runs_come_from_ena_in_request_order(self):
        r, body = self.search('run_accession IN ("SRR7","ERR1")')
        self.assertEqual(r.status, 200)
        header, rows = self.table(body)
        self.assertEqual(header, FIELDS.split(","))
        self.assertEqual([x["run_accession"] for x in rows], ["SRR7", "ERR1"])
        self.assertEqual(rows[0]["fastq_ftp"], "ftp.sra.ebi.ac.uk/vol1/fastq/SRR7_1.fastq.gz")
        self.assertTrue(rows[1]["fastq_ftp"].startswith(f"http://127.0.0.1:{self.port}/f/ERR1/1/"))
        (q,) = FakeEna.queries
        self.assertEqual(q["query"], ['run_accession IN ("SRR7")'])
        self.assertEqual(q["fields"], [FIELDS])
        self.assertEqual(self.audit_lines()[-1]["proxied"], 1)

    def test_mapped_runs_never_reach_ena(self):
        self.assertEqual(self.search('run_accession="ERR1"')[0].status, 200)
        self.assertEqual(FakeEna.queries, [])

    def test_ena_error_is_502(self):
        FakeEna.fail = True
        r, body = self.search('run_accession IN ("SRR7","ERR1")')
        self.assertEqual(r.status, 502)
        self.assertIn(b"500", body)


class TestFiles(Served):
    path1 = "/f/ERR1/1/e1_R1_001.fastq.gz"

    def data1(self):
        return self.files["e1_R1_001.fastq.gz"][1]

    def test_head(self):
        r, body = self.req("HEAD", self.path1)
        self.assertEqual(r.status, 200)
        self.assertEqual(body, b"")
        self.assertEqual(r.getheader("Content-Length"), "1024")
        self.assertEqual(r.getheader("Accept-Ranges"), "bytes")
        self.assertEqual(r.getheader("Content-Type"), "application/octet-stream")
        self.assertIsNone(r.getheader("Content-Encoding"))
        self.assertTrue(r.getheader("ETag"))
        self.assertTrue(r.getheader("Last-Modified"))

    def test_ranges(self):
        d = self.data1()
        cases = {
            "bytes=0-9": (0, 10),
            "bytes=10-": (10, len(d)),
            "bytes=-5": (len(d) - 5, len(d)),
            f"bytes=1000-{len(d) + 500}": (1000, len(d)),
        }
        for hdr, (a, b) in cases.items():
            with self.subTest(hdr=hdr):
                r, body = self.req("GET", self.path1, {"Range": hdr})
                self.assertEqual(r.status, 206)
                self.assertEqual(body, d[a:b])
                self.assertEqual(r.getheader("Content-Range"), f"bytes {a}-{b - 1}/{len(d)}")
                self.assertEqual(r.getheader("Content-Length"), str(b - a))

    def test_unsatisfiable_and_malformed_ranges_are_416(self):
        for hdr in ("bytes=1024-", "bytes=5000-6000", "bytes=0-1,5-9", "bytes=abc", "bytes=-", "bytes=9-3", "items=0-1"):
            with self.subTest(hdr=hdr):
                r, _ = self.req("GET", self.path1, {"Range": hdr})
                self.assertEqual(r.status, 416)
                self.assertEqual(r.getheader("Content-Range"), "bytes */1024")

    def test_no_range_is_200_and_flagged_full_get(self):
        r, body = self.req("GET", self.path1)
        self.assertEqual((r.status, body), (200, self.data1()))
        self.assertTrue(self.audit_lines()[-1]["full_get"])

    def test_range_above_4gib(self):
        p = os.path.join(self.dir, "big_R1_001.fastq.gz")
        off = (4 << 30) + 12345
        with open(p, "wb") as f:
            f.truncate(5 << 30)
            f.seek(off)
            f.write(b"MARKER")
        rows = self.rows + [dict(self.rows[2], run_accession="ERRBIG", r1_path=p, r1_bytes=5 << 30)]
        write_manifest(self.manifest, rows)
        self.stop()
        self.server = shim.make_server(shim.load_manifest(self.manifest), self.audit, self.dir)
        self.port = self.server.server_address[1]
        threading.Thread(target=self.server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True).start()
        self.addCleanup(self.stop)
        r, body = self.req("GET", "/f/ERRBIG/1/big_R1_001.fastq.gz", {"Range": f"bytes={off}-{off + 5}"})
        self.assertEqual((r.status, body), (206, b"MARKER"))
        self.assertEqual(r.getheader("Content-Range"), f"bytes {off}-{off + 5}/{5 << 30}")

    def test_changed_file_is_500(self):
        p = self.files["e1_R1_001.fastq.gz"][0]
        with open(p, "ab") as f:
            f.write(b"more")
        r, body = self.req("GET", self.path1)
        self.assertEqual(r.status, 500)
        self.assertIn(f"shim: {p} size 1028 != manifest 1024".encode(), body)
        with open(p, "r+b") as f:
            f.truncate(10)
        self.assertEqual(self.req("HEAD", self.path1)[0].status, 500)
        os.remove(p)
        self.assertEqual(self.req("GET", self.path1)[0].status, 500)

    def test_404s(self):
        for path in (
            "/f/ERR9/1/e1_R1_001.fastq.gz",
            "/f/ERR1/3/e1_R1_001.fastq.gz",
            "/f/ERR1/2/e1_R1_001.fastq.gz",
            "/f/ERR1/1/other.fastq.gz",
            "/f/ERR3/2/s1_R1_001.fastq.gz",
            "/f/ERR1/1/../e1_R1_001.fastq.gz",
            "/f/ERR1/1/%2e%2e%2fe1_R1_001.fastq.gz",
            "/f/..%2fERR1/1/e1_R1_001.fastq.gz",
            "/f/ERR1/1",
        ):
            with self.subTest(path=path):
                self.assertEqual(self.req("GET", path)[0].status, 404)

    def test_keep_alive_serves_two_requests_on_one_connection(self):
        c = http.client.HTTPConnection("127.0.0.1", self.port, timeout=10)
        self.addCleanup(c.close)
        r1, _ = self.req("HEAD", self.path1, conn=c)
        sock = c.sock
        r2, body = self.req("GET", self.path1, {"Range": "bytes=0-3"}, conn=c)
        self.assertEqual(r1.version, 11)
        self.assertEqual((r2.status, body), (206, self.data1()[:4]))
        self.assertIs(c.sock, sock)

    def test_healthz(self):
        r, body = self.req("GET", "/healthz")
        self.assertEqual((r.status, body), (200, b"ok"))
        os.rename(self.dir, self.dir + ".gone")
        self.addCleanup(os.rename, self.dir + ".gone", self.dir)
        self.assertEqual(self.req("GET", "/healthz")[0].status, 503)

    def test_audit_has_one_line_per_request(self):
        self.req("HEAD", self.path1)
        self.req("GET", self.path1, {"Range": "bytes=0-9"})
        self.req("GET", "/f/ERR9/1/x.gz")
        lines = self.audit_lines()
        self.assertEqual([(ln["method"], ln["status"]) for ln in lines], [("HEAD", 200), ("GET", 206), ("GET", 404)])
        self.assertEqual(lines[1]["range"], "bytes=0-9")
        self.assertEqual(lines[1]["bytes"], 10)
        self.assertEqual(lines[1]["path"], self.path1)
        self.assertNotIn("full_get", lines[1])
        for key in ("time", "client", "duration"):
            self.assertIn(key, lines[0])


class TestMain(Env):
    def test_main_writes_env_and_runs_and_keeps_them(self):
        import subprocess
        import sys
        env_out, runs_out = os.path.join(self.dir, "portal.env"), os.path.join(self.dir, "portal.runs")
        p = subprocess.Popen(
            [sys.executable, os.path.join(os.path.dirname(__file__), "ena_shim.py"), "--manifest", self.manifest,
             "--mount-root", self.dir, "--audit", os.path.join(self.dir, "a.audit.jsonl"),
             "--env-out", env_out, "--runs-out", runs_out, "--bind", "127.0.0.1"],
            stdout=subprocess.PIPE, text=True)
        try:
            p.stdout.readline()
            with open(env_out) as f:
                line = f.read().strip()
            key, _, url = line.partition("=")
            self.assertEqual(key, "MIINT_ENA_PORTAL_URL_BASE")
            self.assertRegex(url, r"^http://[^:/]+:\d+/portal/api$")
            with open(runs_out) as f:
                self.assertEqual(f.read().split(), ["ERR1", "ERR2", "ERR3"])
        finally:
            p.terminate()
            p.stdout.close()
            p.wait(10)
        self.assertTrue(os.path.exists(env_out) and os.path.exists(runs_out))

    def test_main_proxy_mode_marks_env(self):
        import subprocess
        import sys
        env_out = os.path.join(self.dir, "portal.env")
        p = subprocess.Popen(
            [sys.executable, os.path.join(os.path.dirname(__file__), "ena_shim.py"), "--manifest", self.manifest,
             "--mount-root", self.dir, "--audit", os.path.join(self.dir, "a.audit.jsonl"), "--mode", "proxy",
             "--env-out", env_out, "--bind", "127.0.0.1"],
            stdout=subprocess.PIPE, text=True)
        try:
            p.stdout.readline()
            with open(env_out) as f:
                lines = f.read().split()
            self.assertTrue(lines[0].startswith("MIINT_ENA_PORTAL_URL_BASE=http://"))
            self.assertEqual(lines[1], "ENA_SHIM_MODE=proxy")
        finally:
            p.terminate()
            p.stdout.close()
            p.wait(10)

    def test_main_refuses_bad_manifest(self):
        write_manifest(self.manifest, [dict(self.rows[0], r1_md5="")])
        self.assertNotEqual(shim.main(["--manifest", self.manifest, "--mount-root", self.dir]), 0)


if __name__ == "__main__":
    unittest.main()
