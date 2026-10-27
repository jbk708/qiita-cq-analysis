"""Strict-mode ENA Portal shim and Range file server for local FASTQs listed in a manifest.

Serves `/portal/api/search` (read_run, tsv) and `/f/<run>/<mate>/<basename>` so miint's
`read_ena_sequences` can fetch files that ENA hosts under other names. Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import socket
import sys
import threading
import time
from email.utils import formatdate
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, unquote, urlsplit

REQUIRED = ("run_accession study_accession sample_accession experiment_accession library_layout "
            "r1_path r1_bytes r1_md5 r2_path r2_bytes r2_md5").split()
BATCH_FIELDS = frozenset((
    "run_accession sample_accession experiment_accession study_accession fastq_ftp fastq_aspera fastq_bytes "
    "fastq_md5 library_layout submitted_ftp submitted_aspera submitted_bytes submitted_format").split())
PARAMS = frozenset(("result", "query", "fields", "limit", "format"))
MAX_ACCESSIONS = 50
MD5 = re.compile(r"^[0-9a-fA-F]{32}$")
EQ_QUERY = re.compile(r'^run_accession\s*=\s*"([A-Za-z0-9]+)"$')
IN_QUERY = re.compile(r"^run_accession\s+IN\s*\((.*)\)$")
QUOTED = re.compile(r'^\s*"([A-Za-z0-9]+)"\s*$')
RANGE = re.compile(r"^bytes=(\d*)-(\d*)$")


class ManifestError(Exception):
    pass


class Run:
    def __init__(self, row):
        self.accession = row["run_accession"]
        self.study = row["study_accession"]
        self.sample = row["sample_accession"]
        self.experiment = row["experiment_accession"]
        self.layout = row["library_layout"]
        self.files = [(row[f"r{m}_path"], int(row[f"r{m}_bytes"]), row[f"r{m}_md5"])
                      for m in (1, 2) if row[f"r{m}_path"]]


def _check_row(row, seen_paths):
    run = row["run_accession"]
    paired = bool(row["r2_path"])
    if row["library_layout"] != ("PAIRED" if paired else "SINGLE"):
        raise ManifestError(f"{run}: layout {row['library_layout']} disagrees with r2 presence")
    if not paired and (row["r2_bytes"] or row["r2_md5"]):
        raise ManifestError(f"{run}: r2 bytes/md5 set without r2_path")
    for m in ((1, 2) if paired else (1,)):
        path = row[f"r{m}_path"]
        if not os.path.isabs(path):
            raise ManifestError(f"{run}: r{m}_path {path!r} is not absolute")
        if not path.endswith(".gz"):
            raise ManifestError(f"{run}: r{m}_path {path!r} does not end in .gz")
        if path in seen_paths:
            raise ManifestError(f"{run}: duplicate path {path}")
        seen_paths.add(path)
        if not MD5.match(row[f"r{m}_md5"]):
            raise ManifestError(f"{run}: r{m}_md5 {row[f'r{m}_md5']!r} is not 32 hex chars")
        try:
            size = os.stat(path).st_size
        except OSError as e:
            raise ManifestError(f"{run}: cannot stat {path}: {e}")
        if not os.access(path, os.R_OK):
            raise ManifestError(f"{run}: {path} is not readable")
        if str(size) != row[f"r{m}_bytes"]:
            raise ManifestError(f"{run}: {path} size {size} != manifest {row[f'r{m}_bytes']}")


def load_manifest(path):
    """Parse and validate the manifest; returns {run_accession: Run} or raises ManifestError listing all problems."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise ManifestError(f"manifest is missing columns: {missing}")
        rows = list(reader)
    runs, seen_paths, errors = {}, set(), []
    for row in rows:
        try:
            if row["run_accession"] in runs:
                raise ManifestError(f"duplicate run {row['run_accession']}")
            _check_row(row, seen_paths)
            runs[row["run_accession"]] = Run(row)
        except (ManifestError, ValueError) as e:
            errors.append(str(e))
    if errors:
        raise ManifestError("\n".join(errors))
    return runs


def parse_range(header, size):
    """Return inclusive (start, end) for a single satisfiable byte range, else None."""
    m = RANGE.match(header.strip())
    if not m or (not m.group(1) and not m.group(2)):
        return None
    a, b = m.groups()
    if not a:
        n = int(b)
        return None if n == 0 or size == 0 else (max(size - n, 0), size - 1)
    start, end = int(a), int(b) if b else size - 1
    if start >= size or end < start:
        return None
    return start, min(end, size - 1)


def parse_accessions(query):
    m = EQ_QUERY.match(query)
    if m:
        return [m.group(1)]
    m = IN_QUERY.match(query)
    if not m:
        raise ValueError(f"unsupported query {query!r}; expected run_accession=\"X\" or run_accession IN (...)")
    accs = []
    for part in m.group(1).split(","):
        q = QUOTED.match(part)
        if not q:
            raise ValueError(f"bad accession {part!r} in query")
        accs.append(q.group(1))
    if not accs or len(accs) > MAX_ACCESSIONS:
        raise ValueError(f"need 1..{MAX_ACCESSIONS} accessions, got {len(accs)}")
    return list(dict.fromkeys(accs))


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    server_version = "ena-shim"

    def log_message(self, *args):
        pass

    def do_GET(self):
        self._serve(head=False)

    def do_HEAD(self):
        self._serve(head=True)

    def _serve(self, head):
        t0 = time.monotonic()
        self._status, self._sent, self._full = 0, 0, False
        try:
            url = urlsplit(self.path)
            if url.path == "/healthz":
                self._healthz()
            elif url.path == "/portal/api/search":
                self._search(url.query, head)
            elif url.path.startswith("/f/"):
                self._file(url.path[3:], head)
            else:
                self._text(404, f"shim: unknown path {url.path}")
        except (BrokenPipeError, ConnectionResetError):
            self.close_connection = True
        finally:
            self.server.audit({
                "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "client": self.client_address[0],
                "method": self.command, "path": self.path, "range": self.headers.get("Range"),
                "status": self._status, "bytes": self._sent, "duration": round(time.monotonic() - t0, 4),
                **({"full_get": True} if self._full else {})})

    def _head(self, status, headers):
        self._status = status
        self.send_response(status)
        for k, v in headers.items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.flush()

    def _text(self, status, text, head=False, ctype="text/plain; charset=utf-8", extra=None):
        body = text.encode()
        self._head(status, {"Content-Type": ctype, "Content-Length": str(len(body)), **(extra or {})})
        if not head:
            self.wfile.write(body)
            self._sent = len(body)

    def _healthz(self):
        try:
            os.stat(self.server.mount_root)
        except OSError as e:
            return self._text(503, f"shim: {e}")
        self._text(200, "ok")

    def _search(self, query, head):
        qs = parse_qs(query, keep_blank_values=True)
        try:
            if set(qs) - PARAMS or any(len(v) != 1 for v in qs.values()):
                raise ValueError("unsupported or repeated parameters")
            if qs.get("result") != ["read_run"] or qs.get("format") != ["tsv"]:
                raise ValueError("only result=read_run and format=tsv are supported")
            fields = qs["fields"][0].split(",") if "fields" in qs else []
            bad = [f for f in fields if f not in BATCH_FIELDS]
            if not fields or bad:
                raise ValueError(f"unsupported fields {bad or 'none given'}")
            accs = parse_accessions(qs.get("query", [""])[0])
            unmapped = [a for a in accs if a not in self.server.runs]
            if unmapped:
                raise ValueError(f"run(s) not in the shim manifest: {', '.join(unmapped)}")
        except ValueError as e:
            return self._text(422, f"shim: {e}", head)
        host = self.headers.get("Host") or f"{socket.gethostname()}:{self.server.server_address[1]}"
        lines = ["\t".join(fields)]
        lines += ["\t".join(self._value(self.server.runs[a], f, host) for f in fields) for a in accs]
        self._text(200, "\n".join(lines) + "\n", head, "text/tab-separated-values; charset=utf-8")

    @staticmethod
    def _value(run, field, host):
        if field == "fastq_ftp":
            return ";".join(f"http://{host}/f/{run.accession}/{i}/{os.path.basename(p)}"
                            for i, (p, _, _) in enumerate(run.files, 1))
        if field == "fastq_bytes":
            return ";".join(str(n) for _, n, _ in run.files)
        if field == "fastq_md5":
            return ";".join(m for _, _, m in run.files)
        return {"run_accession": run.accession, "sample_accession": run.sample, "study_accession": run.study,
                "experiment_accession": run.experiment, "library_layout": run.layout}.get(field, "")

    def _file(self, rest, head):
        parts = [unquote(p) for p in rest.split("/")]
        run = self.server.runs.get(parts[0]) if len(parts) == 3 else None
        if run is None or parts[1] not in ("1", "2") or int(parts[1]) > len(run.files) \
                or "/" in parts[2] or os.path.basename(run.files[int(parts[1]) - 1][0]) != parts[2]:
            return self._text(404, f"shim: no such file {rest}", head)
        path, size, _ = run.files[int(parts[1]) - 1]
        try:
            st = os.stat(path)
        except OSError as e:
            return self._text(500, f"shim: cannot stat {path}: {e}", head)
        if st.st_size != size:
            return self._text(500, f"shim: {path} size {st.st_size} != manifest {size}", head)
        info = {
            "Accept-Ranges": "bytes", "Content-Type": "application/octet-stream",
            "ETag": f'"{st.st_ino:x}-{st.st_size:x}-{st.st_mtime_ns:x}"',
            "Last-Modified": formatdate(st.st_mtime, usegmt=True),
        }
        rng = self.headers.get("Range")
        if head or rng is None:
            self._full = not head
            self._head(200, {**info, "Content-Length": str(size)})
            start, count = 0, 0 if head else size
        else:
            span = parse_range(rng, size)
            if span is None:
                return self._text(416, f"shim: unsatisfiable range {rng!r}", extra={"Content-Range": f"bytes */{size}"})
            start, count = span[0], span[1] - span[0] + 1
            self._head(206, {**info, "Content-Length": str(count), "Content-Range": f"bytes {span[0]}-{span[1]}/{size}"})
        if count:
            with open(path, "rb") as f:
                self._sent = self.connection.sendfile(f, offset=start, count=count)


class Server(ThreadingHTTPServer):
    daemon_threads = True
    request_queue_size = 64

    def __init__(self, addr, runs, audit_path, mount_root):
        super().__init__(addr, Handler)
        self.runs, self.mount_root = runs, mount_root
        self._audit_fh = open(audit_path, "a", buffering=1)
        self._audit_lock = threading.Lock()

    def audit(self, record):
        with self._audit_lock:
            self._audit_fh.write(json.dumps(record) + "\n")

    def server_close(self):
        super().server_close()
        self._audit_fh.close()


def make_server(runs, audit_path, mount_root, bind="127.0.0.1", port=0):
    return Server((bind, port), runs, audit_path, mount_root)


def _write_atomic(path, text):
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as f:
        f.write(text)
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--mode", choices=("strict", "proxy"), default="strict")
    ap.add_argument("--mount-root", default="/qmounts/qiita_data", help="stat'ed by /healthz")
    ap.add_argument("--audit", default="ena-shim.audit.jsonl")
    ap.add_argument("--bind", default="0.0.0.0")
    ap.add_argument("--port", type=int, default=0)
    ap.add_argument("--env-out", help="write MIINT_ENA_PORTAL_URL_BASE=... here once listening")
    ap.add_argument("--runs-out", help="write the manifest's run accessions here, one per line")
    a = ap.parse_args(argv)
    if a.mode == "proxy":
        print("ena_shim: proxy mode is not implemented", file=sys.stderr)
        return 2
    try:
        runs = load_manifest(a.manifest)
    except ManifestError as e:
        print(f"ena_shim: refusing to start:\n{e}", file=sys.stderr)
        return 2
    server = make_server(runs, a.audit, a.mount_root, a.bind, a.port)
    port = server.server_address[1]
    # Left in place on exit: a stale pointer fails jobs loudly, deleting it would route them to real ENA.
    if a.runs_out:
        _write_atomic(a.runs_out, "".join(f"{r}\n" for r in sorted(runs)))
    if a.env_out:
        host = socket.gethostname().split(".")[0]
        _write_atomic(a.env_out, f"MIINT_ENA_PORTAL_URL_BASE=http://{host}:{port}/portal/api\n")
    print(f"ena_shim: serving {len(runs)} runs on port {port}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
