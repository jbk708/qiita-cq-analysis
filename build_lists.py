#!/usr/bin/env python3
"""M1-OPS4: normalize MMC2 accession codes into ENA study target lists.

Ops tooling only -- not part of the Qiita PR.
"""

# Rules mirrored from the Qiita fork at upstream/main c849da1d:
#   platform table:      qiita-control-plane/src/qiita_control_plane/ena_import/platform_mapping.py
#   protocol dispatch:   qiita-control-plane/src/qiita_control_plane/ena_import/protocol_mapping.py:87
#   accession prefixes:  qiita-control-plane/src/qiita_control_plane/ena_import/accession.py
#   large-study / resubmit-skips-stored-runs: workflows/download-ena-study/1.0.0.yaml:65-66,
#                        qiita-control-plane/src/qiita_control_plane/ena_import/registration.py
#   concurrency:         qiita-control-plane/src/qiita_control_plane/dispatch.py:62 (_DISPATCH_CONCURRENCY=8),
#                        qiita-compute-orchestrator/.../jobs/ingest_ena_reads.py:56 (_CONCURRENCY=4)

from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import NamedTuple

PORTAL_BASE = "https://www.ebi.ac.uk/ena/portal/api"

SIZE_FIELDS = [
    "run_accession", "study_accession", "secondary_study_accession",
    "instrument_platform", "library_strategy", "library_source",
    "library_layout", "library_selection", "fastq_ftp", "fastq_bytes",
    "read_count", "scientific_name", "host_tax_id", "target_gene",
    "isolation_source", "tissue_type", "host_body_site", "environmental_medium",
]

# per-job throughput: 0.8 MB/s/conn * _CONCURRENCY=4 (ingest_ena_reads.py:56)
THROUGHPUT_BYTES_PER_S = 0.8e6 * 4
PASS_BYTES = int(THROUGHPUT_BYTES_PER_S * 86400)  # one 24h pass, ~276.48 GB
LARGE_THRESHOLD_BYTES = 200_000_000_000  # 200 GB
# 6 concurrent studies (_DISPATCH_CONCURRENCY=8 minus 2 headroom) * per-job throughput
TRANCHE_BUDGET_BYTES = int(THROUGHPUT_BYTES_PER_S * 6 * 86400)  # ~1.66 TB/day


# ---------------------------------------------------------------------------
# Accession normalization
# ---------------------------------------------------------------------------

STUDY_PREFIXES = ("PRJNA", "PRJEB", "PRJDB", "ERP", "SRP", "DRP")
SAMPLE_PREFIXES = ("SAMN", "SAME", "SAMD", "ERS", "DRS")
RUN_PREFIXES = ("SRR", "ERR", "DRR")
EXPERIMENT_PREFIXES = ("SRX", "ERX", "DRX")
SUBMISSION_PREFIXES = ("SRA", "ERA", "DRA")

QUERY_FIELD_BY_CATEGORY = {
    "sample": "sample_accession",
    "run": "run_accession",
    "experiment": "experiment_accession",
    "submission": "submission_accession",
}

_DROP_EXACT = {"ENA_NOT_FOUND": "ena_not_found"}

# Ordered; first prefix match wins. Chinese-archive (GSA/CNGB) family and
# restricted-access (EGA/dbGaP) family are each collapsed to one reason --
# out of scope per ARCHITECTURE.md D3 and access-restricted respectively.
_DROP_PREFIXES: tuple[tuple[str, str], ...] = (
    ("GSE", "geo"),
    ("subCRA", "cngb_gsa"),
    ("CRA", "cngb_gsa"),
    ("HRA", "cngb_gsa"),
    ("PRJCA", "cngb_gsa"),
    ("CNP", "cngb_gsa"),
    ("OEP", "cngb_gsa"),
    ("CNS", "cngb_gsa"),
    ("CNR", "cngb_gsa"),
    ("CNX", "cngb_gsa"),
    ("CRR", "cngb_gsa"),
    ("SAMC", "cngb_gsa"),
    ("EGAS", "restricted_access"),
    ("EGAD", "restricted_access"),
    ("phs", "restricted_access"),
    ("MTBLS", "metabolomics"),
    ("E-MTAB", "arrayexpress"),
    ("NMDC", "nmdc"),
    ("GCF_", "assembly"),
    ("GCA_", "assembly"),
    ("PXD", "other_non_insdc"),
)

# GenBank/DDBJ/EMBL nucleotide accessions: 1-2 letters + 5-8 digits, optional
# version suffix. Checked only after every named-database rule above, so e.g.
# MTBLS1971 (metabolomics) is never miscaught here.
_GENBANK_NUCLEOTIDE_RE = re.compile(r"^[A-Za-z]{1,2}\d{5,8}(\.\d+)?$")

# A token that doesn't even look like a bare accession (spaces, parens, a
# "(+N more)" summary artifact from the source spreadsheet) is truncated or
# malformed rather than an unrecognized-but-real prefix.
_BARE_ACCESSION_RE = re.compile(r"^[A-Za-z0-9._-]+$")


class TokenClass(NamedTuple):
    category: str  # "study" | "sample" | "run" | "experiment" | "submission" | "drop"
    reason: str | None = None
    query_field: str | None = None


def classify_token(token: str) -> TokenClass:
    """Classify one ';'-split, stripped Accession Code token. Never returns
    an unclassified result -- an unrecognized-but-well-formed token drops
    into "unrecognized_prefix" rather than being silently skipped."""
    if not token:
        raise ValueError("classify_token requires a non-empty token")

    if token in _DROP_EXACT:
        return TokenClass("drop", reason=_DROP_EXACT[token])
    if token.startswith(STUDY_PREFIXES):
        return TokenClass("study")
    if token.startswith(SAMPLE_PREFIXES):
        return TokenClass("sample", query_field="sample_accession")
    if token.startswith(RUN_PREFIXES):
        return TokenClass("run", query_field="run_accession")
    if token.startswith(EXPERIMENT_PREFIXES):
        return TokenClass("experiment", query_field="experiment_accession")
    if token.startswith(SUBMISSION_PREFIXES):
        return TokenClass("submission", query_field="submission_accession")
    for prefix, reason in _DROP_PREFIXES:
        if token.startswith(prefix):
            return TokenClass("drop", reason=reason)
    if _GENBANK_NUCLEOTIDE_RE.match(token):
        return TokenClass("drop", reason="genbank_nucleotide")
    if not _BARE_ACCESSION_RE.match(token):
        return TokenClass("drop", reason="truncated_or_malformed")
    return TokenClass("drop", reason="unrecognized_prefix")


def sum_fastq_bytes(field_value: str) -> int:
    """ENA's fastq_bytes is ';'-separated per mate for paired runs. Sum every
    non-empty part; a non-numeric part fails loud rather than being skipped."""
    total = 0
    for part in field_value.split(";"):
        part = part.strip()
        if not part:
            continue
        total += int(part)
    return total


# ---------------------------------------------------------------------------
# Platform / protocol mapping (mirrors platform_mapping.py / protocol_mapping.py)
# ---------------------------------------------------------------------------

_ENA_PLATFORM_TO_QIITA = {
    "ILLUMINA": "ILLUMINA",
    "PACBIO_SMRT": "PACBIO_SMRT",
    "OXFORD_NANOPORE": "OXFORD_NANOPORE",
    "BGISEQ": "DNBSEQ",
    "DNBSEQ": "DNBSEQ",
    "LS454": "LS454",
    "ION_TORRENT": "ION_TORRENT",
    "COMPLETE_GENOMICS": "COMPLETE_GENOMICS",
}


def map_platform(instrument_platform: str | None) -> str | None:
    """Closed-table platform lookup. Unlike the CP's map_ena_platform this
    never raises -- an unmapped value is counted (unmapped_platform:<v>:N),
    not fatal to the whole batch."""
    if not instrument_platform:
        return None
    return _ENA_PLATFORM_TO_QIITA.get(instrument_platform.strip().upper())


def protocol_category_for_run(library_strategy: str | None, library_source: str | None) -> str:
    """Category dispatch mirrored from protocol_mapping.py:87 (read-length
    bucket omitted -- irrelevant to this report)."""
    strategy = (library_strategy or "").strip().upper()
    source = (library_source or "").strip().upper()
    if strategy == "AMPLICON":
        return "amplicon"
    if strategy == "WGS" or source == "METAGENOMIC":
        return "metagenomics"
    if strategy == "RNA-SEQ" or source in {"TRANSCRIPTOMIC", "METATRANSCRIPTOMIC"}:
        return "transcriptomics"
    return "unmapped"


# ---------------------------------------------------------------------------
# Per-study aggregate + flags
# ---------------------------------------------------------------------------


@dataclass
class StudyMetrics:
    n_runs: int = 0
    n_runs_with_fastq: int = 0
    fastq_bytes: int = 0
    total_read_count: int = 0
    platforms: Counter = field(default_factory=Counter)
    strategies: Counter = field(default_factory=Counter)
    sources: Counter = field(default_factory=Counter)
    layouts: Counter = field(default_factory=Counter)
    platform_pool_bytes: dict = field(default_factory=dict)
    max_pool_bytes: int = 0
    n_platform_pools: int = 0
    unmapped_platform_runs: dict = field(default_factory=dict)  # raw value -> N
    unmapped_strategy_runs: int = 0
    no_fastq_runs: int = 0
    amplicon_runs: int = 0
    metagenomics_runs: int = 0
    transcriptomics_runs: int = 0
    genomic_source_runs: int = 0  # strategy==WGS and source!=METAGENOMIC
    host_tax_ids: Counter = field(default_factory=Counter)
    scientific_names: Counter = field(default_factory=Counter)
    target_genes: set = field(default_factory=set)
    isolation_sources: set = field(default_factory=set)
    tissue_types: set = field(default_factory=set)
    host_body_sites: set = field(default_factory=set)
    environmental_mediums: set = field(default_factory=set)


def derive_study_flags(metrics: StudyMetrics, resolution_status: str) -> tuple[list[str], bool]:
    """Return (flags, is_failure); failure means not_found/ena_error/no_public_runs
    or no_fastq/unmapped_platform covering ALL runs -- everything else is informational."""
    if resolution_status in ("not_found", "ena_error"):
        return [resolution_status], True
    if resolution_status == "no_public_runs" or metrics.n_runs == 0:
        return ["no_public_runs"], True

    flags: list[str] = []
    is_failure = False

    if metrics.no_fastq_runs:
        flags.append(f"no_fastq:{metrics.no_fastq_runs}")
        if metrics.no_fastq_runs == metrics.n_runs:
            is_failure = True

    for raw_value, n in sorted(metrics.unmapped_platform_runs.items()):
        flags.append(f"unmapped_platform:{raw_value}:{n}")
    total_unmapped_platform = sum(metrics.unmapped_platform_runs.values())
    if total_unmapped_platform:
        if total_unmapped_platform == metrics.n_runs:
            is_failure = True
        else:
            flags.append("unmapped_platform partial")

    if metrics.unmapped_strategy_runs:
        flags.append(f"unmapped_strategy:{metrics.unmapped_strategy_runs}")

    if metrics.amplicon_runs:
        flags.append(f"amplicon:{metrics.amplicon_runs}")
    if metrics.amplicon_runs and (metrics.metagenomics_runs or metrics.transcriptomics_runs):
        flags.append("mixed_assay")
    if metrics.genomic_source_runs:
        flags.append(f"genomic_source:{metrics.genomic_source_runs}")

    return flags, is_failure


def study_coarse_assay(metrics: StudyMetrics) -> str:
    buckets = {
        "AMPLICON": metrics.amplicon_runs,
        "WGS": metrics.metagenomics_runs,
        "RNA-SEQ": metrics.transcriptomics_runs,
    }
    nonzero = [name for name, n in buckets.items() if n]
    if len(nonzero) == 1:
        return nonzero[0]
    return "MIXED"


def is_single_assay(metrics: StudyMetrics) -> bool:
    if metrics.n_runs == 0:
        return False
    return metrics.amplicon_runs == metrics.n_runs or metrics.metagenomics_runs == metrics.n_runs


def est_job_hours(max_pool_bytes: int) -> float:
    if max_pool_bytes <= 0:
        return 0.0
    return max_pool_bytes / THROUGHPUT_BYTES_PER_S / 3600


def est_passes(max_pool_bytes: int) -> int:
    if max_pool_bytes <= 0:
        return 0
    return math.ceil(max_pool_bytes / PASS_BYTES)


# ---------------------------------------------------------------------------
# Sample-hint / manual-review heuristic (target_candidates.tsv)
# ---------------------------------------------------------------------------

_FECAL_TERMS = ("fecal", "faecal", "feces", "stool")
_TISSUE_TERMS = ("tissue", "biopsy", "mucosa", "tumor", "tumour")
_16S_TERMS = ("16s",)


def _hint_values(hints: dict) -> list[str]:
    values = []
    for v in hints.values():
        values.extend(x.lower() for x in v)
    return values


def sample_hint(hints: dict) -> str:
    values = " ".join(_hint_values(hints))
    is_fecal = any(t in values for t in _FECAL_TERMS)
    is_tissue = any(t in values for t in _TISSUE_TERMS)
    if is_fecal and is_tissue:
        return "mixed"
    if is_fecal:
        return "fecal"
    if is_tissue:
        return "tissue"
    return "unknown"


def needs_manual_review(coarse_assay: str, hints: dict, target_genes: set) -> bool:
    if sample_hint(hints) == "unknown":
        return True
    if coarse_assay == "AMPLICON":
        genes = " ".join(g.lower() for g in target_genes)
        if not any(t in genes for t in _16S_TERMS):
            return True
    return False


# ---------------------------------------------------------------------------
# Tranche assignment (provisional; final pass runs after user curation)
# ---------------------------------------------------------------------------


@dataclass
class TrancheResult:
    tranches: list
    hold: list


def assign_tranches(
    studies: list[tuple[str, int]],
    budget_bytes: int = TRANCHE_BUDGET_BYTES,
    large_threshold: int = LARGE_THRESHOLD_BYTES,
    pass_bytes: int = PASS_BYTES,
) -> TrancheResult:
    """Largest-first-into-lightest-tranche bin balance for small studies;
    large studies (> large_threshold) get their own consecutive tranche
    block, one pass per tranche, at most one large study per tranche.
    Zero/negative-byte studies (no route) go to `hold`, never a tranche."""
    large = [(acc, b) for acc, b in studies if b > large_threshold]
    small = [(acc, b) for acc, b in studies if 0 < b <= large_threshold]
    hold = [acc for acc, b in studies if b <= 0]

    large_passes = {acc: math.ceil(b / pass_bytes) for acc, b in large}
    total_large_tranches = sum(large_passes.values())
    k_small = math.ceil(sum(b for _, b in small) / budget_bytes) if small else 0
    total_tranches = max(total_large_tranches, k_small, 1 if (large or small) else 0)

    tranches: list[list[str]] = [[] for _ in range(total_tranches)]

    idx = 0
    for acc, _ in large:
        for _ in range(large_passes[acc]):
            tranches[idx].append(acc)
            idx += 1

    heap = [(0, i) for i in range(total_tranches)]
    heapq.heapify(heap)
    for acc, b in sorted(small, key=lambda x: -x[1]):
        load, i = heapq.heappop(heap)
        tranches[i].append(acc)
        heapq.heappush(heap, (load + b, i))

    return TrancheResult(tranches=tranches, hold=hold)


_T1_BUDGET_BYTES = 300_000_000_000
_T1_SIZE_BANDS = (
    ("under_5gb", 5_000_000_000),
    ("5_to_50gb", 50_000_000_000),
    ("50_to_250gb", 250_000_000_000),
)


def size_band(fastq_bytes: int) -> str:
    for name, upper in _T1_SIZE_BANDS:
        if fastq_bytes < upper:
            return name
    return _T1_SIZE_BANDS[-1][0]


def suggest_tranche_1(candidates: list[dict], exclude: frozenset = frozenset()) -> list[dict]:
    """Heuristic draft pick for the comment.md T1 suggestion: at least one
    study per size band (<5GB / 5-50GB / 50-250GB) where candidates exist,
    >=1 AMPLICON + >=1 WGS, a non-Illumina platform if any exists in the
    pool, <=300GB total, 3-5 studies, review-clean candidates preferred.
    Final pick is the user's."""
    pool = [c for c in candidates if c["study"] not in exclude]
    clean_by_band: dict[str, list[dict]] = {b: [] for b, _ in _T1_SIZE_BANDS}
    review_by_band: dict[str, list[dict]] = {b: [] for b, _ in _T1_SIZE_BANDS}
    for c in pool:
        bucket = review_by_band if c["needs_manual_review"] else clean_by_band
        bucket[size_band(c["fastq_bytes"])].append(c)

    picked: list[dict] = []
    total = 0

    def _try_add(c: dict) -> bool:
        nonlocal total
        if c in picked or total + c["fastq_bytes"] > _T1_BUDGET_BYTES:
            return False
        picked.append(c)
        total += c["fastq_bytes"]
        return True

    def _pick(band: str, *, prefer=lambda c: True) -> bool:
        for bucket in (clean_by_band[band], review_by_band[band]):
            for c in sorted(bucket, key=lambda c: c["fastq_bytes"]):
                if prefer(c) and _try_add(c):
                    return True
        return False

    def _ensure(predicate) -> None:
        if any(predicate(c) for c in picked):
            return
        for band, _ in _T1_SIZE_BANDS:
            if _pick(band, prefer=predicate):
                return

    for band, _ in _T1_SIZE_BANDS:
        _pick(band)

    _ensure(lambda c: c["coarse_assay"] == "AMPLICON")
    _ensure(lambda c: c["coarse_assay"] == "WGS")
    if any(c["platform"] and c["platform"] != "ILLUMINA" for c in pool):
        _ensure(lambda c: c["platform"] and c["platform"] != "ILLUMINA")

    while len(picked) < 3 and any(_pick(band) for band, _ in _T1_SIZE_BANDS):
        pass

    return sorted(picked, key=lambda c: c["fastq_bytes"])[:5]


# ---------------------------------------------------------------------------
# ENA HTTP client: cached, rate-limited (<=3 req/s), 3 retries on 5xx/timeout
# ---------------------------------------------------------------------------


class EnaRequestError(RuntimeError):
    def __init__(self, url: str, cause: Exception | None):
        self.url = url
        self.cause = cause
        super().__init__(f"ENA request failed for {url!r}: {cause!r}")


class EnaClient:
    def __init__(self, cache_dir: Path, min_interval: float = 1.0 / 3.0, max_retries: int = 3):
        cache_dir.mkdir(parents=True, exist_ok=True)
        self._cache_dir = cache_dir
        self._min_interval = min_interval
        self._max_retries = max_retries
        self._last_request_t = 0.0

    def _cache_path(self, url: str) -> Path:
        return self._cache_dir / (hashlib.sha256(url.encode()).hexdigest() + ".tsv")

    def _throttle(self) -> None:
        now = time.monotonic()
        wait = self._min_interval - (now - self._last_request_t)
        if wait > 0:
            time.sleep(wait)
        self._last_request_t = time.monotonic()

    def get(self, url: str) -> str:
        cache_path = self._cache_path(url)
        if cache_path.exists():
            return cache_path.read_text()

        last_exc: Exception | None = None
        for attempt in range(self._max_retries + 1):
            self._throttle()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "miint-ops-m1-ops4/1.0"})
                with urllib.request.urlopen(req, timeout=30) as resp:
                    body = resp.read().decode("utf-8")
                cache_path.write_text(body)
                return body
            except urllib.error.HTTPError as exc:
                if exc.code < 500:
                    raise EnaRequestError(url, exc) from exc
                last_exc = exc
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                last_exc = exc
            if attempt < self._max_retries:
                time.sleep(2 ** attempt)
        raise EnaRequestError(url, last_exc)


def parse_tsv(body: str) -> list[dict]:
    lines = [ln for ln in body.splitlines() if ln != ""]
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def build_study_search_url(accession: str) -> str:
    query = f'study_accession="{accession}" OR secondary_study_accession="{accession}"'
    params = {
        "result": "study",
        "query": query,
        "fields": "study_accession,secondary_study_accession,study_title",
        "limit": "0",
        "format": "tsv",
    }
    return f"{PORTAL_BASE}/search?{urllib.parse.urlencode(params)}"


def build_map_search_url(query_field: str, accession: str) -> str:
    params = {
        "result": "read_run",
        "query": f'{query_field}="{accession}"',
        "fields": "study_accession,secondary_study_accession",
        "limit": "0",
        "format": "tsv",
    }
    return f"{PORTAL_BASE}/search?{urllib.parse.urlencode(params)}"


def build_filereport_url(accession: str, fields: list[str] = SIZE_FIELDS) -> str:
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(fields),
        "format": "tsv",
    }
    return f"{PORTAL_BASE}/filereport?{urllib.parse.urlencode(params)}"


def resolve_study_token(client: EnaClient, accession: str) -> tuple[str | None, str]:
    try:
        body = client.get(build_study_search_url(accession))
    except EnaRequestError:
        return None, "ena_error"
    rows = parse_tsv(body)
    if not rows:
        return None, "not_found"
    return rows[0]["study_accession"], "resolved"


def resolve_map_token(client: EnaClient, query_field: str, accession: str) -> tuple[str | None, str]:
    try:
        body = client.get(build_map_search_url(query_field, accession))
    except EnaRequestError:
        return None, "ena_error"
    rows = parse_tsv(body)
    if not rows:
        return None, "not_found"
    return rows[0]["study_accession"], "resolved"


def size_study(client: EnaClient, accession: str) -> tuple[StudyMetrics, str]:
    try:
        body = client.get(build_filereport_url(accession))
    except EnaRequestError:
        return StudyMetrics(), "ena_error"
    rows = parse_tsv(body)
    if not rows:
        return StudyMetrics(), "no_public_runs"

    m = StudyMetrics()
    for r in rows:
        m.n_runs += 1
        platform_raw = r.get("instrument_platform", "").strip()
        strategy = r.get("library_strategy", "").strip()
        source = r.get("library_source", "").strip()
        m.platforms[platform_raw or "(blank)"] += 1
        m.strategies[strategy or "(blank)"] += 1
        m.sources[source or "(blank)"] += 1
        m.layouts[r.get("library_layout", "").strip() or "(blank)"] += 1

        fastq_ftp = r.get("fastq_ftp", "").strip()
        fastq_bytes_field = r.get("fastq_bytes", "").strip()
        run_bytes = sum_fastq_bytes(fastq_bytes_field) if fastq_bytes_field else 0
        has_fastq = bool(fastq_ftp) and run_bytes > 0
        if has_fastq:
            m.n_runs_with_fastq += 1
            m.fastq_bytes += run_bytes
        else:
            m.no_fastq_runs += 1

        read_count = r.get("read_count", "").strip()
        if read_count.isdigit():
            m.total_read_count += int(read_count)

        mapped = map_platform(platform_raw)
        if mapped is None:
            key = platform_raw or "(blank)"
            m.unmapped_platform_runs[key] = m.unmapped_platform_runs.get(key, 0) + 1
        elif has_fastq:
            pool_key = mapped
            m.platform_pool_bytes[pool_key] = m.platform_pool_bytes.get(pool_key, 0) + run_bytes

        category = protocol_category_for_run(strategy, source)
        if category == "amplicon":
            m.amplicon_runs += 1
        elif category == "metagenomics":
            m.metagenomics_runs += 1
        elif category == "transcriptomics":
            m.transcriptomics_runs += 1
        else:
            m.unmapped_strategy_runs += 1
        if strategy.upper() == "WGS" and source.upper() != "METAGENOMIC":
            m.genomic_source_runs += 1

        sci_name = r.get("scientific_name", "").strip()
        if sci_name:
            m.scientific_names[sci_name] += 1
        host_tax_id = r.get("host_tax_id", "").strip()
        if host_tax_id:
            m.host_tax_ids[host_tax_id] += 1
        for field_name, target in (
            ("target_gene", m.target_genes),
            ("isolation_source", m.isolation_sources),
            ("tissue_type", m.tissue_types),
            ("host_body_site", m.host_body_sites),
            ("environmental_medium", m.environmental_mediums),
        ):
            # some depositors pack several terms into one field (';'-joined,
            # same convention as paired fastq_bytes) -- split before adding.
            for part in r.get(field_name, "").split(";"):
                part = part.strip()
                if part:
                    target.add(part)

    m.max_pool_bytes = max(m.platform_pool_bytes.values()) if m.platform_pool_bytes else 0
    m.n_platform_pools = len(m.platform_pool_bytes)
    return m, "resolved"


# ---------------------------------------------------------------------------
# CSV -> token pool
# ---------------------------------------------------------------------------


@dataclass
class Citation:
    doi: str
    pmid: str
    title: str


def load_rows(csv_path: Path) -> list[dict]:
    with csv_path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def tokenize(rows: list[dict]) -> list[tuple[str, dict]]:
    """One (token, source_row) pair per ';'-split Accession Code entry."""
    out = []
    for r in rows:
        raw = (r.get("Accession Code") or "").strip()
        if not raw:
            continue
        for tok in raw.split(";"):
            tok = tok.strip()
            if tok:
                out.append((tok, r))
    return out


def normalize(client: EnaClient, rows: list[dict]) -> tuple[Counter, dict[str, dict]]:
    """Classify + resolve every token. Returns (bucket_counts, pool) where
    pool maps resolved PRJ accession -> {"tokens": [...], "citations": [Citation, ...]}."""
    token_rows = tokenize(rows)
    counts: Counter = Counter()
    pool: dict[str, dict] = {}

    resolved_cache: dict[tuple[str, str], tuple[str | None, str]] = {}

    for token, row in token_rows:
        counts["total_tokens"] += 1
        cls = classify_token(token)
        if cls.category == "drop":
            counts[f"drop:{cls.reason}"] += 1
            continue

        if cls.category == "study":
            cache_key = ("study", token)
            resolve_fn = lambda: resolve_study_token(client, token)
        else:
            cache_key = (cls.category, token)
            resolve_fn = lambda cat=cls.category, tok=token: resolve_map_token(
                client, QUERY_FIELD_BY_CATEGORY[cat], tok
            )

        if cache_key not in resolved_cache:
            resolved_cache[cache_key] = resolve_fn()
        study_acc, status = resolved_cache[cache_key]

        counts[f"lookup:{cls.category}:{status}"] += 1
        if status != "resolved" or not study_acc:
            counts["unresolved_after_lookup"] += 1
            continue

        entry = pool.setdefault(study_acc, {"tokens": [], "citations": []})
        entry["tokens"].append(token)
        entry["citations"].append(
            Citation(
                doi=(row.get("doi") or "").strip(),
                pmid=(row.get("pubmed_id") or "").strip(),
                title=(row.get("title") or "").strip(),
            )
        )

    counts["resolved_distinct_studies"] = len(pool)
    return counts, pool


# ---------------------------------------------------------------------------
# TSV writers
# ---------------------------------------------------------------------------


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def write_normalization_report(path: Path, counts: Counter) -> None:
    header = ["bucket", "count"]
    rows = [[k, v] for k, v in sorted(counts.items())]
    write_tsv(path, header, rows)


_FIELD_DISPLAY_LIMIT = 8  # cap wide multi-value fields (e.g. EMP500's ~800 host_tax_ids)


def _counter_str(c: Counter, limit: int = _FIELD_DISPLAY_LIMIT) -> str:
    items = c.most_common()
    shown = sorted(items[:limit])
    text = ";".join(f"{k}:{v}" for k, v in shown)
    if len(items) > limit:
        text += f";(+{len(items) - limit} more distinct)"
    return text


def _dict_str(d: dict, limit: int = _FIELD_DISPLAY_LIMIT) -> str:
    items = sorted(d.items(), key=lambda kv: -kv[1])
    shown = sorted(items[:limit])
    text = ";".join(f"{k}:{v}" for k, v in shown)
    if len(items) > limit:
        text += f";(+{len(items) - limit} more distinct)"
    return text


def _set_str(s: set, limit: int = _FIELD_DISPLAY_LIMIT) -> str:
    items = sorted(s)
    text = ";".join(items[:limit])
    if len(items) > limit:
        text += f";(+{len(items) - limit} more distinct)"
    return text


def build_study_row(accession: str, entry: dict, metrics: StudyMetrics, status: str) -> dict:
    flags, is_failure = derive_study_flags(metrics, status)
    citations = entry.get("citations", []) if entry else []
    dois = sorted({c.doi for c in citations if c.doi})
    pmids = sorted({c.pmid for c in citations if c.pmid})
    titles = sorted({c.title for c in citations if c.title})
    hints = {
        "isolation_source": metrics.isolation_sources,
        "tissue_type": metrics.tissue_types,
        "host_body_site": metrics.host_body_sites,
        "environmental_medium": metrics.environmental_mediums,
    }
    return {
        "study": accession,
        "n_citing_tokens": len(entry.get("tokens", [])) if entry else 0,
        "citing_dois": ";".join(dois),
        "citing_pmids": ";".join(pmids),
        "citing_titles": ";".join(titles),
        "resolution_status": status,
        "n_runs": metrics.n_runs,
        "n_runs_with_fastq": metrics.n_runs_with_fastq,
        "fastq_bytes": metrics.fastq_bytes,
        "total_read_count": metrics.total_read_count,
        "platforms": _counter_str(metrics.platforms),
        "strategies": _counter_str(metrics.strategies),
        "sources": _counter_str(metrics.sources),
        "layouts": _counter_str(metrics.layouts),
        "platform_pool_bytes": _dict_str(metrics.platform_pool_bytes),
        "max_pool_bytes": metrics.max_pool_bytes,
        "n_platform_pools": metrics.n_platform_pools,
        "flags": ";".join(flags),
        "is_failure": is_failure,
        "est_job_hours": round(est_job_hours(metrics.max_pool_bytes), 2),
        "est_passes": est_passes(metrics.max_pool_bytes),
        "coarse_assay": study_coarse_assay(metrics),
        "target_gene": _set_str(metrics.target_genes),
        "isolation_source": _set_str(metrics.isolation_sources),
        "tissue_type": _set_str(metrics.tissue_types),
        "host_body_site": _set_str(metrics.host_body_sites),
        "environmental_medium": _set_str(metrics.environmental_mediums),
        "host_tax_id": _counter_str(metrics.host_tax_ids),
        "scientific_name": _counter_str(metrics.scientific_names),
        "hints": hints,
    }


TARGET_POOL_COLUMNS = [
    "study", "n_citing_tokens", "citing_dois", "citing_pmids", "citing_titles",
    "resolution_status", "n_runs", "n_runs_with_fastq", "fastq_bytes",
    "total_read_count", "platforms", "strategies", "sources", "layouts",
    "platform_pool_bytes", "max_pool_bytes", "n_platform_pools", "flags",
    "is_failure", "est_job_hours", "est_passes", "coarse_assay", "target_gene",
    "isolation_source", "tissue_type", "host_body_site", "environmental_medium",
    "host_tax_id", "scientific_name",
]


def write_target_pool(path: Path, study_rows: list[dict]) -> None:
    rows = [[r[c] for c in TARGET_POOL_COLUMNS] for r in study_rows]
    write_tsv(path, TARGET_POOL_COLUMNS, rows)


CANDIDATE_COLUMNS = [
    "study", "paper_titles", "paper_dois", "fastq_bytes", "est_job_hours",
    "platforms", "strategies", "sources", "layouts", "target_gene",
    "isolation_source", "tissue_type", "host_body_site", "environmental_medium",
    "host_tax_id", "coarse_assay", "sample_hint", "needs_manual_review",
]


def build_candidates(study_rows: list[dict]) -> list[dict]:
    candidates = []
    for r in study_rows:
        if r["is_failure"]:
            continue
        # single assay: all AMPLICON, or all WGS/METAGENOMIC
        n_runs = r["n_runs"]
        if n_runs == 0:
            continue
        if r["coarse_assay"] not in ("AMPLICON", "WGS"):
            continue
        if "mixed_assay" in r["flags"]:
            continue
        if r["est_job_hours"] > 24:
            continue
        host_tax_id = r["host_tax_id"]
        sci_name = r["scientific_name"].lower()
        is_human_tax_id = "9606:" in host_tax_id
        is_human_metagenome_name = "human" in sci_name and "metagenome" in sci_name
        if not (is_human_tax_id or is_human_metagenome_name):
            continue

        hints = r["hints"]
        target_genes = set(r["target_gene"].split(";")) if r["target_gene"] else set()
        hint = sample_hint(hints)
        review = needs_manual_review(r["coarse_assay"], hints, target_genes)
        first_pool = r["platform_pool_bytes"].split(";")[0] if r["platform_pool_bytes"] else ""
        platform = first_pool.split(":")[0]
        candidates.append({
            "study": r["study"],
            "paper_titles": r["citing_titles"],
            "paper_dois": r["citing_dois"],
            "fastq_bytes": r["fastq_bytes"],
            "est_job_hours": r["est_job_hours"],
            "platforms": r["platforms"],
            "strategies": r["strategies"],
            "sources": r["sources"],
            "layouts": r["layouts"],
            "target_gene": r["target_gene"],
            "isolation_source": r["isolation_source"],
            "tissue_type": r["tissue_type"],
            "host_body_site": r["host_body_site"],
            "environmental_medium": r["environmental_medium"],
            "host_tax_id": r["host_tax_id"],
            "coarse_assay": r["coarse_assay"],
            "sample_hint": hint,
            "needs_manual_review": review,
            "platform": platform,
        })
    candidates.sort(key=lambda c: (c["coarse_assay"], c["sample_hint"], c["fastq_bytes"]))
    return candidates


def write_candidates(path: Path, candidates: list[dict]) -> None:
    rows = [[c[col] for col in CANDIDATE_COLUMNS] for c in candidates]
    write_tsv(path, CANDIDATE_COLUMNS, rows)


# ---------------------------------------------------------------------------
# references.tsv
# ---------------------------------------------------------------------------

REFERENCE_STUDIES = [
    # (panel, accession_or_None, note, tranche_override)
    ("emp500_shotgun", "PRJEB42019", "local_copy_on_barnacle", "local"),
    ("mouse_gut", "PRJEB7759", "not_held_locally; full download; pending T17-3; large lane", "large"),
    ("hmp_body_site", "PRJNA48479", "", None),
    ("hmp_body_site", "PRJNA43021", "", None),
    ("hmp_body_site", "PRJNA48333", "HMP 16S 454", None),
    ("hmp_body_site", None, "unresolved (shotgun)", None),
    ("emp_16s", None, "unresolved", None),
]


def reference_tranche(
    accession: str | None, tranche_override: str | None, is_failure: bool, fastq_bytes: int
) -> str:
    """Same hold/local/large vocabulary as assign_tranches' TrancheResult.
    Explicit overrides win; otherwise hold (no-FASTQ) or large (>200GB) is
    derived from the sizing outcome, and an unresolved row is blank."""
    if accession is None:
        return ""
    if tranche_override:
        return tranche_override
    if is_failure:
        return "hold"
    if fastq_bytes > LARGE_THRESHOLD_BYTES:
        return "large"
    return ""


def build_references(client: EnaClient) -> list[dict]:
    today = date.today().isoformat()
    out = []
    for panel, accession, note, tranche_override in REFERENCE_STUDIES:
        if accession is None:
            out.append({
                "panel": panel, "study": "", "tranche": "", "paper_titles": "", "paper_dois": "",
                "fastq_bytes": "", "est_job_hours": "", "platforms": "", "strategies": "",
                "sources": "", "layouts": "", "target_gene": "", "isolation_source": "",
                "tissue_type": "", "host_body_site": "", "environmental_medium": "",
                "host_tax_id": "", "coarse_assay": "", "sample_hint": "", "needs_manual_review": "",
                "verified": "no", "note": note,
            })
            continue
        metrics, status = size_study(client, accession)
        flags, is_failure = derive_study_flags(metrics, status)
        note_full = ";".join(x for x in (note, ";".join(flags)) if x)
        out.append({
            "panel": panel, "study": accession,
            "tranche": reference_tranche(accession, tranche_override, is_failure, metrics.fastq_bytes),
            "paper_titles": "", "paper_dois": "",
            "fastq_bytes": metrics.fastq_bytes, "est_job_hours": round(est_job_hours(metrics.max_pool_bytes), 2),
            "platforms": _counter_str(metrics.platforms), "strategies": _counter_str(metrics.strategies),
            "sources": _counter_str(metrics.sources), "layouts": _counter_str(metrics.layouts),
            "target_gene": _set_str(metrics.target_genes),
            "isolation_source": _set_str(metrics.isolation_sources),
            "tissue_type": _set_str(metrics.tissue_types),
            "host_body_site": _set_str(metrics.host_body_sites),
            "environmental_medium": _set_str(metrics.environmental_mediums),
            "host_tax_id": _counter_str(metrics.host_tax_ids),
            "coarse_assay": study_coarse_assay(metrics),
            "sample_hint": "", "needs_manual_review": "",
            "verified": today, "note": note_full,
        })
    return out


REFERENCE_COLUMNS = ["panel", "study", "tranche", "paper_titles", "paper_dois", "fastq_bytes",
                      "est_job_hours", "platforms", "strategies", "sources", "layouts", "target_gene",
                      "isolation_source", "tissue_type", "host_body_site", "environmental_medium",
                      "host_tax_id", "coarse_assay", "sample_hint", "needs_manual_review",
                      "verified", "note"]


def write_references(path: Path, refs: list[dict]) -> None:
    rows = [[r[c] for c in REFERENCE_COLUMNS] for r in refs]
    write_tsv(path, REFERENCE_COLUMNS, rows)


# ---------------------------------------------------------------------------
# comment.md
# ---------------------------------------------------------------------------


def write_comment(path: Path, counts: Counter, candidates: list[dict], refs_tsv_text: str,
                   suggested_t1: list[dict]) -> None:
    n_16s_fecal = sum(1 for c in candidates if c["coarse_assay"] == "AMPLICON" and c["sample_hint"] == "fecal")
    n_16s_tissue = sum(1 for c in candidates if c["coarse_assay"] == "AMPLICON" and c["sample_hint"] == "tissue")
    n_wgs_fecal = sum(1 for c in candidates if c["coarse_assay"] == "WGS" and c["sample_hint"] == "fecal")
    n_wgs_tissue = sum(1 for c in candidates if c["coarse_assay"] == "WGS" and c["sample_hint"] == "tissue")

    t1_lines = "\n".join(f"- `{c['study']}` ({c['coarse_assay']}, {c['sample_hint']}, "
                          f"{c['fastq_bytes']:,} bytes)" for c in suggested_t1) or "- (none met the T1 heuristic)"
    counts_lines = "\n".join(f"{k}\t{v}" for k, v in sorted(counts.items()))

    body = f"""# M1-OPS4: ENA target normalization (provisional, unverified)

**Status:** draft. Throughput/tranche arithmetic and the T1 suggestion below are
provisional until the owner curates `target_candidates.tsv` into a real
`targets.tsv`.

## Arithmetic (provisional)

- Per-job throughput: 0.8 MB/s/conn * `_CONCURRENCY=4`
  (`qiita-compute-orchestrator/.../jobs/ingest_ena_reads.py:56`) = 3.2 MB/s ~=
  {PASS_BYTES:,} bytes/24h pass (~276 GB).
- Large-study threshold: > 200 GB. `est_passes = ceil(max_pool_bytes / {PASS_BYTES:,})`.
- Tranche budget: 6 concurrent studies (`_DISPATCH_CONCURRENCY=8` at
  `qiita-control-plane/src/qiita_control_plane/dispatch.py:62`, minus 2
  headroom) * 3.2 MB/s = 19.2 MB/s ~= {TRANCHE_BUDGET_BYTES:,} bytes/day (~1.66 TB).

## Normalization counts

```
{counts_lines}
```

Per-prefix drop/lookup counts sum to `total_tokens` (see `normalization_report.tsv`).

## Candidate pool (target_candidates.tsv)

| assay | sample_hint | n |
|---|---|---|
| AMPLICON | fecal | {n_16s_fecal} |
| AMPLICON | tissue | {n_16s_tissue} |
| WGS | fecal | {n_wgs_fecal} |
| WGS | tissue | {n_wgs_tissue} |

Full breakdown, including `mixed`/`unknown` sample_hint and
`needs_manual_review`, is in `target_candidates.tsv`.

## targets.tsv (placeholder)

<!-- TODO(owner): after curating ~2 studies per {{16S, shotgun}} x {{fecal,
tissue}} from target_candidates.tsv, this section gets replaced with the
final targets.tsv contents and the curated study list. -->

## Large-study protocol

Studies with `fastq_bytes > 200 GB` go in a "large" lane: `est_passes =
ceil(max_pool_bytes / ~276 GB)`, one 24h pass per successive tranche, at
most one large study per tranche. Retries/re-imports skip runs already
stored (`workflows/download-ena-study/1.0.0.yaml:65-66`,
`ena_import/registration.py`); a study is complete when stored read counts
match ENA's `read_count`. Never two concurrent imports of the same study
(qiita#578).

## References (references.tsv)

<details><summary>references.tsv</summary>

```tsv
{refs_tsv_text}
```

</details>

## Suggested Tranche 1 (draft; final pick after curation)

{t1_lines}

**Tranche 1: TBD after curation.**
"""
    path.write_text(body)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()

    out_dir = args.out_dir
    cache_dir = out_dir / "cache"
    client = EnaClient(cache_dir)

    rows = load_rows(args.csv)
    counts, pool = normalize(client, rows)

    study_rows = []
    for accession, entry in sorted(pool.items()):
        metrics, status = size_study(client, accession)
        study_rows.append(build_study_row(accession, entry, metrics, status))

    write_normalization_report(out_dir / "normalization_report.tsv", counts)
    write_target_pool(out_dir / "target_pool.tsv", study_rows)

    candidates = build_candidates(study_rows)
    write_candidates(out_dir / "target_candidates.tsv", candidates)

    refs = build_references(client)
    write_references(out_dir / "references.tsv", refs)
    refs_tsv_text = (out_dir / "references.tsv").read_text().rstrip("\n")

    suggested_t1 = suggest_tranche_1(candidates, exclude=frozenset({"PRJDB13464", "PRJDB35947"}))
    write_comment(out_dir / "comment.md", counts, candidates, refs_tsv_text, suggested_t1)

    print(f"wrote {len(study_rows)} target_pool rows, {len(candidates)} candidates, "
          f"{len(refs)} reference rows to {out_dir}")


if __name__ == "__main__":
    main()
