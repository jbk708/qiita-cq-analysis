"""Per-study 16S check: region, primers, trim length and multiplexing, from each run's Qiita parquet.

Feeds the `amplicon` submit args (#244): `primer` + `orient_primer` (strip the forward primer when
present) and `trim` (R1 truncation length; #244 denoises R1 only). GG2 V4 exact-match needs reads
anchored at 515F and trim 150 (90/100 also exist in GG2 but cover less of the catalog).
Usage: amplicon_check.py PRJ... [--manifest M] [--reads N] [--out study.tsv]
"""
import argparse
import os
import sys

import duckdb
import pandas as pd

MIINT_REPO = "https://ftp.microbio.me/pub/miint"
MANIFEST = os.path.join(os.environ.get("QDEV_SHARE", f"/ddn_scratch/{os.environ.get('USER')}/qiita-pilot/qiita-dev-share"), "manifest.tsv")

# region -> (forward primer, reverse primer). Several share 515F/341F; the reverse primer on R2
# tells them apart, and single-end runs report every region their forward primer fits.
REGIONS = {
    "V1-V2": ("AGAGTTTGATYMTGGCTCAG", "TGCTGCCTCCCGTAGGAGT"),
    "V3": ("CCTACGGGAGGCAGCAG", "ATTACCGCGGCTGCTGG"),
    "V3-V4": ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"),
    "V4": ("GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"),
    "V4-V5": ("GTGYCAGCMGCCGCGGTAA", "CCGYCAATTYMTTTRAGTTT"),
    "V5-V6": ("AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"),
}
# read starts just past each forward primer, for runs whose primers were already removed.
ANCHORS = {"V1-V2": "^GA[CT]GAACGC", "V3/V3-V4": "^T[AG]GGGAAT", "V4": "^TAC[GA]"}
HIT = 0.70
SLACK = 30  # bases of barcode/linker allowed before a primer
GG2_TRIMS = (150, 100, 90)

RUN_SQL = """
WITH r AS (
  SELECT upper(sequence1) AS s1, upper(sequence2) AS s2, qual1, qual2
  FROM read_parquet($path) LIMIT $n
), f AS (
  SELECT *,
    regexp_extract(left(s1, $flen + {slack}), $fwd_core) AS fhit,
    regexp_matches(left(s1, $flen + {slack}), $fwd) AS ffull,
    regexp_matches(left(s2, $rlen + {slack}), $rev) AS rhit,
    regexp_matches(left(s1, $rlen + {slack}), $rev) AS r_on_r1
  FROM r
), g AS (
  SELECT *, CASE WHEN fhit <> '' THEN strpos(s1, fhit) - 1 END AS off,
         CASE WHEN fhit <> '' THEN length(s1) - (strpos(s1, fhit) - 1 + $flen) END AS post
  FROM f
)
SELECT count(*) AS n,
       count(s2) > 0 AS paired,
       avg((fhit <> '')::INT) AS fwd_rate,
       avg(ffull::INT) AS fwd_full_rate,
       avg(rhit::INT) AS rev_rate,
       avg(r_on_r1::INT) AS rev_on_r1_rate,
       median(off) AS off_median,
       quantile_disc(post, 0.10) AS post_p10,
       quantile_disc(length(s1), 0.10) AS len_p10,
       {anchors}
       count(DISTINCT left(s1, off)) FILTER (WHERE off >= 6) AS prefixes
FROM g
""".format(slack=SLACK, anchors="".join(
    f"avg(regexp_matches(s1, '{rx}')::INT) AS \"anchor:{r}\",\n" for r, rx in ANCHORS.items()))

MERGE_SQL = """
SELECT avg((merge_pairs_vsearch(sequence1, qual1, sequence2, qual2)).merged::INT)
FROM (SELECT * FROM read_parquet($path) WHERE sequence2 IS NOT NULL LIMIT 2000)
"""


def check_run(con, path, n):
    rx = lambda p: con.execute("SELECT sequence_dna_as_regexp($p)", {"p": p}).fetchone()[0]
    out = {}
    for region, (fwd, rev) in REGIONS.items():
        # the primer's last 2 bases may not match off-target templates; the core still anchors it.
        row = con.execute(RUN_SQL, {"path": path, "n": n, "fwd": rx(fwd), "fwd_core": rx(fwd[:-2]),
                                    "rev": rx(rev),
                                    "flen": len(fwd), "rlen": len(rev)}).fetchone()
        cols = [d[0] for d in con.description]
        out[region] = dict(zip(cols, row))
    any_r = next(iter(out.values()))
    merge = con.execute(MERGE_SQL, {"path": path}).fetchone()[0] if any_r["paired"] else None
    return out, merge


def verdict(per_region, merge):
    paired = next(iter(per_region.values()))["paired"]
    fits = [r for r, m in per_region.items()
            if m["fwd_rate"] >= HIT and (not paired or m["rev_rate"] >= HIT)]
    flags = []
    if fits:
        primer_present = True
        best = per_region[fits[0]]
        trim90 = best["post_p10"]
        fwd = REGIONS[fits[0]][0]
        if (best["off_median"] or 0) > 0:
            flags.append("prefix-before-primer")
        if best["prefixes"] >= 8:
            flags.append("multiplexed?")
        if best["rev_on_r1_rate"] >= 0.2:
            flags.append("mixed-orientation")
        if best["fwd_full_rate"] < best["fwd_rate"] - 0.2:
            flags.append(f"primer-3'-mismatch:{1 - best['fwd_full_rate']:.0%}")
    else:
        primer_present = False
        fwd = None
        m = per_region["V4"]
        trim90 = m["len_p10"]
        fits = [r + "*" for r in ANCHORS if m[f"anchor:{r}"] >= HIT]  # * = primers already removed
        if not fits:
            flags.append("no-region")
    if len(fits) > 1:
        flags.append("ambiguous:" + "/".join(fits))
    anchored_515f = fwd == REGIONS["V4"][0] or "V4*" in fits
    gg2 = next((t for t in GG2_TRIMS if anchored_515f and trim90 and trim90 >= t), None)
    if anchored_515f and gg2 is None:
        flags.append("too-short-for-gg2")
    return {
        "region": fits[0] if fits else "",
        "paired": paired,
        "primer_present": primer_present,
        "orient_primer": primer_present,
        "primer": fwd or "",
        "trim90": trim90,
        "trim": gg2 or trim90,
        "gg2_v4": gg2 is not None,
        "merge_rate": None if merge is None else round(merge, 3),
        "flags": ",".join(flags),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("studies", nargs="+")
    p.add_argument("--manifest", default=MANIFEST)
    p.add_argument("--reads", type=int, default=20000)
    p.add_argument("--out")
    p.add_argument("--memory", default="4GB")
    a = p.parse_args()

    con = duckdb.connect(config={"allow_unsigned_extensions": "true", "memory_limit": a.memory})
    con.execute(f"INSTALL miint FROM '{MIINT_REPO}'; LOAD miint")
    runs = con.execute(
        "SELECT study_accession, run_accession, prep_protocol, platform, path "
        "FROM read_csv($m, delim='\t', header=true) WHERE study_accession IN (SELECT unnest($s)) "
        "ORDER BY 1, 2", {"m": a.manifest, "s": a.studies}).fetchall()
    missing = set(a.studies) - {r[0] for r in runs}
    if missing:
        sys.exit(f"not in manifest (not imported yet?): {' '.join(sorted(missing))}")

    rows = []
    for study, run, protocol, platform, path in runs:
        per_region, merge = check_run(con, path, a.reads)
        rows.append({"study": study, "run": run, "protocol": protocol, "platform": platform,
                     **verdict(per_region, merge)})

    runs_df = pd.DataFrame(rows)
    con.register("runs", runs_df)
    runs_df.to_csv((a.out or "amplicon_check.tsv").removesuffix(".tsv") + ".runs.tsv", sep="\t", index=False)
    # study summary: the args to submit, plus any disagreement across runs flagged for review.
    study = con.sql("""
        SELECT study, count(*) AS runs,
               mode(region) AS region, count(DISTINCT region) > 1 AS region_mixed,
               bool_or(paired) AS paired,
               mode(primer) AS primer, bool_and(orient_primer) AS orient_primer,
               count(DISTINCT orient_primer) > 1 AS primer_mixed,
               min(trim) AS trim, bool_and(gg2_v4) AS gg2_v4,
               round(avg(merge_rate), 3) AS merge_rate,
               string_agg(DISTINCT nullif(flags, ''), ';') AS flags
        FROM runs GROUP BY study ORDER BY study""")
    study.show(max_width=250)
    if a.out:
        study.write_csv(a.out, sep="\t")


if __name__ == "__main__":
    main()
