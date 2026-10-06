# qiita-cq-analysis

Study accession lists for batch ENA imports into [Qiita](https://github.com/the-miint/Qiita)
(`POST /api/v1/ena-import-batch`), sized from the ENA portal filereport, plus the scripts used to run the
import pilot and process the imported reads.

## Lists

- [`lists/targets.tsv`](lists/targets.tsv): target studies with their tranche and 16S primer region.
- [`lists/references.tsv`](lists/references.tsv): reference-panel candidates (EMP500, HMP, mouse gut, EMP 16S).
- [`lists/normalization_report.tsv`](lists/normalization_report.tsv): how the source accessions were kept, mapped to a study, or dropped.

Context and the first tranche: [the-miint/Qiita#601](https://github.com/the-miint/Qiita/issues/601).

## Tranches

The tranches are provisional until the pilot measures real throughput. The current sizing assumes:

- 3.2 MB/s per study download (4 connections at about 0.8 MB/s each), about 276 GB per 24 h job.
- 6 studies downloading at once, about 1.66 TB/day per tranche.

Studies over 200 GB go in their own lane, one per tranche. They're resubmitted in each following tranche until the stored read counts match ENA.

## Rebuilding

```
python3 build_lists.py --csv <screen.csv>
python3 -m unittest tests_build_lists
```

The script uses the standard library only. It caches raw ENA responses under `cache/` (not committed) and writes the full candidate pool locally. Only the curated lists are committed.

## Pilot scripts (barnacle)

Scripts used with a personal Qiita dev stack ([jbk708/qiita-barnacle-dev](https://github.com/jbk708/qiita-barnacle-dev))
for the import pilot and for processing the imported reads outside Qiita. Paths come from [`config.sh`](config.sh)
(`QDEV_ROOT`, `QDEV_SHARE`, `QDEV_CONDA`, `QDEV_CLI`); override them in the environment.

| Script | Does |
|---|---|
| `pilot/watch.sh <batch_idx>`, `pilot/progress.sh` | status of an ENA import batch and per-ticket download progress |
| `pilot/portal-watch.sh` | waits for the ENA Portal API to answer, then submits a study list once |
| `pilot/wave-chain.sh <list>…` | submits study lists one after another, each once the previous one's downloads finish |
| `pilot/auto-redrive.sh <first_batch>` | redrives tickets that failed on a transient ENA error, in rounds (max 6 each) |
| `pilot/profile-watch.sh <out>` | rolling profiling: when a study's downloads all complete, submits its new metagenomic runs to `submit_profile.sh` |
| `pilot/amplicon-run.sh <plan.tsv>` | runs the Rapid 16S `amplicon` workflow pool by pool from a plan file |
| `pilot/redrive-seq.sh <idx>…` | redrives failed download tickets one at a time, stopping on a failure |
| `pilot/verify_counts.py PRJ…` | compares stored read counts (staging and lake) with ENA `read_count` |
| `pilot/share-refresh.sh` | rebuilds the shared read manifest and permissions, and copies `pipeline/` and [`docs/reads-share.md`](docs/reads-share.md) into the share |
| `pipeline/submit_profile.sh <out> PRJ…` | **pilot profiling:** one Slurm job running `miint_profile.py` |
| `pipeline/miint_profile.py <out> PRJ…` | miint `sylph_profile` straight from each run's parquet (no FASTQ, no host depletion) vs GTDB r220; resumable; writes `profile.tsv` and `species.tsv` |
| `pipeline/submit_runs.sh <out> PRJ…` | Slurm array over a study's metagenomic runs → `deplete_sketch.sh` |
| `pipeline/deplete_sketch.sh` | parquet → miint FASTQ stream → deacon (panhuman-1) → sylph sketch |
| `pipeline/profile.sh <out>` | sylph profile vs GTDB r232, then sylph-tax |
| `pipeline/export_tables.py <out> …` | combined metagenomic and 16S V4 feature tables (TSV + BIOM) and per-sample metadata, read from the dev stack |
| `pipeline/amplicon_check.py PRJ…` | per-study 16S region, primers, trim and multiplexing → `amplicon` workflow args |
| `pipeline/jupyter.sbatch`, `pipeline/qiita_reads.ipynb` | JupyterLab on the `jupyter` partition |

miint embeds a sylph 0.9 fork, so `miint_profile.py` needs the GTDB r220 `.syldb` (sylph 1.0's r232 `.syl2db` does not load).
The deacon + sylph-CLI scripts stay for benchmarking host-depleted against raw profiles.

Databases (panhuman-1, the GTDB r220 `.syldb` and r232 `.syl2db` sylph dbs, sylph-tax) are expected in `$QDEV_SHARE/db/`; the conda env
(DuckDB 1.5.4, deacon, sylph, sylph-tax, JupyterLab) in `$QDEV_CONDA`.
