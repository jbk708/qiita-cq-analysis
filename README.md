# qiita-ena-manifests

Study accession lists for batch ENA imports into [Qiita](https://github.com/the-miint/Qiita)
(`POST /api/v1/ena-import-batch`), sized from the ENA portal filereport.

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
