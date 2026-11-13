# ENA Portal shim

Lets Qiita's `ingest_ena_reads` fetch a study from files already on disk, for runs whose ENA copy has
been re-encoded or renamed (PRJEB76409 is stored under `per_sample_FASTQ/178353/`). It answers the two
things miint's `read_ena_sequences` asks ENA for: the Portal `read_run` search and the FASTQ files.

- `shim/ena_shim.py`: the server. Strict (default): a run missing from the manifest gets a 422, never a fallback to real ENA.
  Proxy (`--mode proxy`): ENA's own Portal row is passed through for such runs, so miint fetches them from ENA; a run ENA
  doesn't return is a 502 naming it.
- `shim/build_manifest.py`: matches ENA `submitted_bytes` to local files and writes the manifest.
- `shim/stage_sra.py`: stages runs from NCBI's SRA bucket and writes a per-study manifest; `--partial` keeps the runs that
  stage and lists the rest in `<study>.skipped`.
- `shim/native-python.sh`, `shim/scope_check.py`: sets `MIINT_ENA_PORTAL_URL_BASE` for `ingest_ena_reads` jobs whose runs are all in the manifest.
- `shim/ena-shim.sbatch`: runs the server as its own Slurm job.

miint reads the variable only from `8f79da38` (2026-10-22) on; older builds ignore it and go to real ENA.
The shim has no authentication and listens on all interfaces: serve public studies only.

## Manifest

TSV with the columns `run_accession, study_accession, sample_accession, experiment_accession, library_layout,
r1_path, r1_bytes, r1_md5, r2_path, r2_bytes, r2_md5` (extra columns are ignored). md5s are ENA's `submitted_md5`.
The server refuses to start on a duplicate run or path, a relative path, a missing or wrong-sized file,
a file without a 32-hex md5, or a layout that disagrees with `r2_path`.

```
curl -s 'https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession%3D%22PRJEB76409%22&fields=run_accession,sample_accession,experiment_accession,run_alias,library_layout,submitted_bytes,submitted_md5&format=tsv&limit=0' > PRJEB76409.tsv
python3 shim/build_manifest.py --ena-tsv PRJEB76409.tsv --qiita-root /qmounts/qiita_data --study PRJEB76409 --out PRJEB76409.manifest.tsv
```

The builder exits non-zero, naming each run, when a size matches no file or several, or a file is claimed twice.
Its only check is size; miint's md5 check at download catches a wrong match.

## Running

State lives in `/ddn_scratch/$USER/qiita-pilot/ena-shim/` (manifests and audit logs are not committed).

1. `sbatch shim/ena-shim.sbatch <manifest> [strict|proxy]` writes `portal.env` (the URL) and `portal.runs` (the manifest's runs) and audits to `<manifest>.audit.jsonl`.
2. `curl http://<node>:<port>/healthz` should print `ok`; check it from another node too.
3. `SLURM_NATIVE_PYTHON` in `compute-orchestrator.env` points at `shim/native-python.sh`. It does nothing unless `portal.env` exists and the job is `ingest_ena_reads`.
   Then `scope_check.py` compares the job's `ena_run_map` runs with `portal.runs`: all in exports the URL, none leaves the job on real ENA, a mix fails the job
   unless `portal.env` has `ENA_SHIM_MODE=proxy`. A job routed to the shim also gets `QIITA_DEV_ENA_PRIORITY_RUNS=portal.runs`,
   which a dev Qiita build uses to start the shim's runs first.
4. Check the audit log: one `/portal/api/search` per run, HEADs, 206 range reads, and no `"full_get": true`.

## Teardown

`portal.env` and `portal.runs` are never deleted by the server, so a stopped shim makes jobs fail with "connection refused"
instead of silently reaching real ENA. When the pilot is final, remove both files, then `scancel` the shim.

## Tests

```
cd shim && python3 -m unittest tests_build_manifest tests_ena_shim tests_scope_wrapper tests_stage_sra
SHIM_MIINT_PYTHON=<orchestrator venv python> MIINT_EXTENSION_DIRECTORY=<staged ext dir> python3 -m unittest tests_miint_integration
```

The integration test runs real miint against fixtures served on 127.0.0.1 (use `srun`, not the login node). It checks
range reads with no full GET, that a same-size one-base corruption fails on md5, and that an unmapped run fails with HTTP 422.
Qiita classifies that 422 as a transient md5 error, because its classifier matches "md5" in the echoed search URL
(`fastq_md5`); the test marks the BAD_INPUT expectation as an expected failure.
