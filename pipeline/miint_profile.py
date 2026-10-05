"""Profile metagenomic runs with miint sylph_profile straight from their Qiita parquet (no FASTQ, no host depletion).

miint embeds a sylph 0.9 fork, so the database must be a .syldb (GTDB r220), not sylph 1.0's .syl2db.
Runs already in <out_dir>/profile.tsv are skipped. Writes profile.tsv (genome level, with lineage) and
species.tsv (run, species, lineage, abundances).
Usage: miint_profile.py <out_dir> PRJ... [--batch 10] [--threads 16] [--memory 48GB] [--estimate-unknown]
"""
import argparse
import os

import duckdb

share = os.environ.get("QDEV_SHARE", f"/ddn_scratch/{os.environ['USER']}/qiita-pilot/qiita-dev-share")
root = os.environ.get("QDEV_ROOT", f"/ddn_scratch/{os.environ['USER']}/qiita-pilot/qiita-dev")
p = argparse.ArgumentParser()
p.add_argument("out")
p.add_argument("studies", nargs="+")
p.add_argument("--db", default=f"{share}/db/gtdb-r220-c200-dbv1.syldb")
p.add_argument("--taxonomy", default=f"{share}/db/sylph-tax/gtdb_r220_metadata.tsv.gz")
p.add_argument("--batch", type=int, default=10, help="runs per sylph_profile call")
p.add_argument("--threads", type=int, default=16)
p.add_argument("--memory", default="48GB")
p.add_argument("--estimate-unknown", action="store_true")
a = p.parse_args()

os.makedirs(f"{a.out}/tmp", exist_ok=True)
prof_path, sp_path = f"{a.out}/profile.tsv", f"{a.out}/species.tsv"
c = duckdb.connect(config={"allow_unsigned_extensions": "true", "extension_directory": f"{root}/duckdb-ext",
                           "memory_limit": a.memory, "threads": a.threads, "temp_directory": f"{a.out}/tmp"})
c.execute("LOAD miint")
c.execute(f"CREATE TABLE tax AS SELECT column0 acc, column1 lineage, regexp_extract(column1, 's__([^;]*)$', 1) species "
          f"FROM read_csv('{a.taxonomy}', delim='\t', header=false)")
runs = c.execute(f"""SELECT study_accession, run_accession, path FROM read_csv('{share}/manifest.tsv', delim='\t', header=true)
                     WHERE study_accession IN (SELECT unnest($s)) AND prep_protocol LIKE '%metagenomics'
                     ORDER BY 1, 2""", {"s": a.studies}).fetchall()
done = set()
if os.path.exists(prof_path):
    done = {r for (r,) in c.execute(f"SELECT DISTINCT run FROM read_csv('{prof_path}', delim='\t', header=true)").fetchall()}
todo = [r for r in runs if r[1] not in done]
print(f"{len(runs)} metagenomic runs, {len(done)} already profiled, {len(todo)} to do", flush=True)

for i in range(0, len(todo), a.batch):
    batch = todo[i:i + a.batch]
    c.execute("CREATE OR REPLACE VIEW reads AS " + " UNION ALL ".join(
        f"SELECT '{s}' AS study, '{r}' AS run, read_id, sequence1, sequence2 FROM read_parquet('{path}')"
        for s, r, path in batch))
    c.execute(f"""CREATE OR REPLACE TABLE prof AS
        SELECT s.study, p.sample_id AS run, p.* EXCLUDE (sample_id), t.lineage, t.species
        FROM sylph_profile('reads', '{a.db}', sample_id := 'run', estimate_unknown := {str(a.estimate_unknown).lower()}) p
        LEFT JOIN (SELECT DISTINCT study, run FROM (VALUES {", ".join(f"('{s}', '{r}')" for s, r, _ in batch)}) v(study, run)) s
               ON s.run = p.sample_id
        LEFT JOIN tax t ON t.acc = regexp_extract(p.genome_name, '(GC[AF]_[0-9]+\\.[0-9]+)', 1)""")
    missing = {r for _, r, _ in batch} - {r for (r,) in c.execute("SELECT DISTINCT run FROM prof").fetchall()}
    header = not os.path.exists(prof_path)
    c.execute(f"COPY prof TO '{a.out}/.batch.tsv' (DELIMITER '\t', HEADER {str(header).lower()})")
    with open(prof_path, "a") as out, open(f"{a.out}/.batch.tsv") as b:
        out.write(b.read())
    print(f"batch {i // a.batch + 1}: {len(batch)} runs, {c.execute('SELECT count(*) FROM prof').fetchone()[0]} genome rows"
          + (f"; no genomes passed for {sorted(missing)}" if missing else ""), flush=True)
os.remove(f"{a.out}/.batch.tsv") if os.path.exists(f"{a.out}/.batch.tsv") else None

c.execute(f"""COPY (SELECT study, run, coalesce(species, genome_name) species, any_value(lineage) lineage,
                    sum(taxonomic_abundance) taxonomic_abundance, sum(sequence_abundance) sequence_abundance,
                    count(*) genomes
                    FROM read_csv('{prof_path}', delim='\t', header=true) GROUP BY ALL ORDER BY run, taxonomic_abundance DESC)
              TO '{sp_path}' (DELIMITER '\t', HEADER)""")
c.sql(f"""WITH s AS (SELECT * FROM read_csv('{sp_path}', delim='\t', header=true)),
               r AS (SELECT study, run, count(*) n FROM s GROUP BY ALL)
          SELECT study, count(*) runs, (SELECT count(DISTINCT species) FROM s WHERE s.study = r.study) species,
                 round(median(n), 0) median_species_per_run, min(n) min_species, max(n) max_species
          FROM r GROUP BY study ORDER BY study""").show()
