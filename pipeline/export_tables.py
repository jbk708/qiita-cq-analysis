"""Export the pilot's combined tables: metagenomic and 16S V4 feature tables plus per-sample metadata for each.

Reads the dev stack directly (its Postgres app DB and DuckLake, read-only), so it runs where the stack's
Postgres is reachable (srun into the stack job). Sample IDs are ENA run accessions throughout.

  metag_feature_table.{tsv,biom}   species x run, sylph taxonomic abundance (%), from miint_profile.py's species.tsv
  16s_v4_feature_table.{tsv,biom}  ASV x run, Rapid 16S (amplicon workflow) counts; ASVs keyed by sequence
  16s_v4_feature_taxonomy.tsv      ASV, length, GG2 2024.09 exact-match taxonomy (blank when no match)
  metag_metadata.tsv   one row per profiled run: run/study fields plus the harmonized sample metadata
  16s_metadata.tsv     the same for every amplicon run, with in_v4_feature_table and the region check
Usage: export_tables.py <out_dir> --species <species.tsv> --v4-studies PRJ,... [--amplicon-check study.runs.tsv]
"""
import argparse
import os

import duckdb

root = os.environ.get("QDEV_ROOT", f"/ddn_scratch/{os.environ['USER']}/qiita-pilot/qiita-dev")
p = argparse.ArgumentParser()
p.add_argument("out")
p.add_argument("--species", required=True)
p.add_argument("--v4-studies", required=True, help="comma-separated study accessions processed with Rapid 16S")
p.add_argument("--amplicon-check", help="amplicon_check.py per-run output, adds region/primer/trim to 16s_metadata")
p.add_argument("--gg2-taxonomy", default="/databases/gg/2024.09/2024.09.taxonomy.asv.tsv.gz")
p.add_argument("--pg", default="dbname=qiita host=localhost port=55432")
p.add_argument("--lake", default="dbname=qiita_ducklake host=localhost port=55432")
a = p.parse_args()
os.makedirs(a.out, exist_ok=True)
v4 = [s for s in a.v4_studies.split(",") if s]

c = duckdb.connect(config={"allow_unsigned_extensions": "true", "memory_limit": "24GB", "threads": 8,
                           "extension_directory": f"{a.out}/.duckdb-ext"})
c.execute("INSTALL postgres; INSTALL ducklake; LOAD postgres; LOAD ducklake")
c.execute(f"LOAD '{root}/duckdb-ext/v1.5.4/linux_amd64/miint.duckdb_extension'")
c.execute(f"ATTACH '{a.pg}' AS pg (TYPE postgres, READ_ONLY)")
c.execute(f"ATTACH 'ducklake:postgres:{a.lake}' AS lake (DATA_PATH '{root}/persistent/ducklake', READ_ONLY)")

# run-level sample frame: every ENA-imported prep_sample
c.execute("""CREATE TABLE runs AS
  SELECT ss.ena_run_accession AS sample_name, ps.idx AS prep_sample_idx, ps.biosample_idx,
         coalesce(st.bioproject_accession, st.ena_study_accession) AS study_accession, st.ena_study_accession AS secondary_study_accession,
         st.title AS study_title, ss.ena_experiment_accession AS experiment_accession, b.biosample_accession,
         b.ena_sample_accession, pp.name AS prep_protocol, sr.platform
  FROM pg.qiita.sequenced_sample ss
  JOIN pg.qiita.prep_sample ps ON ps.idx = ss.prep_sample_idx
  JOIN pg.qiita.prep_protocol pp ON pp.idx = ps.prep_protocol_idx
  JOIN pg.qiita.biosample b ON b.idx = ps.biosample_idx
  JOIN pg.qiita.sequenced_pool sp ON sp.idx = ss.sequenced_pool_idx
  JOIN pg.qiita.sequencing_run sr ON sr.idx = sp.sequencing_run_idx
  JOIN pg.qiita.prep_sample_to_study pst ON pst.prep_sample_idx = ps.idx
  JOIN pg.qiita.study st ON st.idx = pst.study_idx
  WHERE ss.ena_run_accession IS NOT NULL""")

VALUE = """coalesce(m.value_text, m.value_numeric::VARCHAR, m.value_boolean::VARCHAR, m.value_date::VARCHAR,
                    tt.label, 'missing: ' || mr.name)"""
c.execute(f"""CREATE TABLE meta_long AS
  SELECT r.sample_name, coalesce(gf.internal_name, sf.display_name) AS field, {VALUE} AS value
  FROM runs r JOIN pg.qiita.biosample_metadata m ON m.biosample_idx = r.biosample_idx
  JOIN pg.qiita.biosample_study_field sf ON sf.idx = m.biosample_study_field_idx
  LEFT JOIN pg.qiita.biosample_global_field gf ON gf.idx = sf.biosample_global_field_idx
  LEFT JOIN pg.qiita.terminology_term tt ON tt.idx = m.value_terminology_term_idx
  LEFT JOIN pg.qiita.missing_value_reason mr ON mr.idx = m.value_missing_reason_idx
  UNION ALL
  SELECT r.sample_name, coalesce(gf.internal_name, sf.display_name), {VALUE}
  FROM runs r JOIN pg.qiita.prep_sample_metadata m ON m.prep_sample_idx = r.prep_sample_idx
  JOIN pg.qiita.prep_sample_study_field sf ON sf.idx = m.prep_sample_study_field_idx
  LEFT JOIN pg.qiita.prep_sample_global_field gf ON gf.idx = sf.prep_sample_global_field_idx
  LEFT JOIN pg.qiita.terminology_term tt ON tt.idx = m.value_terminology_term_idx
  LEFT JOIN pg.qiita.missing_value_reason mr ON mr.idx = m.value_missing_reason_idx""")


def write_metadata(path, samples_sql, extra_join="", extra_cols=""):
    c.execute(f"CREATE OR REPLACE TEMP TABLE s AS {samples_sql}")
    fields = [f for (f,) in c.execute("SELECT DISTINCT field FROM meta_long WHERE sample_name IN (SELECT sample_name FROM s) "
                                      "AND field IS NOT NULL ORDER BY 1").fetchall()]
    c.execute("CREATE OR REPLACE TEMP TABLE w AS PIVOT (SELECT * FROM meta_long WHERE sample_name IN (SELECT sample_name FROM s)) "
              "ON field USING first(value) GROUP BY sample_name")
    cols = ", ".join(f'w."{f}"' for f in fields) or "NULL AS no_metadata"
    c.execute(f"COPY (SELECT s.*{extra_cols}, {cols} FROM s LEFT JOIN w USING (sample_name) {extra_join} ORDER BY study_accession, sample_name) "
              f"TO '{path}' (DELIMITER '\t', HEADER)")
    n = c.execute(f"SELECT count(*) FROM read_csv('{path}', delim='\t', header=true)").fetchone()[0]
    print(f"{os.path.basename(path)}: {n} samples, {len(fields)} metadata fields")


# metagenomics
c.execute(f"CREATE TABLE metag AS SELECT run AS sample_name, species, any_value(lineage) lineage, sum(taxonomic_abundance) value "
          f"FROM read_csv('{a.species}', delim='\t', header=true) GROUP BY ALL")
c.execute(f"""COPY (PIVOT (SELECT species AS feature_id, lineage, sample_name, value FROM metag)
                   ON sample_name USING sum(value) GROUP BY feature_id, lineage ORDER BY feature_id)
              TO '{a.out}/metag_feature_table.tsv' (DELIMITER '\t', HEADER)""")
c.execute(f"COPY (SELECT species AS feature_id, sample_name AS sample_id, value FROM metag) "
          f"TO '{a.out}/metag_feature_table.biom' (FORMAT BIOM, COMPRESSION 'gzip', ID 'qiita-pilot-metag-sylph-gtdb-r220')")
print("metag_feature_table:", c.execute("SELECT count(DISTINCT species), count(DISTINCT sample_name) FROM metag").fetchone(), "(species, samples)")
write_metadata(f"{a.out}/metag_metadata.tsv", "SELECT * FROM runs WHERE sample_name IN (SELECT sample_name FROM metag)")

# 16S V4 (Rapid 16S): latest processing per prep_sample, ASV bytes from the chunk table
c.execute("""CREATE TABLE asv AS SELECT feature_idx, string_agg(chunk_data, '' ORDER BY chunk_index) AS sequence
             FROM lake.amplicon_sequence_chunks GROUP BY feature_idx""")
c.execute(f"""CREATE TABLE v4 AS
  WITH m AS (SELECT * FROM lake.amplicon_membership
             QUALIFY processing_idx = max(processing_idx) OVER (PARTITION BY prep_sample_idx))
  SELECT r.sample_name, r.study_accession, x.sequence, m.count
  FROM m JOIN runs r USING (prep_sample_idx) JOIN asv x USING (feature_idx)
  WHERE r.study_accession IN (SELECT unnest($v4))""", {"v4": v4})
c.execute(f"""COPY (PIVOT (SELECT sequence AS feature_id, sample_name, count FROM v4)
                   ON sample_name USING sum(count) GROUP BY feature_id ORDER BY feature_id)
              TO '{a.out}/16s_v4_feature_table.tsv' (DELIMITER '\t', HEADER)""")
c.execute(f"COPY (SELECT sequence AS feature_id, sample_name AS sample_id, count::DOUBLE AS value FROM v4) "
          f"TO '{a.out}/16s_v4_feature_table.biom' (FORMAT BIOM, COMPRESSION 'gzip', ID 'qiita-pilot-16s-v4-rapid16s')")
c.execute(f"""COPY (SELECT DISTINCT v.sequence AS feature_id, length(v.sequence) AS length_bp, g.Taxon AS gg2_2024_09_taxonomy
                    FROM v4 v LEFT JOIN read_csv('{a.gg2_taxonomy}', delim='\t', header=true) g ON g."Feature ID" = v.sequence
                    ORDER BY 1) TO '{a.out}/16s_v4_feature_taxonomy.tsv' (DELIMITER '\t', HEADER)""")
print("16s_v4_feature_table:", c.execute("SELECT count(DISTINCT sequence), count(DISTINCT sample_name), sum(count) FROM v4").fetchone(),
      "(ASVs, samples, reads)")
print("GG2 exact matches:", c.execute(f"SELECT count(*) FILTER (WHERE gg2_2024_09_taxonomy IS NOT NULL), count(*) "
                                      f"FROM read_csv('{a.out}/16s_v4_feature_taxonomy.tsv', delim='\t', header=true)").fetchone())
extra, extra_cols = "", ""
if a.amplicon_check:
    c.execute(f"CREATE TABLE chk AS SELECT run AS sample_name, region AS amplicon_region, primer AS amplicon_fwd_primer, "
              f"orient_primer AS amplicon_primer_in_reads, trim AS amplicon_trim, flags AS amplicon_check_flags "
              f"FROM read_csv('{a.amplicon_check}', delim='\t', header=true)")
    extra = "LEFT JOIN chk USING (sample_name)"
    extra_cols = ", chk.* EXCLUDE (sample_name)"
write_metadata(f"{a.out}/16s_metadata.tsv",
               "SELECT *, sample_name IN (SELECT sample_name FROM v4) AS in_v4_feature_table FROM runs "
               "WHERE prep_protocol LIKE '%amplicon'", extra, extra_cols)
