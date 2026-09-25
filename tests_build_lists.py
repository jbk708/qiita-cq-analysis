"""Unit tests for build_lists.py. Pure-logic only, no network."""

from __future__ import annotations

import unittest

import build_lists as bl


class TestClassifyToken(unittest.TestCase):
    def test_study_prefixes_kept(self):
        for tok in ["PRJNA810087", "PRJEB29237", "PRJDB1234", "ERP012345", "SRP012345", "DRP001234"]:
            with self.subTest(tok=tok):
                c = bl.classify_token(tok)
                self.assertEqual(c.category, "study")

    def test_sample_prefixes_mapped(self):
        for tok in ["SAMN12345678", "SAME1234567", "SAMD00123456", "ERS1234567", "DRS001948"]:
            with self.subTest(tok=tok):
                c = bl.classify_token(tok)
                self.assertEqual(c.category, "sample")
                self.assertEqual(c.query_field, "sample_accession")

    def test_run_prefixes_mapped(self):
        for tok in ["SRR1234567", "ERR1234567", "DRR1234567"]:
            c = bl.classify_token(tok)
            self.assertEqual(c.category, "run")
            self.assertEqual(c.query_field, "run_accession")

    def test_experiment_prefixes_mapped(self):
        for tok in ["SRX1234567", "ERX1234567", "DRX275800"]:
            c = bl.classify_token(tok)
            self.assertEqual(c.category, "experiment")
            self.assertEqual(c.query_field, "experiment_accession")

    def test_submission_prefixes_mapped(self):
        for tok in ["SRA000001", "ERA3198568", "DRA000316"]:
            c = bl.classify_token(tok)
            self.assertEqual(c.category, "submission")
            self.assertEqual(c.query_field, "submission_accession")

    def test_drop_ena_not_found(self):
        c = bl.classify_token("ENA_NOT_FOUND")
        self.assertEqual(c.category, "drop")
        self.assertEqual(c.reason, "ena_not_found")

    def test_drop_geo(self):
        c = bl.classify_token("GSE154918")
        self.assertEqual((c.category, c.reason), ("drop", "geo"))

    def test_drop_cngb_gsa_family(self):
        for tok in ["CRA009631", "HRA000001", "PRJCA001234", "CNP0001234", "OEP001234",
                    "CNS2022", "subCRA028124"]:
            with self.subTest(tok=tok):
                c = bl.classify_token(tok)
                self.assertEqual((c.category, c.reason), ("drop", "cngb_gsa"))

    def test_drop_restricted_access(self):
        for tok in ["EGAS00001000001", "EGAD50000001154", "phs000001"]:
            with self.subTest(tok=tok):
                c = bl.classify_token(tok)
                self.assertEqual((c.category, c.reason), ("drop", "restricted_access"))

    def test_drop_metabolomics(self):
        c = bl.classify_token("MTBLS1971")
        self.assertEqual((c.category, c.reason), ("drop", "metabolomics"))

    def test_drop_arrayexpress(self):
        c = bl.classify_token("E-MTAB-6940")
        self.assertEqual((c.category, c.reason), ("drop", "arrayexpress"))

    def test_drop_nmdc(self):
        c = bl.classify_token("NMDC12345")
        self.assertEqual((c.category, c.reason), ("drop", "nmdc"))

    def test_drop_assembly(self):
        for tok in ["GCF_000001405.33", "GCA_900066225"]:
            with self.subTest(tok=tok):
                c = bl.classify_token(tok)
                self.assertEqual((c.category, c.reason), ("drop", "assembly"))

    def test_drop_genbank_nucleotide(self):
        for tok in ["MH108987", "GU413083", "MT780937", "AB971822", "AF183403"]:
            with self.subTest(tok=tok):
                c = bl.classify_token(tok)
                self.assertEqual((c.category, c.reason), ("drop", "genbank_nucleotide"))

    def test_mtbls_not_caught_by_genbank_regex(self):
        # MTBLS1971 starts with two letters (MT) but the remainder isn't a bare
        # digit run -- the explicit mtbls rule must win, not the genbank regex.
        c = bl.classify_token("MTBLS1971")
        self.assertEqual(c.reason, "metabolomics")

    def test_drop_other_non_insdc(self):
        c = bl.classify_token("PXD062630")
        self.assertEqual((c.category, c.reason), ("drop", "other_non_insdc"))

    def test_drop_truncated_or_malformed(self):
        c = bl.classify_token("(+22 more)")
        self.assertEqual((c.category, c.reason), ("drop", "truncated_or_malformed"))

    def test_drop_unrecognized_catch_all(self):
        c = bl.classify_token("ZZZ12345")
        self.assertEqual((c.category, c.reason), ("drop", "unrecognized_prefix"))

    def test_empty_token_raises(self):
        with self.assertRaises(ValueError):
            bl.classify_token("")


class TestSumFastqBytes(unittest.TestCase):
    def test_paired(self):
        self.assertEqual(bl.sum_fastq_bytes("12345;67890"), 80235)

    def test_single(self):
        self.assertEqual(bl.sum_fastq_bytes("12345"), 12345)

    def test_empty(self):
        self.assertEqual(bl.sum_fastq_bytes(""), 0)

    def test_trailing_semicolon_and_spaces(self):
        self.assertEqual(bl.sum_fastq_bytes(" 100 ; 200 ;"), 300)

    def test_non_numeric_raises(self):
        with self.assertRaises(ValueError):
            bl.sum_fastq_bytes("100;notanumber")


class TestPlatformMapping(unittest.TestCase):
    def test_known_platforms(self):
        cases = {
            "ILLUMINA": "ILLUMINA",
            "BGISEQ": "DNBSEQ",
            "DNBSEQ": "DNBSEQ",
            "PACBIO_SMRT": "PACBIO_SMRT",
            "OXFORD_NANOPORE": "OXFORD_NANOPORE",
            "LS454": "LS454",
            "ION_TORRENT": "ION_TORRENT",
            "COMPLETE_GENOMICS": "COMPLETE_GENOMICS",
        }
        for raw, mapped in cases.items():
            with self.subTest(raw=raw):
                self.assertEqual(bl.map_platform(raw), mapped)

    def test_case_insensitive(self):
        self.assertEqual(bl.map_platform("illumina"), "ILLUMINA")

    def test_unmapped_returns_none(self):
        self.assertIsNone(bl.map_platform("CAPILLARY"))
        self.assertIsNone(bl.map_platform(""))
        self.assertIsNone(bl.map_platform(None))


class TestProtocolCategory(unittest.TestCase):
    """Mirrors ena_import/protocol_mapping.py's dispatch (category only, no
    platform-driven read-length bucket)."""

    def test_amplicon(self):
        self.assertEqual(bl.protocol_category_for_run("AMPLICON", "GENOMIC"), "amplicon")

    def test_wgs_strategy(self):
        self.assertEqual(bl.protocol_category_for_run("WGS", "GENOMIC"), "metagenomics")

    def test_metagenomic_source(self):
        self.assertEqual(bl.protocol_category_for_run("OTHER", "METAGENOMIC"), "metagenomics")

    def test_transcriptomics(self):
        self.assertEqual(bl.protocol_category_for_run("RNA-SEQ", "TRANSCRIPTOMIC"), "transcriptomics")
        self.assertEqual(
            bl.protocol_category_for_run("OTHER", "METATRANSCRIPTOMIC"), "transcriptomics"
        )

    def test_unmapped(self):
        self.assertEqual(bl.protocol_category_for_run("CHIP-SEQ", "GENOMIC"), "unmapped")


def _blank_metrics(**overrides) -> bl.StudyMetrics:
    m = bl.StudyMetrics()
    for k, v in overrides.items():
        setattr(m, k, v)
    return m


class TestDeriveStudyFlags(unittest.TestCase):
    def test_not_found_short_circuits(self):
        flags, failure = bl.derive_study_flags(_blank_metrics(), "not_found")
        self.assertEqual(flags, ["not_found"])
        self.assertTrue(failure)

    def test_ena_error_short_circuits(self):
        flags, failure = bl.derive_study_flags(_blank_metrics(), "ena_error")
        self.assertEqual(flags, ["ena_error"])
        self.assertTrue(failure)

    def test_no_public_runs(self):
        flags, failure = bl.derive_study_flags(_blank_metrics(n_runs=0), "resolved")
        self.assertEqual(flags, ["no_public_runs"])
        self.assertTrue(failure)

    def test_no_fastq_all_runs_is_failure(self):
        m = _blank_metrics(n_runs=5, n_runs_with_fastq=0, no_fastq_runs=5)
        flags, failure = bl.derive_study_flags(m, "resolved")
        self.assertIn("no_fastq:5", flags)
        self.assertTrue(failure)

    def test_no_fastq_partial_is_not_failure(self):
        m = _blank_metrics(n_runs=5, n_runs_with_fastq=3, no_fastq_runs=2)
        flags, failure = bl.derive_study_flags(m, "resolved")
        self.assertIn("no_fastq:2", flags)
        self.assertFalse(failure)

    def test_unmapped_platform_all_runs_is_failure(self):
        m = _blank_metrics(n_runs=3, n_runs_with_fastq=3,
                            unmapped_platform_runs={"CAPILLARY": 3})
        flags, failure = bl.derive_study_flags(m, "resolved")
        self.assertIn("unmapped_platform:CAPILLARY:3", flags)
        self.assertTrue(failure)

    def test_unmapped_platform_partial_is_not_failure(self):
        m = _blank_metrics(n_runs=5, n_runs_with_fastq=5,
                            unmapped_platform_runs={"CAPILLARY": 2})
        flags, failure = bl.derive_study_flags(m, "resolved")
        self.assertIn("unmapped_platform:CAPILLARY:2", flags)
        self.assertIn("unmapped_platform partial", flags)
        self.assertFalse(failure)

    def test_unmapped_strategy_never_a_failure(self):
        m = _blank_metrics(n_runs=3, n_runs_with_fastq=3, unmapped_strategy_runs=3)
        flags, failure = bl.derive_study_flags(m, "resolved")
        self.assertIn("unmapped_strategy:3", flags)
        self.assertFalse(failure)

    def test_amplicon_and_mixed_assay(self):
        m = _blank_metrics(n_runs=4, n_runs_with_fastq=4, amplicon_runs=2, metagenomics_runs=2)
        flags, failure = bl.derive_study_flags(m, "resolved")
        self.assertIn("amplicon:2", flags)
        self.assertIn("mixed_assay", flags)
        self.assertFalse(failure)

    def test_amplicon_and_transcriptomics_is_mixed(self):
        m = _blank_metrics(n_runs=4, n_runs_with_fastq=4, amplicon_runs=2, transcriptomics_runs=2)
        flags, _ = bl.derive_study_flags(m, "resolved")
        self.assertIn("mixed_assay", flags)

    def test_genomic_source(self):
        m = _blank_metrics(n_runs=2, n_runs_with_fastq=2, metagenomics_runs=2, genomic_source_runs=1)
        flags, failure = bl.derive_study_flags(m, "resolved")
        self.assertIn("genomic_source:1", flags)
        self.assertFalse(failure)


class TestStudyCoarseAssay(unittest.TestCase):
    def test_pure_amplicon(self):
        m = _blank_metrics(amplicon_runs=5)
        self.assertEqual(bl.study_coarse_assay(m), "AMPLICON")

    def test_pure_metagenomics(self):
        m = _blank_metrics(metagenomics_runs=5)
        self.assertEqual(bl.study_coarse_assay(m), "WGS")

    def test_pure_transcriptomics(self):
        m = _blank_metrics(transcriptomics_runs=5)
        self.assertEqual(bl.study_coarse_assay(m), "RNA-SEQ")

    def test_mixed(self):
        m = _blank_metrics(amplicon_runs=2, metagenomics_runs=2)
        self.assertEqual(bl.study_coarse_assay(m), "MIXED")

    def test_single_assay_wgs_tolerates_genomic_source_informational_flag(self):
        # genomic_source:N is informational (strategy==WGS runs not confirmed
        # METAGENOMIC by source) -- it does not itself break the WGS/METAGENOMIC
        # category bucket that single-assay candidacy checks.
        m = _blank_metrics(n_runs=5, metagenomics_runs=5, genomic_source_runs=1)
        self.assertTrue(bl.is_single_assay(m))
        self.assertEqual(bl.study_coarse_assay(m), "WGS")

    def test_mixed_is_not_single_assay(self):
        m = _blank_metrics(n_runs=4, amplicon_runs=2, metagenomics_runs=2)
        self.assertFalse(bl.is_single_assay(m))


class TestEstJobHoursAndPasses(unittest.TestCase):
    def test_zero(self):
        self.assertEqual(bl.est_job_hours(0), 0.0)
        self.assertEqual(bl.est_passes(0), 0)

    def test_one_pass_exact(self):
        self.assertAlmostEqual(bl.est_job_hours(bl.PASS_BYTES), 24.0, places=6)
        self.assertEqual(bl.est_passes(bl.PASS_BYTES), 1)

    def test_just_over_one_pass(self):
        self.assertEqual(bl.est_passes(bl.PASS_BYTES + 1), 2)


class TestSampleHintAndReview(unittest.TestCase):
    def test_fecal(self):
        hints = {"isolation_source": {"human stool sample"}}
        self.assertEqual(bl.sample_hint(hints), "fecal")

    def test_tissue(self):
        hints = {"tissue_type": {"colon biopsy"}}
        self.assertEqual(bl.sample_hint(hints), "tissue")

    def test_mixed(self):
        hints = {"isolation_source": {"stool"}, "tissue_type": {"tumour"}}
        self.assertEqual(bl.sample_hint(hints), "mixed")

    def test_unknown(self):
        hints = {"isolation_source": {"unspecified"}}
        self.assertEqual(bl.sample_hint(hints), "unknown")

    def test_needs_review_when_unknown(self):
        hints = {"isolation_source": {"unspecified"}}
        self.assertTrue(bl.needs_manual_review("WGS", hints, set()))

    def test_amplicon_without_16s_needs_review(self):
        hints = {"isolation_source": {"stool"}}
        self.assertTrue(bl.needs_manual_review("AMPLICON", hints, {"ITS"}))

    def test_amplicon_with_16s_no_review(self):
        hints = {"isolation_source": {"stool"}}
        self.assertFalse(bl.needs_manual_review("AMPLICON", hints, {"16S rRNA V4"}))

    def test_wgs_with_known_hint_no_review(self):
        hints = {"tissue_type": {"biopsy"}}
        self.assertFalse(bl.needs_manual_review("WGS", hints, set()))


class TestAssignTranches(unittest.TestCase):
    def test_small_studies_fit_in_one_tranche(self):
        studies = [("A", 1_000_000_000), ("B", 2_000_000_000), ("C", 500_000_000)]
        result = bl.assign_tranches(studies)
        self.assertEqual(len(result.tranches), 1)
        self.assertEqual(set(result.tranches[0]), {"A", "B", "C"})
        self.assertEqual(result.hold, [])

    def test_hold_bucket_for_zero_byte_studies(self):
        studies = [("A", 0), ("B", 1_000_000_000)]
        result = bl.assign_tranches(studies)
        self.assertIn("A", result.hold)
        self.assertNotIn("A", [s for t in result.tranches for s in t])

    def test_large_study_spans_multiple_consecutive_tranches(self):
        large_bytes = bl.PASS_BYTES * 2 + 1  # needs 3 passes
        studies = [("BIG", large_bytes)]
        result = bl.assign_tranches(studies)
        self.assertEqual(len(result.tranches), 3)
        for t in result.tranches:
            self.assertEqual(t, ["BIG"])

    def test_at_most_one_large_per_tranche(self):
        big_bytes = bl.PASS_BYTES + 1  # 2 passes each
        studies = [("BIG1", big_bytes), ("BIG2", big_bytes)]
        result = bl.assign_tranches(studies)
        for t in result.tranches:
            large_in_tranche = [s for s in t if s in ("BIG1", "BIG2")]
            self.assertLessEqual(len(large_in_tranche), 1)

    def test_largest_first_into_lightest_tranche(self):
        # Two tranches' worth of small studies; the balancer should keep
        # per-tranche totals roughly even rather than filling one first.
        studies = [("A", 900_000_000_000), ("B", 900_000_000_000), ("C", 900_000_000_000)]
        result = bl.assign_tranches(studies)
        totals = [sum(dict(studies)[s] for s in t) for t in result.tranches]
        self.assertLessEqual(max(totals) - min(totals), 900_000_000_000)


class TestReferenceTranche(unittest.TestCase):
    def test_unresolved_row_is_blank(self):
        self.assertEqual(bl.reference_tranche(None, None, False, 0), "")

    def test_explicit_override_wins(self):
        self.assertEqual(bl.reference_tranche("PRJEB42019", "local", False, 1_000), "local")
        self.assertEqual(bl.reference_tranche("PRJEB7759", "large", False, 1_000), "large")

    def test_no_fastq_failure_is_hold(self):
        self.assertEqual(bl.reference_tranche("PRJNA48479", None, True, 0), "hold")

    def test_large_by_size_when_no_override(self):
        self.assertEqual(
            bl.reference_tranche("PRJNA43021", None, False, bl.LARGE_THRESHOLD_BYTES + 1), "large"
        )

    def test_small_clean_study_is_blank(self):
        self.assertEqual(bl.reference_tranche("PRJNA1", None, False, 1_000_000_000), "")


def _candidate(study, coarse_assay, fastq_bytes, needs_manual_review=False, platform="ILLUMINA"):
    return {
        "study": study,
        "coarse_assay": coarse_assay,
        "fastq_bytes": fastq_bytes,
        "needs_manual_review": needs_manual_review,
        "platform": platform,
    }


class TestSuggestTranche1(unittest.TestCase):
    def test_size_band_helper(self):
        self.assertEqual(bl.size_band(1_000_000_000), "under_5gb")
        self.assertEqual(bl.size_band(20_000_000_000), "5_to_50gb")
        self.assertEqual(bl.size_band(100_000_000_000), "50_to_250gb")

    def test_spans_all_three_size_bands_when_available(self):
        # Plenty of tiny AMPLICON/WGS candidates plus one real candidate per
        # larger band. The old smallest-first-only implementation always
        # grabbed the tiniest amplicon/wgs pair and never advanced past
        # `under_5gb`, so this fails on that behavior.
        candidates = [
            _candidate("AMP_TINY_1", "AMPLICON", 1_00_000_000),
            _candidate("AMP_TINY_2", "AMPLICON", 1_20_000_000),
            _candidate("WGS_TINY_1", "WGS", 1_00_000_000),
            _candidate("WGS_TINY_2", "WGS", 1_20_000_000),
            _candidate("WGS_MED", "WGS", 20_000_000_000),
            _candidate("WGS_LARGE", "WGS", 100_000_000_000),
        ]
        result = bl.suggest_tranche_1(candidates)
        bands = {bl.size_band(c["fastq_bytes"]) for c in result}
        self.assertEqual(bands, {"under_5gb", "5_to_50gb", "50_to_250gb"})

    def test_includes_amplicon_and_wgs(self):
        candidates = [
            _candidate("AMP_SMALL", "AMPLICON", 1_000_000_000),
            _candidate("WGS_MED", "WGS", 20_000_000_000),
            _candidate("WGS_LARGE", "WGS", 100_000_000_000),
        ]
        result = bl.suggest_tranche_1(candidates)
        assays = {c["coarse_assay"] for c in result}
        self.assertIn("AMPLICON", assays)
        self.assertIn("WGS", assays)

    def test_prefers_non_illumina_when_available(self):
        candidates = [
            _candidate("AMP_SMALL", "AMPLICON", 1_000_000_000, platform="ILLUMINA"),
            _candidate("WGS_MED", "WGS", 20_000_000_000, platform="ILLUMINA"),
            _candidate("WGS_LARGE", "WGS", 100_000_000_000, platform="ILLUMINA"),
            _candidate("WGS_NANOPORE", "WGS", 30_000_000_000, platform="OXFORD_NANOPORE"),
        ]
        result = bl.suggest_tranche_1(candidates)
        platforms = {c["platform"] for c in result}
        self.assertIn("OXFORD_NANOPORE", platforms)

    def test_prefers_review_clean_candidates(self):
        candidates = [
            _candidate("AMP_CLEAN", "AMPLICON", 1_000_000_000, needs_manual_review=False),
            _candidate("AMP_REVIEW", "AMPLICON", 1_100_000_000, needs_manual_review=True),
            _candidate("WGS_MED", "WGS", 20_000_000_000),
            _candidate("WGS_LARGE", "WGS", 100_000_000_000),
        ]
        result = bl.suggest_tranche_1(candidates)
        picked_studies = {c["study"] for c in result}
        self.assertIn("AMP_CLEAN", picked_studies)
        self.assertNotIn("AMP_REVIEW", picked_studies)

    def test_stays_within_budget_and_size_bound(self):
        candidates = [
            _candidate("AMP_SMALL", "AMPLICON", 1_000_000_000),
            _candidate("WGS_MED", "WGS", 20_000_000_000),
            _candidate("WGS_LARGE", "WGS", 100_000_000_000),
        ]
        result = bl.suggest_tranche_1(candidates)
        self.assertLessEqual(sum(c["fastq_bytes"] for c in result), 300_000_000_000)
        self.assertLessEqual(len(result), 5)

    def test_excludes_given_accessions(self):
        candidates = [
            _candidate("AMP_SMALL", "AMPLICON", 1_000_000_000),
            _candidate("WGS_MED", "WGS", 20_000_000_000),
        ]
        result = bl.suggest_tranche_1(candidates, exclude=frozenset({"WGS_MED"}))
        self.assertNotIn("WGS_MED", {c["study"] for c in result})


if __name__ == "__main__":
    unittest.main()
