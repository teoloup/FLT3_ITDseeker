"""Regression tests for ambiguous support, missing WT, and breakpoint validation."""
import logging
import unittest

import pandas as pd

from itdseeker.GMM_peaks import fit_gmm_itds, refine_peak_substructure_once
from itdseeker.Helper_functions import classify_read_support, validate_itd_supporting_reads
from itdseeker.Pairwise_aligment_toolkit import align_reads_multi_ref_parallel
from itdseeker.haplotype_split import assignments_to_result

WT = (
    "CTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCG"
    "TTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTT"
    "TAAAGATTTTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAAT"
    "TTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCT"
    "GAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTAAATTGCT"
)


def fit_length(length):
    return fit_gmm_itds(
        {f"r{i}": "A" * length for i in range(40)},
        min_gmm_fraction=.01, max_itds_detected=2,
        min_ggmm_peak_distance=10, max_peak_sd=5, assign_width_factor=2,
        assign_mode="manual", prob_threshold=.85, wt_amplicon_length=336,
        force_k=1,
    )


def cluster_fixture(alias, lengths):
    reads = pd.DataFrame([
        dict(read_id=f"r{i}", read_len=n, read_seq="A" * n,
             gmm_peak_alias=alias, is_ambiguous=False, strand="+")
        for i, n in enumerate(lengths)
    ])
    comps = pd.DataFrame([dict(peak_alias=alias, mean_bp=sum(lengths)/len(lengths),
                              sd_bp=4., fraction=1., effective_read_count=len(lengths))])
    return comps, reads, {alias: reads.read_id.tolist()}


class ValidationCorrectnessTests(unittest.TestCase):
    def test_margin_rejects_ties_and_keeps_clear_winners(self):
        for margin, expected in [(0., "Ambiguous"), (.01, "Ambiguous"), (.2, "ITD-supporting")]:
            with self.subTest(margin=margin):
                row = dict(ref_alias="ITD_1", metric_value=2., delta=margin)
                self.assertEqual(classify_read_support(row, "z_score"), expected)

    def test_probability_margin_overrides_large_z_delta(self):
        row = dict(ref_alias="ITD_1", metric_value=2., delta=1., prob_delta=.01)
        self.assertEqual(classify_read_support(row, "z_score"), "Ambiguous")

    def test_binary_clear_winner_is_still_supported(self):
        row = dict(ref_alias="ITD_1", metric_value=.99, delta=.98)
        self.assertEqual(classify_read_support(row, "prob_delta"), "ITD-supporting")

    def test_missing_wt_keeps_actual_insertion_size(self):
        row = fit_length(396).comps.iloc[0]
        self.assertFalse(row.is_wt)
        self.assertTrue(row.peak_alias.startswith("ITD"))
        self.assertAlmostEqual(row.putative_itd_size, 60.)

    def test_true_wt_is_preserved(self):
        row = fit_length(336).comps.iloc[0]
        self.assertTrue(row.is_wt)
        self.assertEqual(row.peak_alias, "WT")

    def test_no_wt_length_refinement_uses_reference_baseline(self):
        comps, reads, subsets = cluster_fixture("ITD_1", [348]*20 + [360]*20)
        result = refine_peak_substructure_once(
            comps, reads, subsets, min_reads_for_refinement=20,
            min_child_fraction=.25, max_subpeak_sd=5, wt_amplicon_length=336,
        )
        self.assertEqual(len(result.comps), 2)
        for got, expected in zip(sorted(result.comps.putative_itd_size), [12., 24.]):
            self.assertAlmostEqual(got, expected)

    def test_no_wt_sequence_refinement_uses_reference_baseline(self):
        comps, reads, subsets = cluster_fixture("ITD_1", [348,348,360,360])
        result = assignments_to_result(
            comps=comps, reads_df=reads, peak_subsets=subsets,
            assignments={"ITD_1": {"r0":"a","r1":"a","r2":"b","r3":"b"}},
            wt_amplicon_length=336, min_child_fraction=.25, min_child_reads=2,
        )
        self.assertEqual(sorted(result.comps.putative_itd_size), [12.,24.])

    def test_wt_split_produces_usable_variant_aliases_and_sizes(self):
        comps, reads, subsets = cluster_fixture("WT", [336,336,348,348])
        result = assignments_to_result(
            comps=comps, reads_df=reads, peak_subsets=subsets,
            assignments={"WT": {"r0":"a","r1":"a","r2":"b","r3":"b"}},
            wt_amplicon_length=336, min_child_fraction=.25, min_child_reads=2,
        )
        aliases = set(result.comps.peak_alias)
        self.assertIn("WT", aliases)
        variant = (aliases - {"WT"}).pop()
        self.assertTrue(variant.startswith("ITD_WT_H"))
        self.assertEqual(set(result.peak_subsets[variant]), {"r2","r3"})
        self.assertEqual(set(result.reads_df.gmm_peak_alias), aliases)
        sizes = result.comps.set_index("peak_alias").putative_itd_size
        self.assertAlmostEqual(sizes[variant], 12.)
        self.assertAlmostEqual(sizes["WT"], 0.)

    def test_deletion_crossing_window_is_counted(self):
        best = pd.DataFrame([dict(read_id="r", ref_alias="ITD_1", pct_identity=1.,
                                 aligned_blocks=(((0,90),(120,400)),((0,90),(90,370))))])
        cons = pd.DataFrame([dict(peak_alias="ITD_1", median_ins_pos_ref=120, consensus_len=30)])
        result = validate_itd_supporting_reads(best, cons, gap_window=15, max_gap_bp=5, min_pid=.9)
        self.assertFalse(result.iloc[0].valid_support)

    def test_actual_alignment_rejects_short_itd_but_keeps_exact_read(self):
        insertion = WT[120:180]
        full = WT[:180] + insertion + WT[180:]
        short = WT[:180] + insertion[12:] + WT[180:]
        refs = {"WT": dict(ref_seq_with_itd=WT), "ITD_1": dict(ref_seq_with_itd=full)}
        reads = pd.DataFrame([dict(read_id=k, read_seq=s, strand="+") for k,s in [("full",full),("short",short)]])
        cons = pd.DataFrame([dict(peak_alias="ITD_1", median_ins_pos_ref=180, consensus_len=60)])
        result = align_reads_multi_ref_parallel(reads, refs, cons, logging.getLogger(__name__), threads=1).set_index("read_id")
        self.assertTrue(result.loc["full","valid_support"])
        self.assertFalse(result.loc["short","valid_support"])
        self.assertEqual(result.loc["short","filter_reason"], "too_many_gaps(12)")

    def test_two_tied_winners_among_five_real_references_are_ambiguous(self):
        full = WT[:180] + WT[120:180] + WT[180:]
        refs = {"WT": dict(ref_seq_with_itd=WT),
                "ITD_1": dict(ref_seq_with_itd=full), "ITD_2": dict(ref_seq_with_itd=full),
                "ITD_3": dict(ref_seq_with_itd="A"*len(full)),
                "ITD_4": dict(ref_seq_with_itd="C"*len(full))}
        reads = pd.DataFrame([dict(read_id="r", read_seq=full, strand="+")])
        cons = pd.DataFrame([dict(peak_alias=a, median_ins_pos_ref=180, consensus_len=60) for a in refs if a!="WT"])
        row = align_reads_multi_ref_parallel(reads, refs, cons, logging.getLogger(__name__), threads=1).iloc[0]
        self.assertGreater(row.metric_value, 1.)
        self.assertAlmostEqual(row.prob_delta, 0.)
        self.assertEqual(row.support_call, "Ambiguous")


if __name__ == "__main__":
    unittest.main()
