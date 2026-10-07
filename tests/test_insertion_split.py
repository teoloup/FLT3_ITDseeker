"""Splitting a peak by insertion length, and keeping every read of a split peak."""
import random
import unittest

import pandas as pd

from itdseeker.Helper_functions import (
    extract_itd_insertions_from_subset_parallel, insertion_window_floor,
)
from itdseeker.haplotype_split import (
    assignments_to_result, insertion_length_modes, merge_duplicate_candidates,
    split_by_insertion_length,
)


def lengths(**counts):
    """lengths(L51=536, L57=172) -> [51]*536 + [57]*172"""
    return [int(k[1:]) for k, n in counts.items() for _ in range(n)]


def modes(values):
    return insertion_length_modes(values, min_reads=20, min_fraction=.15)


class InsertionLengthModeTests(unittest.TestCase):
    def test_two_separated_groups_are_two_itds(self):
        # the 12808 sample: 51 and 57 bp with a near-empty gap between
        values = lengths(L49=10, L50=12, L51=536, L52=17, L53=4, L54=5, L55=2, L56=7, L57=172)
        self.assertEqual(modes(values), [51, 57])

    def test_noise_shoulders_stay_one_itd(self):
        self.assertEqual(modes(lengths(L49=10, L50=40, L51=536, L52=60, L53=12)), [51])

    def test_minor_group_below_guardrails_is_ignored(self):
        self.assertEqual(modes(lengths(L51=536, L57=19)), [51])          # under 20 reads
        self.assertEqual(modes(lengths(L51=536, L57=60)), [51])          # under 15%

    def test_groups_closer_than_the_gap_are_not_split(self):
        self.assertEqual(modes(lengths(L51=300, L52=80, L53=200)), [51])

    def test_shoulder_without_a_dip_is_not_split(self):
        self.assertEqual(modes(lengths(L51=300, L52=150, L53=100, L54=110)), [51])


def peak_fixture():
    """ITD_1 holds 30 reads with a 51 bp insertion, 12 with 57 bp, 5 with none."""
    groups = [("a", 30, 387, 51), ("b", 12, 393, 57), ("n", 5, 386, None)]
    reads = pd.DataFrame([
        dict(read_id=f"{tag}{i}", read_len=read_len, read_seq="A" * read_len,
             gmm_peak_alias="ITD_1", is_ambiguous=False, strand="+")
        for tag, n, read_len, _ in groups for i in range(n)
    ])
    insertions = pd.DataFrame([
        dict(peak_alias="ITD_1", read_id=f"{tag}{i}", ins_pos_ref=200, ins_len=ins_len,
             ins_seq="C" * ins_len)
        for tag, n, _, ins_len in groups if ins_len for i in range(n)
    ])
    comps = pd.DataFrame([dict(peak_alias="ITD_1", mean_bp=388.5, sd_bp=3.,
                               fraction=1., effective_read_count=len(reads))])
    return comps, reads, {"ITD_1": reads.read_id.tolist()}, insertions


class SplitByInsertionLengthTests(unittest.TestCase):
    def test_peak_splits_and_every_read_is_kept(self):
        comps, reads, subsets, insertions = peak_fixture()
        result, relabelled = split_by_insertion_length(
            comps=comps, reads_df=reads, peak_subsets=subsets, insertions_df=insertions,
            wt_amplicon_length=336, min_child_fraction=.15, min_child_reads=10,
        )
        self.assertEqual(sorted(result.peak_subsets), ["ITD_1_H1", "ITD_1_H2"])
        kept = [r for ids in result.peak_subsets.values() for r in ids]
        self.assertEqual(sorted(kept), sorted(reads.read_id))
        major, minor = (set(result.peak_subsets[a]) for a in ("ITD_1_H1", "ITD_1_H2"))
        self.assertEqual(minor, {f"b{i}" for i in range(12)})
        self.assertTrue({f"n{i}" for i in range(5)} <= major)
        alias_of = dict(zip(result.reads_df.read_id, result.reads_df.gmm_peak_alias))
        self.assertTrue(all(alias_of[r] == a for r, a in zip(relabelled.read_id, relabelled.peak_alias)))
        sizes = result.comps.set_index("peak_alias").putative_itd_size
        self.assertAlmostEqual(sizes["ITD_1_H2"], 57.)

    def test_single_length_leaves_the_peak_alone(self):
        comps, reads, subsets, insertions = peak_fixture()
        insertions = insertions[insertions.ins_len == 51]
        result, relabelled = split_by_insertion_length(
            comps=comps, reads_df=reads, peak_subsets=subsets, insertions_df=insertions,
            wt_amplicon_length=336, min_child_fraction=.15, min_child_reads=10,
        )
        self.assertEqual(list(result.peak_subsets), ["ITD_1"])
        self.assertIs(relabelled, insertions)


class ExtractionWindowTests(unittest.TestCase):
    def test_merge_distance_floor_keeps_a_second_itd_in_the_peak(self):
        # a 60 bp ITD sitting in a peak whose GMM size is 51 +/- 2 bp
        ref = "".join(random.Random(0).choices("ACGT", k=336))
        seq = ref[:200] + ref[140:200] + ref[200:]
        reads = pd.DataFrame([dict(read_id=f"r{i}", read_seq=seq, strand="+") for i in range(3)])
        comps = pd.DataFrame([dict(peak_alias="ITD_1", putative_itd_size=51., sd_bp=2.)])
        kw = dict(reads_df=reads, read_ids_subset=list(reads.read_id), ref_seq=ref,
                  peak_alias="ITD_1", comps=comps, threads=1, itd_sd_factor=1.5,
                  min_itd_size=12, max_itd_size=300)
        self.assertTrue(extract_itd_insertions_from_subset_parallel(**kw).empty)
        found = extract_itd_insertions_from_subset_parallel(**kw, min_half_window=10)
        self.assertEqual(set(found.ins_len), {60})


class WindowFloorTests(unittest.TestCase):
    def comps(self, **sizes):
        rows = [dict(peak_alias="WT", putative_itd_size=0.)]
        return pd.DataFrame(rows + [dict(peak_alias=a, putative_itd_size=s) for a, s in sizes.items()])

    def test_lone_peak_gets_the_merge_distance(self):
        self.assertEqual(insertion_window_floor("ITD_1", self.comps(ITD_1=51.), 10), 10.)

    def test_floor_stops_halfway_to_the_next_itd_peak(self):
        comps = self.comps(ITD_1=51., ITD_2=62.)
        self.assertEqual(insertion_window_floor("ITD_1", comps, 10), 5.5)
        self.assertEqual(insertion_window_floor("ITD_2", comps, 10), 5.5)

    def test_same_size_haplotypes_get_no_floor(self):
        self.assertEqual(insertion_window_floor("ITD_1_H1", self.comps(ITD_1_H1=30., ITD_1_H2=30.), 10), 0.)


class MergeDuplicateTests(unittest.TestCase):
    REF = "".join(random.Random(1).choices("ACGT", k=336))

    def fixture(self, second_pos, second_seq):
        reads = pd.DataFrame([dict(read_id=f"{a}{i}", read_len=396, gmm_peak_alias=f"ITD_{a}")
                              for a, n in (("A", 30), ("B", 10)) for i in range(n)])
        subsets = {f"ITD_{a}": reads.read_id[reads.gmm_peak_alias == f"ITD_{a}"].tolist() for a in "AB"}
        comps = pd.DataFrame([
            dict(peak_alias="WT", mean_bp=336., sd_bp=2., fraction=.9, effective_allele_freq=.9),
            dict(peak_alias="ITD_A", mean_bp=396., sd_bp=2., fraction=.075, effective_allele_freq=.075),
            dict(peak_alias="ITD_B", mean_bp=396., sd_bp=2., fraction=.025, effective_allele_freq=.025),
        ])
        insertions = pd.DataFrame([dict(peak_alias=a, read_id=r, ins_len=60)
                                   for a, ids in subsets.items() for r in ids])
        cons = pd.DataFrame([
            dict(peak_alias="ITD_A", consensus_seq=self.REF[140:200], consensus_ins_pos_ref=200),
            dict(peak_alias="ITD_B", consensus_seq=second_seq, consensus_ins_pos_ref=second_pos),
        ])
        return merge_duplicate_candidates(comps=comps, reads_df=reads, peak_subsets=subsets,
                                          insertions_df=insertions, df_cons=cons,
                                          ref_seq=self.REF, wt_amplicon_length=336)

    def test_same_allele_written_at_another_position_is_merged(self):
        # the same 60 bp duplication, placed before its source copy instead of after
        result, insertions, cons = self.fixture(140, self.REF[140:200])
        self.assertEqual(list(cons.peak_alias), ["ITD_A"])
        self.assertEqual(sorted(result.peak_subsets), ["ITD_A"])
        self.assertEqual(len(result.peak_subsets["ITD_A"]), 40)
        self.assertEqual(set(result.reads_df.gmm_peak_alias), {"ITD_A"})
        self.assertEqual(set(insertions.peak_alias), {"ITD_A"})
        merged = result.comps.set_index("peak_alias").loc["ITD_A"]
        self.assertAlmostEqual(merged.effective_allele_freq, .1)

    def test_different_alleles_stay_apart(self):
        result, _, cons = self.fixture(200, self.REF[141:201])
        self.assertEqual(sorted(cons.peak_alias), ["ITD_A", "ITD_B"])
        self.assertEqual(sorted(result.peak_subsets), ["ITD_A", "ITD_B"])


class UnassignedReadTests(unittest.TestCase):
    def test_reads_without_a_cluster_fold_into_the_largest_child(self):
        comps, reads, subsets, _ = peak_fixture()
        labels = {r: ("b" if r.startswith("b") else "a") for r in reads.read_id if not r.startswith("n")}
        result = assignments_to_result(
            comps=comps, reads_df=reads, peak_subsets=subsets, assignments={"ITD_1": labels},
            wt_amplicon_length=336, min_child_fraction=.15, min_child_reads=10,
        )
        kept = [r for ids in result.peak_subsets.values() for r in ids]
        self.assertEqual(sorted(kept), sorted(reads.read_id))
        self.assertTrue({f"n{i}" for i in range(5)} <= set(result.peak_subsets["ITD_1_H1"]))


if __name__ == "__main__":
    unittest.main()
