"""Splitting a peak by insertion length, and keeping every read of a split peak."""
import unittest

import pandas as pd

from itdseeker.haplotype_split import (
    assignments_to_result, insertion_length_modes, split_by_insertion_length,
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
