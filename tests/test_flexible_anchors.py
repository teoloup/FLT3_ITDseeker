import random
import unittest

import pandas as pd

from itdseeker.observed_haplotypes import (
    anchor_coordinates, anchor_matches, locate_anchor_pair,
    discover_observed_haplotypes, add_observed_candidates,
)
from scripts.simulate_itd_data import DEFAULT_REF_WT as WT


def observations(seq):
    return {str(i):dict(seq=seq,qual='I'*len(seq),strand='+' if i%2 else '-') for i in range(24)}


def mutate(seq, positions):
    bases=list(seq)
    for pos in positions:
        bases[pos]='A' if bases[pos]!='A' else 'T'
    return ''.join(bases)


class FlexibleAnchorTests(unittest.TestCase):
    def test_shift_coordinates_and_invalid_settings(self):
        self.assertEqual(anchor_coordinates(336,20,10),(15,325))
        self.assertEqual(anchor_coordinates(336,20,15),(10,330))
        for length,shift in [(0,10),(20,-1),(20,22),(200,10)]:
            with self.assertRaises(ValueError):
                anchor_coordinates(336,length,shift)
        with self.assertRaises(ValueError):
            discover_observed_haplotypes({},WT,anchor_max_mismatches=20)

    def test_seed_partition_matches_exhaustive_substitution_search(self):
        rng=random.Random(3021)
        for mismatches in range(4):
            anchor=''.join(rng.choices('ACGT',k=20))
            changed=mutate(anchor,[0,6,13][:mismatches])
            seq=''.join(rng.choices('ACGT',k=12))+changed+''.join(rng.choices('ACGT',k=25))
            expected=[(i,sum(a!=b for a,b in zip(seq[i:i+20],anchor)))
                      for i in range(len(seq)-19)
                      if sum(a!=b for a,b in zip(seq[i:i+20],anchor))<=mismatches]
            self.assertEqual(anchor_matches(seq,anchor,mismatches),expected)
            self.assertIn((12,mismatches),expected)

    def test_three_substitutions_per_anchor_are_normalized_and_audited(self):
        allele=WT[:210]+WT[189:210]+WT[210:]
        changed=mutate(allele,[15,20,25,305+21,309+21,313+21])
        seeds,audit,stats=discover_observed_haplotypes(observations(changed),WT)
        self.assertEqual(len(seeds),1)
        seed=seeds[0]
        self.assertEqual(seed['allele'],allele)
        self.assertEqual(seed['normalized_flank_reads'],24)
        self.assertEqual(seed['exact_flank_reads'],0)
        self.assertEqual(seed['left_anchor_substitutions'],72)
        self.assertEqual(seed['right_anchor_substitutions'],72)
        self.assertEqual(stats['normalized_flank_reads'],24)
        consensus=add_observed_candidates(pd.DataFrame(),seeds,audit,WT)
        self.assertEqual(consensus.iloc[0].normalized_flank_reads,24)
        self.assertIn('not called flank variants',consensus.iloc[0].discovery_note)

    def test_interior_is_not_repaired_by_anchor_normalization(self):
        allele=WT[:210]+WT[189:210]+WT[210:]
        seeds,audit,_=discover_observed_haplotypes(observations(mutate(allele,[36])),WT)
        self.assertFalse(seeds)
        self.assertEqual(audit[0]['status'],'not_single_complete_insertion')

    def test_multiple_plausible_ordered_pairs_are_rejected_even_with_exact_hit(self):
        left='ACGTCAGTGACTGACCTAGC'
        right='TGCAATCGGCTAGTCACGAT'
        seq=left+'AAAA'+mutate(left,[7])+'CCCC'+right
        pair,status=locate_anchor_pair(seq,left,right,1)
        self.assertIsNone(pair)
        self.assertEqual(status,'ambiguous_anchor_pair')
        # Reversed order is not an acceptable flank pair.
        self.assertEqual(locate_anchor_pair(right+'CCCC'+left,left,right,0)[1],
                         'no_ordered_anchor_pair')

    def test_four_substitutions_and_anchor_indels_are_not_permitted(self):
        anchor='ACGTCAGTGACTGACCTAGC'
        self.assertFalse(anchor_matches(mutate(anchor,[0,5,10,15]),anchor,3))
        self.assertFalse(anchor_matches(anchor[:7]+anchor[8:],anchor,3))

    def test_low_duplication_identity_is_annotation_not_filter(self):
        allele=WT[:210]+'CGCGCCGCGCCGCGCCGCGCC'+WT[210:]
        seeds,audit,_=discover_observed_haplotypes(observations(allele),WT)
        self.assertEqual(len(seeds),1)
        self.assertLess(seeds[0]['duplication_identity'],.8)
        consensus=add_observed_candidates(pd.DataFrame(),seeds,audit,WT)
        self.assertEqual(consensus.iloc[0].duplication_identity,seeds[0]['duplication_identity'])
        self.assertIn('annotation only',consensus.iloc[0].discovery_note)
