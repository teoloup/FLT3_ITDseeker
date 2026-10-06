"""Context-consensus regression tests; requires the pipeline's MUSCLE dependency."""
import os
os.environ["MPLBACKEND"] = "Agg"
import importlib.util
import tempfile
import unittest
from pathlib import Path
import pandas as pd

AVAILABLE = importlib.util.find_spec('pymuscle5') is not None
if AVAILABLE:
    from itdseeker.Multiple_seq_aligment_toolkit import build_itd_consensus_sequences
    from scripts.simulate_itd_data import DEFAULT_REF_WT as WT

@unittest.skipUnless(AVAILABLE, 'pymuscle5 is required')
class ContextConsensusTests(unittest.TestCase):
    def consensus(self, rows, **kwargs):
        with tempfile.TemporaryDirectory() as directory:
            return build_itd_consensus_sequences(
                pd.DataFrame(rows),pd.DataFrame([dict(peak_alias='ITD_1',putative_itd_size=60,sd_bp=2)]),
                ref_seq=WT,out_dir=directory,sample_name='test',min_col_coverage=.7,
                **kwargs,
            )

    def test_equivalent_rotations_recover_exact_allele_and_anchor(self):
        truth=WT[:180]+WT[120:180]+WT[180:]
        rows=[]
        for p in [119,120,135,150,170,180]:
            payload=truth[p:p+60]
            self.assertEqual(WT[:p]+payload+WT[p:],truth)
            rows.extend([dict(peak_alias='ITD_1',ins_pos_ref=p,ins_seq=payload)]*4)
        result=self.consensus(rows).iloc[0]
        p=int(result.consensus_ins_pos_ref)
        self.assertEqual(result.consensus_len,60)
        self.assertEqual(WT[:p]+result.consensus_seq+WT[p:],truth)
        self.assertEqual(result.median_ins_pos_ref,p)
        self.assertNotEqual(result.raw_median_ins_pos_ref,p)
        self.assertEqual(result.n_unique,1)

    def test_payload_is_not_forced_to_match_a_tandem_repeat(self):
        payload='AGTCGATCGCTAGCATGTCAGATCCGTAAC'
        rows=[dict(peak_alias='ITD_1',ins_pos_ref=160,ins_seq=payload)]*20
        result=self.consensus(rows).iloc[0]
        p=int(result.consensus_ins_pos_ref)
        self.assertEqual(WT[:p]+result.consensus_seq+WT[p:],WT[:160]+payload+WT[160:])

    def test_consensus_size_limits_are_enforced(self):
        rows=[dict(peak_alias='ITD_1',ins_pos_ref=180,ins_seq=WT[120:180])]*20
        self.assertTrue(self.consensus(rows,max_itd_size=50).empty)

    def test_unresolved_bases_cannot_become_a_variant(self):
        rows=[dict(peak_alias='ITD_1',ins_pos_ref=180,ins_seq='N'*60)]*20
        self.assertTrue(self.consensus(rows).empty)

    def test_invalid_boundary_is_rejected(self):
        rows=[dict(peak_alias='ITD_1',ins_pos_ref=999,ins_seq='ACGT'*15)]
        with self.assertRaises(ValueError):
            self.consensus(rows)

if __name__=='__main__': unittest.main()
