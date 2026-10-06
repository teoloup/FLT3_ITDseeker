"""Context-consensus regression tests; requires the pipeline's MUSCLE dependency."""
import os
os.environ["MPLBACKEND"] = "Agg"
import importlib.util
import tempfile
import unittest
from pathlib import Path
import pandas as pd

from unittest import mock

AVAILABLE = importlib.util.find_spec('pymuscle5') is not None
if AVAILABLE:
    from itdseeker import Multiple_seq_aligment_toolkit as toolkit
    from itdseeker.Multiple_seq_aligment_toolkit import build_itd_consensus_sequences, boundary_interval, context_window
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

    def test_boundary_interval_spans_every_equivalent_position(self):
        truth=WT[:180]+WT[120:180]+WT[180:]
        a,b=boundary_interval(truth,WT)
        self.assertLessEqual(a,120)
        self.assertGreaterEqual(b,180)
        for p in (a,b):
            self.assertEqual(WT[:p]+truth[p:p+60]+WT[p:],truth)

    def test_window_starts_at_the_shared_position_and_widens_for_a_second_site(self):
        self.assertEqual(context_window([(120,180),(150,200)],[10,10],len(WT),flank=30),(120,len(WT)-180))
        self.assertEqual(context_window([(10,20),(len(WT)-5,len(WT)-5)],[1,1],len(WT),flank=30),(0,0))

    def test_muscle_sees_only_the_window_and_result_is_unchanged(self):
        rows=[dict(peak_alias='ITD_1',ins_pos_ref=180,ins_seq=WT[120:180])]*20
        real=toolkit.run_muscle5_on_pairs
        with mock.patch.object(toolkit,'run_muscle5_on_pairs',side_effect=real) as spy:
            result=self.consensus(rows).iloc[0]
        lengths={len(s) for s,_ in spy.call_args.args[0]}
        self.assertEqual(lengths,{60+2*toolkit.CONTEXT_FLANK_BP})
        p=int(result.consensus_ins_pos_ref)
        self.assertEqual(WT[:p]+result.consensus_seq+WT[p:],WT[:180]+WT[120:180]+WT[180:])
        self.assertEqual(result.allele_consensus_len,len(WT)+60)

    def test_stray_boundary_is_left_out_of_the_window(self):
        rows=[dict(peak_alias='ITD_1',ins_pos_ref=180,ins_seq=WT[120:180])]*40
        rows+=[dict(peak_alias='ITD_1',ins_pos_ref=40,ins_seq=WT[120:180])]
        real=toolkit.run_muscle5_on_pairs
        with mock.patch.object(toolkit,'run_muscle5_on_pairs',side_effect=real) as spy:
            result=self.consensus(rows).iloc[0]
        lengths={len(s) for s,_ in spy.call_args.args[0]}
        self.assertEqual(lengths,{60+2*toolkit.CONTEXT_FLANK_BP})
        p=int(result.consensus_ins_pos_ref)
        self.assertEqual(WT[:p]+result.consensus_seq+WT[p:],WT[:180]+WT[120:180]+WT[180:])

    def test_plot_is_written_whether_or_not_the_insertion_resolves(self):
        for payload,name in [(WT[120:180],'resolved'),('N'*60,'unresolved')]:
            with tempfile.TemporaryDirectory() as directory:
                build_itd_consensus_sequences(
                    pd.DataFrame([dict(peak_alias='ITD_1',ins_pos_ref=180,ins_seq=payload)]*20),
                    pd.DataFrame([dict(peak_alias='ITD_1',putative_itd_size=60,sd_bp=2)]),
                    ref_seq=WT,out_dir=directory,sample_name=name,min_col_coverage=.7,
                )
                self.assertTrue((Path(directory)/f'{name}_ITD_1_MSA_consensus.png').exists(),name)

if __name__=='__main__': unittest.main()
