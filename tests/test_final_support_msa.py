import tempfile
import unittest
from pathlib import Path
import pandas as pd
from itdseeker.final_support_msa import plot_final_support_msas


class FinalSupportMSATests(unittest.TestCase):
    def test_each_allele_has_only_its_own_passed_members_and_own_asset(self):
        validation=pd.DataFrame([
            dict(read_id='a1',ref_alias='A',valid_support=True,filter_reason='pass'),
            dict(read_id='a2',ref_alias='A',valid_support=True,filter_reason='pass'),
            dict(read_id='a3',ref_alias='A',valid_support=True,filter_reason='pass'),
            dict(read_id='b1',ref_alias='B',valid_support=True,filter_reason='pass'),
            dict(read_id='bad',ref_alias='A',valid_support=False,filter_reason='low_score'),
            dict(read_id='wt',ref_alias='WT',valid_support=True,filter_reason='pass')])
        validation['support_call']='ITD-supporting'
        sequences={'a1':'ACGTAAAACGTC','a2':'ACGTAAAACGTC','a3':'ACGTATAACGTC',
                   'b1':'ACGTCCCACGTC','bad':'TTTT','wt':'ACGTACGTC'}
        summary=pd.DataFrame([dict(ref_alias='A',n_itd_reads=3),dict(ref_alias='B',n_itd_reads=1)])
        references={alias:dict(ref_seq_with_itd=seq,local_insertion_pos=4,itd_length=3)
                    for alias,seq in [('A',sequences['a1']),('B',sequences['b1'])]}
        with tempfile.TemporaryDirectory() as folder:
            result=plot_final_support_msas(validation,sequences,summary,references,folder,'sample',max_unique=1)
            self.assertEqual(result['A']['read_ids'],['a1','a2','a3'])
            self.assertEqual(result['A']['panel_read_ids'],['a1','a2'])
            self.assertEqual(result['A']['represented_reads'],2)
            self.assertEqual(result['B']['read_ids'],['b1'])
            self.assertNotEqual(result['A']['plot'],result['B']['plot'])
            for info in result.values():
                self.assertTrue(Path(info['plot']).is_file())
