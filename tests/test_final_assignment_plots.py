import tempfile
import unittest
from unittest import mock

import pandas as pd

from itdseeker.final_call_plots import final_support_lengths, plot_final_call_distributions


class FinalAssignmentPlotTests(unittest.TestCase):
    def setUp(self):
        self.calls={'ITD_OBS_A':dict(itd_len=45,af=.2,reads=2),
                    'ITD_OBS_B':dict(itd_len=45,af=.1,reads=1)}
        self.validation=pd.DataFrame([
            dict(read_id='a',ref_alias='ITD_OBS_A',support_call='ITD-supporting',valid_support=True,filter_reason='pass'),
            dict(read_id='b',ref_alias='ITD_OBS_A',support_call='ITD-supporting',valid_support=True,filter_reason='pass'),
            dict(read_id='c',ref_alias='ITD_OBS_B',support_call='ITD-supporting',valid_support=True,filter_reason='pass'),
            dict(read_id='w',ref_alias='WT',support_call='WT-supporting',valid_support=True,filter_reason='pass'),
            dict(read_id='bad',ref_alias='ITD_OBS_A',support_call='ITD-supporting',valid_support=False,filter_reason='gaps'),
            dict(read_id='old',ref_alias='ITD_OLD_H2',support_call='ITD-supporting',valid_support=True,filter_reason='pass'),
        ])
        self.lengths=dict(a=377,b=380,c=379,w=334,bad=600,old=440)

    def test_only_final_passed_reported_reads_supply_actual_lengths(self):
        frame=final_support_lengths(self.validation,self.lengths,self.calls,336)
        self.assertEqual(set(frame.read_id),{'a','b','c','w'})
        a=frame.loc[frame.ref_alias.eq('ITD_OBS_A')]
        self.assertEqual(list(a.read_length),[377,380])
        self.assertEqual(list(a.net_read_length),[41,44])
        self.assertEqual(list(a.called_insertion_size),[45,45])

    def test_plot_histograms_use_measured_lengths_not_a_consensus_spike(self):
        frame=final_support_lengths(self.validation,self.lengths,self.calls,336)
        with tempfile.TemporaryDirectory() as folder, mock.patch('matplotlib.axes.Axes.hist') as hist, \
                mock.patch('matplotlib.axes.Axes.legend'):
            result=plot_final_call_distributions(frame,self.calls,folder,'test',336)
        supplied=[tuple(c.args[0]) for c in hist.call_args_list]
        self.assertIn((377,380),supplied)
        self.assertIn((41,44),supplied)
        self.assertIn((379,),supplied)
        self.assertNotIn((381,381),supplied)
        self.assertNotIn((45,45),supplied)
        self.assertEqual(result['profiles']['ITD_OBS_A']['net_read_length_counts'],{41:1,44:1})

    def test_missing_actual_length_or_wrong_final_count_is_rejected(self):
        missing=dict(self.lengths)
        del missing['a']
        with self.assertRaisesRegex(ValueError,'missing'):
            final_support_lengths(self.validation,missing,self.calls,336)
        calls={a:dict(c) for a,c in self.calls.items()}
        calls['ITD_OBS_A']['reads']=3
        with self.assertRaisesRegex(ValueError,'disagrees'):
            final_support_lengths(self.validation,self.lengths,calls,336)

    def test_reporting_filter_removes_lower_af_allele_distribution(self):
        calls={'ITD_OBS_A':self.calls['ITD_OBS_A']}
        frame=final_support_lengths(self.validation,self.lengths,calls,336)
        self.assertNotIn('ITD_OBS_B',set(frame.ref_alias))
        self.assertNotIn('ITD_OLD_H2',set(frame.ref_alias))
