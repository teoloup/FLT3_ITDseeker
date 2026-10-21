import tempfile
import unittest
from pathlib import Path
from unittest import mock

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from itdseeker.call_plots import reported_calls, draw_call_panel
from itdseeker.GMM_peaks import plot_gmm_itds
from itdseeker.Helper_functions import plot_itd_size_distribution


class CalledPeakPlotTests(unittest.TestCase):
    def test_report_filters_and_same_size_sequences_have_separate_rows(self):
        calls={'ITD_OBS_A':dict(itd_len=45,af=.0067,reads=137),
               'ITD_OBS_B':dict(itd_len=45,af=.0036,reads=74),
               'ITD_1':dict(itd_len=48,af=.0116,reads=238)}
        calls['WT']=dict(itd_len=0,af=.97,reads=20008)
        self.assertEqual([a for a,_ in reported_calls(calls,{'ITD_OBS_A','ITD_1'})],
                         ['ITD_OBS_A','ITD_1'])
        fig,ax=plt.subplots()
        positions=draw_call_panel(ax,reported_calls(calls),wt_length=336,full_allele=True)
        np.testing.assert_array_equal(positions,[381,381,384])
        self.assertEqual(list(ax.lines[1].get_ydata()),[0])
        self.assertEqual(list(ax.lines[2].get_ydata()),[1])
        self.assertIn('74 reads',ax.get_yticklabels()[1].get_text())
        plt.close(fig)

    def test_size_plot_includes_rescue_without_any_candidate_observations(self):
        calls={'ITD_OBS_A':dict(itd_len=21,af=.0259,reads=278)}
        with tempfile.TemporaryDirectory() as folder, mock.patch(
                'matplotlib.figure.Figure.savefig',autospec=True) as save:
            plot_itd_size_distribution(pd.DataFrame(),folder,'sample',validated=calls,
                                       reported_aliases=set(calls))
            fig=save.call_args.args[0]
            self.assertEqual(len(fig.axes),2)
            self.assertEqual(list(fig.axes[-1].lines[1].get_xdata()),[21])
            self.assertIn('ITD_OBS_A',fig.axes[-1].get_yticklabels()[0].get_text())
            self.assertFalse(fig.axes[0].patches)

    def test_gmm_rescue_without_retained_components_does_not_invent_fit(self):
        calls={'ITD_OBS_A':dict(itd_len=21,af=.0259,reads=278)}
        with tempfile.TemporaryDirectory() as folder, mock.patch(
                'matplotlib.figure.Figure.savefig',autospec=True) as save:
            plot_gmm_itds(pd.DataFrame(columns=['read_len']),pd.DataFrame(),
                title='No GMM components',out_prefix=str(Path(folder)/'sample'),
                all_read_lengths=np.array([334]*40+[357]*10),validated=calls,
                reported_aliases=set(calls),wt_reference_length=336)
            fig=save.call_args.args[0]
            self.assertFalse(fig.axes[0].lines)
            self.assertEqual(list(fig.axes[-1].lines[1].get_xdata()),[357])

    def test_insertion_coordinate_is_not_fitted_read_length(self):
        fig,ax=plt.subplots()
        positions=draw_call_panel(ax,[('ITD_OBS_A',dict(itd_len=21,af=.0259,reads=278))])
        np.testing.assert_array_equal(positions,[21])
        self.assertIn('WT = 0',ax.get_xlabel())
        plt.close(fig)
