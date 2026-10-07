"""Regressions for narrow integer-length peaks merged before manual assignment."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from itdseeker.GMM_peaks import _root_alias, fit_gmm_itds, plot_gmm_itds


class MergedPeakTests(unittest.TestCase):
    def test_merged_integer_peaks_retain_reads(self):
        # An average of near-zero child SDs excludes every integer read length
        # around the fractional merged mean. Total variance must include means.
        lengths = [330, 331, 332, 333, 336]
        reads = {f'{n}_{i}': 'A' * n for n in lengths for i in range(40)}
        result = fit_gmm_itds(
            reads, min_gmm_fraction=.01, max_itds_detected=5,
            min_ggmm_peak_distance=10, max_peak_sd=5,
            assign_width_factor=2, assign_mode='manual', prob_threshold=.85,
            force_k=5,
        )
        self.assertEqual(len(result.comps), 1)
        self.assertAlmostEqual(result.comps.iloc[0].sd_bp,
                               np.sqrt(np.var(lengths) + .001), places=5)
        self.assertEqual(len(result.reads_df), len(reads))

    def test_empty_plot_does_not_reduce_empty_array(self):
        with tempfile.TemporaryDirectory() as folder:
            plot_gmm_itds(pd.DataFrame(columns=['read_len']), pd.DataFrame(),
                          title='No assigned reads', out_prefix=str(Path(folder) / 'empty'))
            self.assertEqual(list(Path(folder).iterdir()), [])

    def test_split_peak_plot_with_zoom_and_reported_afs(self):
        reads = pd.DataFrame(dict(read_len=[336] * 200 + [387] * 40 + [393] * 12))
        comps = pd.DataFrame([
            dict(peak_alias='WT', mean_bp=336., sd_bp=2., fraction=.79,
                 effective_allele_freq=.79, is_wt=True),
            dict(peak_alias='ITD_1_H1', mean_bp=387., sd_bp=2., fraction=.16,
                 effective_allele_freq=.16, is_refined_child=True, split_by='insertion length'),
            dict(peak_alias='ITD_1_H2', mean_bp=393., sd_bp=2., fraction=.05,
                 effective_allele_freq=.05, is_refined_child=True, split_by='insertion length'),
        ])
        with tempfile.TemporaryDirectory() as folder:
            plot_gmm_itds(reads, comps, title='split', out_prefix=str(Path(folder) / 'p'),
                          reported_aliases={'ITD_1_H1'},
                          validated={'ITD_1_H1': dict(af=.049, reads=600, itd_len=51)},
                          all_read_lengths=[336] * 210 + [387] * 45 + [393] * 12 + [370] * 5)
            self.assertTrue((Path(folder) / 'p_plot.png').exists())

    def test_split_children_trace_back_to_their_gmm_peak(self):
        for alias, root in [('ITD_1_H2_H1', 'ITD_1'), ('ITD_1_H1_2', 'ITD_1'),
                            ('ITD_WT_H2', 'ITD_WT'), ('ITD_2', 'ITD_2')]:
            self.assertEqual(_root_alias(alias), root)


if __name__ == '__main__':
    unittest.main()
