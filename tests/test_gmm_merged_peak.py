"""Regressions for narrow integer-length peaks merged before manual assignment."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from itdseeker.GMM_peaks import fit_gmm_itds, plot_gmm_itds


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


if __name__ == '__main__':
    unittest.main()
