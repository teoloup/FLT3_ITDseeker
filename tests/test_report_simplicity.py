import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from itdseeker.Write_output import generate_itd_html_report


class ReportSimplicityTests(unittest.TestCase):
    def render(self,folder,uncertain=False):
        note=('Observed junction sequence alternatives; distinct biological clones are not established. '
              if uncertain else 'Recurrent observed interior haplotype; GMM-independent candidate. ')
        note+='Experimental flank locator: 3/92 observations normalized. Adjacent duplication identity 100% is annotation only.'
        summary=pd.DataFrame([dict(ref_alias='ITD_OBS_A',allele_frequency=.0067,n_total_reads=20457,
                                   n_itd_reads=137,plus_reads=73,minus_reads=64)])
        refs={'ITD_OBS_A':dict(itd_seq='A'*45,itd_length=45,genomic_insertion_pos=28034116,chr='chr13')}
        consensus=pd.DataFrame([dict(peak_alias='ITD_OBS_A',consensus_seq='A'*45,discovery_note=note)])
        path=generate_itd_html_report('test','hg38',22115,22115,summary,refs,folder,folder,consensus,haplotype_method='dada2')
        return Path(path).read_text()

    def test_generic_discovery_mechanisms_are_not_shown(self):
        with tempfile.TemporaryDirectory() as folder:
            report=self.render(folder)
        for term in ['GMM-independent','Experimental flank locator','duplication identity','3/92',
                     'annotation only','Two supported sequence versions']:
            self.assertNotIn(term,report)
        self.assertIn('0.67%',report)
        self.assertIn('137',report)
        self.assertNotIn('<details open>',report)

    def test_only_plain_essential_sequence_uncertainty_is_shown(self):
        with tempfile.TemporaryDirectory() as folder:
            report=self.render(folder,True)
        self.assertEqual(report.count('Two supported sequence versions differ by one base; this does not establish separate clones.'),1)
        self.assertNotIn('Experimental flank locator',report)

    def test_redundant_net_plot_is_not_embedded_even_when_available(self):
        def image(path):
            if path.endswith('_itd_gmm_fit_plot.png'): return 'PRIMARY_READ_LENGTH_IMAGE'
            if path.endswith('_itd_size_distribution.png'): return 'REDUNDANT_NET_IMAGE'
            return ''
        with tempfile.TemporaryDirectory() as folder, mock.patch(
                'itdseeker.Write_output.img_to_base64',side_effect=image):
            report=self.render(folder)
        self.assertEqual(report.count('PRIMARY_READ_LENGTH_IMAGE'),1)
        self.assertNotIn('REDUNDANT_NET_IMAGE',report)
        self.assertNotIn('Observed net read-length change',report)

    def test_requested_clinical_details_and_original_summary_are_restored(self):
        with tempfile.TemporaryDirectory() as folder, mock.patch(
                'itdseeker.Write_output.img_to_base64',side_effect=lambda path: 'MSA_IMAGE' if path.endswith('_MSA_consensus.png') else ''):
            report=self.render(folder)
        for text in ['<span class="chip">45 bp</span>', 'in frame', 'chr13:28,034,116',
                     'Strand balance', '73 + / 64 -', '20,457 reads used for AF',
                     'ITDs reported', 'Haplotype method', 'dada2', 'Multiple sequence alignment (MSA)', 'MSA_IMAGE']:
            self.assertIn(text,report)

        for column in ['ITD','Size','Frame','Position','AF','Reads','Depth']:
            self.assertRegex(report,r'<th(?: class="num")?>'+column+'</th>')

    def test_prior_candidate_msa_is_not_used_for_final_call(self):
        with tempfile.TemporaryDirectory() as folder:
            Path(folder,'test_ITD_1_H2_MSA_consensus.png').write_bytes(b'candidate msa')
            report=self.render(folder)
        self.assertNotIn('Additional candidate MSA plots',report)
        self.assertNotIn('candidate msa',report)
