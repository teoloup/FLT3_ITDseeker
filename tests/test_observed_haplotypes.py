import unittest
import tempfile
from pathlib import Path
import pandas as pd
from itdseeker.observed_haplotypes import discover_observed_haplotypes, add_observed_candidates
from itdseeker.Pairwise_aligment_toolkit import _align_reads_chunk
from scripts.simulate_itd_data import DEFAULT_REF_WT as WT
from itdseeker.Write_output import generate_itd_html_report, export_itd_vcf


def reads(allele, n=24):
    return {str(i):dict(seq=allele,qual='I'*len(allele),strand='+' if i%2 else '-')
            for i in range(n)}


class ObservedDiscoveryTests(unittest.TestCase):
    def test_no_gmm_needed_and_duplicate_candidate_is_not_added(self):
        allele = WT[:210] + WT[189:210] + WT[210:]
        candidates,audit,stats = discover_observed_haplotypes(reads(allele),WT)
        self.assertEqual(len(candidates),1)
        self.assertEqual(stats['tested_seeds'],1)
        cons = add_observed_candidates(pd.DataFrame(),candidates,audit,WT)
        self.assertEqual(cons.iloc[0].consensus_len,21)
        again = add_observed_candidates(cons,candidates,audit,WT)
        self.assertEqual(len(again),1)
        self.assertEqual(audit[0]['status'],'existing_candidate')

    def test_short_isolated_gap_inside_longer_event_is_not_a_short_call(self):
        allele = WT[:210] + WT[195:210] + WT[210:]
        allele = allele[:240] + 'GAGA' + allele[240:]
        candidates,audit,_ = discover_observed_haplotypes(reads(allele),WT,min_size=3)
        self.assertFalse(candidates)
        self.assertEqual(audit[0]['status'],'not_single_complete_insertion')

    def test_weak_and_one_strand_observations_do_not_seed(self):
        allele = WT[:210] + WT[189:210] + WT[210:]
        self.assertFalse(discover_observed_haplotypes(reads(allele,18),WT)[0])
        data=reads(allele)
        for read in data.values(): read['strand']='+'
        self.assertFalse(discover_observed_haplotypes(data,WT)[0])

    def test_reconstructed_context_must_equal_observation(self):
        allele = WT[:210] + WT[189:210] + WT[210:]
        allele = allele[:250] + ('A' if allele[250]!='A' else 'T') + allele[251:]
        self.assertFalse(discover_observed_haplotypes(reads(allele),WT)[0])

    def test_ambiguous_observations_do_not_become_resolved_seeds(self):
        allele=WT[:210]+WT[189:199]+'N'+WT[200:210]+WT[210:]
        candidates,audit,_=discover_observed_haplotypes(reads(allele),WT)
        self.assertFalse(candidates)
        self.assertEqual(audit[0]['status'],'unresolved_bases')

    def test_candidate_cap_is_audited(self):
        a=WT[:210]+WT[189:210]+WT[210:]
        candidates,audit,_=discover_observed_haplotypes(reads(a),WT)
        self.assertTrue(add_observed_candidates(pd.DataFrame(),candidates,audit,WT,0).empty)
        self.assertEqual(audit[0]['status'],'candidate_cap')

    def test_junction_alternatives_are_retained_without_clone_claim(self):
        # Both complete observations are exact pure insertions with one base
        # difference inside the inserted copy, so neither is replaced by an N.
        a=WT[:210]+WT[189:210]+WT[210:]
        b=a[:220]+('A' if a[220]!='A' else 'T')+a[221:]
        data=reads(a)
        data.update({'b'+k:v for k,v in reads(b).items()})
        candidates,audit,_=discover_observed_haplotypes(data,WT)
        cons=add_observed_candidates(pd.DataFrame(),candidates,audit,WT)
        self.assertEqual(len(cons),2)
        self.assertTrue(all('not established' in n for n in cons.discovery_note))

    def test_cached_alignment_preserves_each_read_and_strand(self):
        refs={'WT':dict(ref_seq_with_itd=WT)}
        result=_align_reads_chunk([('a',WT,'+'),('b',WT,'-')],refs)
        self.assertEqual([(r['read_id'],r['strand']) for r in result],[('a','+'),('b','-')])
        for name in ['score','adjusted_score','pct_identity','span_identity']:
            self.assertEqual(result[0][name],result[1][name])

    def test_report_explains_observed_sequence_uncertainty(self):
        with tempfile.TemporaryDirectory() as directory:
            generate_itd_html_report('test','hg38',100,100,
                pd.DataFrame([dict(ref_alias='ITD_OBS_1',allele_frequency=.1,
                    n_total_reads=100,n_itd_reads=10)]),
                {'ITD_OBS_1':dict(itd_seq='ACGT'*3,itd_length=12,
                    genomic_insertion_pos=28034000,chr='chr13')},directory,directory,
                pd.DataFrame([dict(peak_alias='ITD_OBS_1',consensus_seq='ACGT'*3,
                    discovery_note='Distinct biological clones are not established.')]))
            self.assertIn('clones are not established',
                          (Path(directory)/'test_itd_report.html').read_text())

    def test_vcf_preserves_junction_uncertainty(self):
        note='Observed junction sequence alternatives; distinct biological clones are not established.'
        with tempfile.TemporaryDirectory() as directory:
            target=Path(directory)/'test.vcf'
            export_itd_vcf(pd.DataFrame([dict(ref_alias='ITD_OBS_1',allele_frequency=.1,
                n_total_reads=100,n_itd_reads=10)]),
                {'WT':dict(ref_seq_with_itd=WT), 'ITD_OBS_1':dict(chr='chr13',itd_length=12)},
                {'ITD_OBS_1':dict(itd_seq='ACGT'*3,local_insertion_pos=210,
                    genomic_insertion_pos=28034090,discovery_note=note)},
                pd.DataFrame(columns=['peak_alias']),str(target),genome_build='hg38')
            line=target.read_text().splitlines()[-1]
            self.assertIn(';SEQUENCE_UNCERTAINTY',line)
