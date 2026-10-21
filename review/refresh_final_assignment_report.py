"""Refresh main plots and the SAME saved report from final validation records.

No candidate discovery, validation alignment, AF change or VCF rewrite is performed.
Report-only MSAs use the retained final supporting reads.
Actual trimmed read lengths come from the retained FASTQ, matched by read ID.
"""
import argparse
import json
from pathlib import Path
import sys

import pandas as pd
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from itdseeker.final_support_msa import plot_final_support_msas
from itdseeker.final_call_plots import final_support_lengths, plot_final_call_distributions
from itdseeker.Helper_functions import calculate_allele_frequencies_and_strand_bias
from itdseeker.GMM_peaks import plot_gmm_itds
from itdseeker.Write_output import generate_itd_html_report

parser=argparse.ArgumentParser()
parser.add_argument('--sample',default='13462')
parser.add_argument('--mode',default='merge_review_af001')
parser.add_argument('--min-af',type=float,default=.001)
args=parser.parse_args()
sample=args.sample
results=ROOT/'review/observed_rescue_benchmark'/sample/args.mode/'results'
data=results/'flt3_data'
audit=ROOT/'FRAGMENT_DATA/recheck_20261019'/sample
validation=pd.read_csv(data/f'{sample}_validation_read_support.tsv',sep='\t')
consensus=pd.read_csv(data/f'{sample}_itd_consensus_seq.tsv',sep='\t')
sequences={r.id:str(r.seq) for r in SeqIO.parse(str(data/f'{sample}_validation_refs.fasta'),'fasta')}
wt_length=len(sequences['WT'])
read_lengths={r.id:len(r.seq) for r in SeqIO.parse(str(audit/'fastq/trimmed.fastq'),'fastq') if len(r.seq)>=330}
summary=calculate_allele_frequencies_and_strand_bias(validation,args.min_af)
sizes=dict(zip(consensus.peak_alias,consensus.consensus_len))
calls={r.ref_alias:dict(itd_len=int(sizes[r.ref_alias]),af=float(r.allele_frequency),reads=int(r.n_itd_reads))
       for r in summary.itertuples()}
frame=final_support_lengths(validation,read_lengths,calls,wt_length)
plots=plot_final_call_distributions(frame,calls,data,sample,wt_length)

components=pd.read_csv(audit/'default_gmm_components.tsv',sep='\t')
assignments=pd.read_csv(audit/'default_gmm_assignments.tsv.gz',sep='\t')
plot_gmm_itds(assignments,components,title=f'{sample}: retained earlier-stage GMM fit before denoising',
              out_prefix=str(data/f'{sample}_candidate_gmm_fit'),all_read_lengths=list(read_lengths.values()))
refs={}
for _,row in consensus.iterrows():
    alias=row['peak_alias']
    pos=int(row['consensus_ins_pos_ref'])
    refs[alias]=dict(alias=alias,itd_seq=row['consensus_seq'],ref_seq_with_itd=sequences[alias],
                    local_insertion_pos=pos,genomic_insertion_pos=28033881+pos-1,
                    chr='chr13',itd_length=int(row['consensus_len']),
                    discovery_note=row.get('discovery_note',''),
                    duplication_identity=row.get('duplication_identity'),
                    normalized_flank_reads=row.get('normalized_flank_reads'))
sequences_by_read={r.id:str(r.seq) for r in SeqIO.parse(str(audit/'fastq/trimmed.fastq'),'fastq')}
plot_final_support_msas(validation,sequences_by_read,summary,refs,data,sample)
generate_itd_html_report(sample,'hg38',len(read_lengths),len(validation),summary,refs,
                         str(results),str(data),consensus,haplotype_method='dada2')
record=dict(sample=sample,mode=args.mode,min_af=args.min_af,wt_reference_length=wt_length,
            report=str(results/f'{sample}_itd_report.html'),
            input_validation_rows=len(validation),matched_final_read_counts=frame.ref_alias.value_counts().to_dict(),
            source='Final passed competitive assignments joined by read ID to measured trimmed FASTQ lengths',
            net_axis='Actual trimmed read length minus WT reference length; not extracted insertion events',
            plots=plots)
target=ROOT/'review'/f'final_assignment_plot_check_{sample}.json'
target.write_text(json.dumps(record,indent=2))
print(record['report'])
print(record['matched_final_read_counts'])
print(target)
