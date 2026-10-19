"""Render actual final-support profiles from retained validation and read lengths."""
import base64
import html
from pathlib import Path
import sys

import pandas as pd
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from itdseeker.Helper_functions import calculate_allele_frequencies_and_strand_bias
from itdseeker.final_call_plots import final_support_lengths, plot_final_call_distributions

OUT=ROOT/'review/observed_rescue_plots'
OUT.mkdir(exist_ok=True)
cards=[]
for sample,af in [('11528',.005),('13462',.005),('13462',.001)]:
    source=ROOT/'review/observed_rescue_benchmark'/sample/'enabled/results/flt3_data'
    audit=ROOT/'FRAGMENT_DATA/recheck_20261019'/sample
    validation=pd.read_csv(source/f'{sample}_validation_read_support.tsv',sep='\t')
    summary=calculate_allele_frequencies_and_strand_bias(validation,af)
    consensus=pd.read_csv(source/f'{sample}_itd_consensus_seq.tsv',sep='\t')
    sizes=dict(zip(consensus.peak_alias,consensus.consensus_len))
    calls={r.ref_alias:dict(itd_len=int(sizes[r.ref_alias]),af=float(r.allele_frequency),
                           reads=int(r.n_itd_reads)) for r in summary.itertuples()}
    lengths={r.id:len(r.seq) for r in SeqIO.parse(str(audit/'fastq/trimmed.fastq'),'fastq') if len(r.seq)>=330}
    wt=next(r for r in SeqIO.parse(str(source/f'{sample}_validation_refs.fasta'),'fasta') if r.id=='WT')
    directory=OUT/f'{sample}_af{af:g}'
    directory.mkdir(exist_ok=True)
    frame=final_support_lengths(validation,lengths,calls,len(wt.seq))
    plots=plot_final_call_distributions(frame,calls,directory,sample,len(wt.seq))
    images=''.join(f'<img alt="{label}" src="data:image/png;base64,{base64.b64encode(Path(path).read_bytes()).decode()}">'
                   for label,path in [('Final supporting read lengths',plots['read_length_plot']),
                                      ('Final supporting net read-length changes',plots['net_read_length_plot'])])
    cards.append(f'<section><h2>{sample} — reporting AF ≥ {af*100:g}%</h2>{images}</section>')
    print(sample,af,[(a,c['itd_len'],c['reads']) for a,c in calls.items()],flush=True)

target=ROOT/'FRAGMENT_DATA/recheck_20261019/called_peak_plot_previews.html'
target.write_text('''<!doctype html><html><meta charset="utf-8"><title>Final FLT3 call plot previews</title>
<style>body{font:16px/1.5 system-ui;max-width:1250px;margin:30px auto;padding:0 20px;color:#203442}
img{width:100%;height:auto}section{border-top:1px solid #ccc;margin-top:30px}h1,h2{color:#183b56}</style>
<h1>Final FLT3 call plot previews</h1><p><a href="FLT3_fragment_recheck.html">Fragment recheck report</a></p>
<p>Main histograms use actual trimmed lengths of reads that pass final competitive validation and support
the displayed reported allele. The full panel includes WT; the zoom contains reported ITD reads only.
The second plot shows actual trimmed read length minus the 336 bp WT reference, not extracted insertion
events. Dotted markers show exact called insertion sizes (or WT reference plus those sizes in the full plot).
Trimming and sequencing indels can shift observed histogram peaks away from exact called sizes.</p>
<p>These previews reuse saved final competitive validation and retained FASTQ lengths matched by read ID.
No GMM fitting, denoising or read alignment was repeated. No earlier GMM cluster assignments supply the main
histogram counts. The two 45 bp sequences differ at one junction
base; separate biological clones are not established. The 0.1% preview includes both alternatives;
the 0.5% preview applies the original review reporting threshold. The caller default remains 5%.</p>
'''+''.join(cards)+'</html>',encoding='utf-8')
print(target,flush=True)
