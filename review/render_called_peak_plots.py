"""Render plot previews from retained evidence, with no fitting or read alignment.

The retained GMM tables precede denoising, so previews explicitly retain that
original fitted model. Final markers come from retained competitive validation.
"""
import base64
import html
from pathlib import Path
import sys

import pandas as pd
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from itdseeker.GMM_peaks import plot_gmm_itds
from itdseeker.Helper_functions import plot_itd_size_distribution, calculate_allele_frequencies_and_strand_bias

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
    components=pd.read_csv(audit/'default_gmm_components.tsv',sep='\t')
    reads=pd.read_csv(audit/'default_gmm_assignments.tsv.gz',sep='\t')
    lengths=[len(r.seq) for r in SeqIO.parse(str(audit/'fastq/trimmed.fastq'),'fastq') if len(r.seq)>=330]
    wt=next(r for r in SeqIO.parse(str(source/f'{sample}_validation_refs.fasta'),'fasta') if r.id=='WT')
    directory=OUT/f'{sample}_af{af:g}'
    directory.mkdir(exist_ok=True)
    gmm=plot_gmm_itds(reads,components,title=f'{sample}: retained GMM fit before denoising',
        out_prefix=str(directory/sample),all_read_lengths=lengths,
        validated=calls,reported_aliases=set(calls),wt_reference_length=len(wt.seq))
    insertions=pd.read_csv(source/f'{sample}_itd_insertions.tsv',sep='\t')
    size=plot_itd_size_distribution(insertions,directory,sample,validated=calls,reported_aliases=set(calls))
    images=''.join(f'<img alt="{label}" src="data:image/png;base64,{base64.b64encode(Path(path).read_bytes()).decode()}">'
                   for label,path in [('Read lengths and final calls',gmm),('Insertion observations and final calls',size)])
    cards.append(f'<section><h2>{sample} — reporting AF ≥ {af*100:g}%</h2>{images}</section>')
    print(sample,af,[(a,c['itd_len'],c['reads']) for a,c in calls.items()],flush=True)

target=ROOT/'FRAGMENT_DATA/recheck_20261019/called_peak_plot_previews.html'
target.write_text('''<!doctype html><html><meta charset="utf-8"><title>Final FLT3 call plot previews</title>
<style>body{font:16px/1.5 system-ui;max-width:1250px;margin:30px auto;padding:0 20px;color:#203442}
img{width:100%;height:auto}section{border-top:1px solid #ccc;margin-top:30px}h1,h2{color:#183b56}</style>
<h1>Final FLT3 call plot previews</h1><p><a href="FLT3_fragment_recheck.html">Fragment recheck report</a></p>
<p>Histograms and curves retain their observed read/GMM meanings. The final-call panels show only
the calls surviving the displayed reporting threshold, one row per sequence. Full allele position
is WT reference length (336 bp) plus insertion size; insertion plots show inserted bases relative to WT (0 bp).
Markers are sequence calls, not fitted Gaussian peaks or extra histogram counts.</p>
<p>These previews reuse the retained original GMM fit before denoising and the final saved competitive
validation. No GMM fitting, denoising or read alignment was repeated. Production plots use the pipeline's current
GMM/refined candidate tables plus these same final-call panels. The two 45 bp sequences differ at one junction
base; separate biological clones are not established. The 0.1% preview includes both alternatives;
the 0.5% preview applies the original review reporting threshold. The caller default remains 5%.</p>
'''+''.join(cards)+'</html>',encoding='utf-8')
print(target,flush=True)
