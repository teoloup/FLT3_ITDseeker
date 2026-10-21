import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import logging
import numpy as np
import pandas as pd
from Helper_functions import softmax, classify_read_support, validate_itd_supporting_reads
from GMM_peaks import fit_gmm_itds
from haplotype_split import assignments_to_result

kw=dict(min_gmm_fraction=.01,max_itds_detected=6,min_ggmm_peak_distance=10,max_peak_sd=5,assign_width_factor=2,assign_mode='manual',prob_threshold=.85)
print('Z SCORE:')
for scores in ([.9,.9,.1], [.9,.9,.1,.1,.1], [.90001,.9,.1,.1,.1]):
    p=softmax(scores,500)
    z=(p-p.mean())/(p.std(ddof=0)+1e-8)
    call=classify_read_support({'ref_alias':'ITD_1','metric_value':z[0],'delta':z[0]-z[1]},'z_score')
    print('scores=',scores,'probabilities=',p.round(5),'zbest=',round(z[0],5),'delta=',round(z[0]-z[1],5),'call=',call)
print('\nABSENT WT:')
r=fit_gmm_itds({f'r{i}':'A'*396 for i in range(200)},force_k=1,**kw)
print(r.comps[['peak_alias','mean_bp','putative_itd_size','effective_read_count']].to_string(index=False))
print('\nLOW DEPTH:')
try:
    fit_gmm_itds({f'r{i}':'A'*(336+i) for i in range(3)},**kw)
except Exception as e:
    print(type(e).__name__,str(e))
print('\nMERGED POSTERIORS:')
rng=np.random.default_rng(123)
lengths=np.rint(np.concatenate([rng.normal(336,2,300),rng.normal(342,2,300)])).astype(int)
reads={f'r{i}':'A'*int(n) for i,n in enumerate(lengths)}
r=fit_gmm_itds(reads,force_k=2,**dict(kw,assign_mode='predict_proba'))
p=r.gmm.predict_proba(lengths.reshape(-1,1))
print('merged_peaks=',len(r.comps),'input_reads=',len(reads),'assigned_reads=',len(r.reads_df),'all merged posteriors=',np.unique(p.sum(axis=1).round(6)))
print('\nPARTIAL CLUSTER ASSIGNMENTS:')
rd=pd.DataFrame([dict(read_id=f'r{i}',read_len=396,read_seq='A'*396,strand='+',gmm_peak_alias='ITD_1',is_ambiguous=False) for i in range(100)])
cp=pd.DataFrame([dict(peak_alias='WT',mean_bp=336,sd_bp=1),dict(peak_alias='ITD_1',mean_bp=396,sd_bp=1)])
a={f'r{i}':('a' if i<40 else 'b') for i in range(80)}
out=assignments_to_result(comps=cp,reads_df=rd,peak_subsets={'WT':[],'ITD_1':rd.read_id.tolist()},assignments={'ITD_1':a},wt_amplicon_length=336,min_child_fraction=.15,min_child_reads=20)
print('subset_reads=',sum(map(len,out.peak_subsets.values())),'dataframe_reads=',len(out.reads_df),'stale aliases=',sorted(set(out.reads_df.gmm_peak_alias)-set(out.comps.peak_alias)))
print('\nWT-DERIVED ALIAS:')
row=dict(read_id='hidden',ref_alias='WT_H2',pct_identity=1.,aligned_blocks=np.array([[[0,350]],[[0,350]]]),metric_value=1.4,delta=1.)
print('classification=',classify_read_support(row,'z_score'))
v=validate_itd_supporting_reads(pd.DataFrame([row]),pd.DataFrame([dict(peak_alias='WT_H2',median_ins_pos_ref=150,consensus_len=14)]),gap_window=15,max_gap_bp=5,min_pid=.9)
print(v[['ref_alias','valid_support','filter_reason']].to_string(index=False))
