import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
from Multiple_seq_aligment_toolkit import dedup_with_counts,select_unique_panel,run_muscle5_on_pairs,weighted_consensus_from_msa
from Helper_functions import build_default_aligner,find_insertions
from simulate_itd_data import DEFAULT_REF_WT as ref
p=Path('review/runs/dada2/fixed/sim_no_wt/flt3_data/sim_no_wt_itd_insertions.tsv')
df=pd.read_csv(p,sep='\t')
print('positions:',df.ins_pos_ref.value_counts().head(12).to_dict())
for mode in ['isolated','context']:
    seqs=df.ins_seq.tolist() if mode=='isolated' else [ref[:int(r.ins_pos_ref)]+r.ins_seq+ref[int(r.ins_pos_ref):] for r in df.itertuples()]
    panel=select_unique_panel(dedup_with_counts(seqs),150,.98)
    msa,w=run_muscle5_on_pairs(panel)
    cons,minor=weighted_consensus_from_msa(msa,w,base_threshold=.7,min_col_coverage=.7)
    if mode=='context':
        aln=build_default_aligner().align(ref,cons)[0]
        print(mode,'length',len(cons),'minor',minor,'insertions',[(p,n,cons[q:q+n]) for p,q,n in find_insertions(aln)],'exact_truth',cons==ref[:180]+ref[120:180]+ref[180:])
    else: print(mode,'length',len(cons),'minor',minor)
