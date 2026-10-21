import sys, re, json, logging
from pathlib import Path
from collections import Counter
import pandas as pd
sys.path.insert(0, '/workspace')
from itdseeker.Helper_functions import process_chunk
from itdseeker.Multiple_seq_aligment_toolkit import build_itd_consensus_sequences

def main():
    logging.basicConfig(level=logging.INFO)
    p=Path('/workspace/review/fastq_51_57_analysis')
    g=pd.read_csv(p/'gmm_read_assignments.tsv',sep='\t').fillna('')
    comps=pd.read_csv(p/'gmm_components.tsv',sep='\t')
    ref=re.search('DEFAULT_REF_WT = """(.*?)"""',Path('/workspace/Nano_ITDseeker.py').read_text()).group(1)
    row=comps[comps.peak_alias=='ITD_1'].iloc[0]
    sub=g[g.gmm_peak_alias=='ITD_1']
    ins=pd.DataFrame(process_chunk(list(sub[['read_id','read_seq','strand']].itertuples(index=False,name=None)),row.putative_itd_size-1.5*row.sd_bp,row.putative_itd_size+1.5*row.sd_bp,ref,'ITD_1'))
    ins.to_csv(p/'pipeline_insertions.tsv',sep='\t',index=False)
    print('PIPELINE_INSERTIONS',Counter(ins.ins_len),flush=True)
    cons=build_itd_consensus_sequences(all_itd_insertions=ins,comps=comps,ref_seq=ref,min_itd_size=12,max_itd_size=300,out_dir=str(p),sample_name='diagnostic',max_unique=150,min_weight_coverage=.98,base_threshold=.7,min_col_coverage=.7,ambiguous='N',threads=4)
    cons.to_csv(p/'pipeline_consensus.tsv',sep='\t',index=False)
    print('CONSENSUS',cons.to_string(index=False),flush=True)
    independent=pd.read_csv(p/'long_read_insertions.tsv',sep='\t')
    c=pd.read_csv(p/'haplotypes/dada2_ITD_1/clusters.tsv',sep='\t',header=None,names=['asv','read_id'])
    ids=set(g.read_id); mapped=dict(zip(c.read_id,c.asv)); extracted=dict(zip(ins.read_id,ins.ins_len))
    summary={}
    for size in (51,57):
        selected=independent[independent.insertions.map(lambda s:any(e[1]==size for e in json.loads(s)))]
        summary[str(size)]={'trimmed_support':len(selected),'strand_counts':dict(Counter(selected.strand)),'retained_in_GMM':int(sum(selected.read_id.isin(ids))),'DADA2_assignments':dict(Counter(mapped.get(r,'not_in_DADA2') for r in selected.read_id)),'pipeline_extracted_lengths':dict(Counter(str(extracted.get(r,'not_extracted')) for r in selected.read_id))}
    (p/'trace_summary.json').write_text(json.dumps(summary,indent=2))
    print('TRACE',json.dumps(summary,indent=2),flush=True)

if __name__=='__main__': main()
