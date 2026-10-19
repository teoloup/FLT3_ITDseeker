import sys,re,logging
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
from Pairwise_aligment_toolkit import align_reads_multi_ref_parallel

def main():
    wt=re.search('DEFAULT_REF_WT = '+chr(34)*3+'(.*?)'+chr(34)*3,Path('Nano_ITDseeker.py').read_text()).group(1)
    insertion=wt[120:180]
    refs={'WT':{'ref_seq_with_itd':wt},'ITD_1':{'ref_seq_with_itd':wt[:180]+insertion+wt[180:]}}
    reads=pd.DataFrame([dict(read_id='exact60',read_seq=refs['ITD_1']['ref_seq_with_itd'],strand='+'),dict(read_id='short48',read_seq=wt[:180]+insertion[12:]+wt[180:],strand='+')])
    cons=pd.DataFrame([dict(peak_alias='ITD_1',median_ins_pos_ref=180,consensus_len=60)])
    out=align_reads_multi_ref_parallel(reads,refs,cons,logging.getLogger('probe'),threads=1)
    print(out[['read_id','ref_alias','pct_identity','adjusted_score','support_call','valid_support','filter_reason']].to_string(index=False))

if __name__=='__main__':
    main()
