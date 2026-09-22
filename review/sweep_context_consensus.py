import json,sys,tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
import simulate_itd_data as sim
from Helper_functions import process_chunk
from Multiple_seq_aligment_toolkit import build_itd_consensus_sequences

def case(args):
    pos,scale=args
    ref=sim.DEFAULT_REF_WT
    truth=ref[:pos]+ref[pos-60:pos]+ref[pos:]
    sim.ERR_MISMATCH=.00569*scale
    sim.ERR_INSERTION=.00298*scale
    sim.ERR_DELETION=.00435*scale
    rng=np.random.default_rng(927+pos)
    reads=[(f'r{i}',sim.apply_errors(truth,rng)[0],'+') for i in range(400)]
    ins=pd.DataFrame(process_chunk(reads,54,66,ref,'ITD_1'))
    comps=pd.DataFrame([dict(peak_alias='ITD_1',putative_itd_size=60,sd_bp=2)])
    with tempfile.TemporaryDirectory() as out:
        cons=build_itd_consensus_sequences(ins,comps,ref_seq=ref,out_dir=out,sample_name='sweep',min_col_coverage=.7,max_unique=150,threads=2)
    row=dict(position=pos,error_scale=scale,input_reads=400,extracted=len(ins),called_length=None,exact_allele=False)
    if not cons.empty:
        c=cons.iloc[0]; p=int(c.consensus_ins_pos_ref)
        row.update(called_length=int(c.consensus_len),anchor=p,exact_allele=ref[:p]+c.consensus_seq+ref[p:]==truth)
    return row

if __name__=='__main__':
    cases=[(p,e) for p in [120,180,240] for e in [0,1,2]]
    rows=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(case,cases):
            rows.append(row); print(json.dumps(row),flush=True)
            (ROOT/'review/runs/context_sweep.json').write_text(json.dumps(rows,indent=2))
    if not all(r['exact_allele'] for r in rows): raise SystemExit(1)
