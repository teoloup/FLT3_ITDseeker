import csv,json,sys
from pathlib import Path
from collections import Counter,defaultdict
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import edlib
from Helper_functions import build_default_aligner,find_insertions
from simulate_itd_data import DEFAULT_REF_WT as WT
p=ROOT/'review/runs/amplicon_sorter/13697_2runs_hg38_RG/sc99__native/summary.json'
r=json.loads(p.read_text());ins=[]
for c in r['clusters']:
    seq=c['consensus'];a=build_default_aligner().align(WT,seq)[0]
    events=[(p,n,seq[q:q+n]) for p,q,n in find_insertions(a) if n>=12]
    print(c['n_reads'],c['length'],events)
    if any(e[1]==72 for e in events):ins.append((seq,events))
print('72bp full allele distance',edlib.align(ins[0][0],ins[1][0])['editDistance'])
print('72bp insertion payloads equal',ins[0][1]==ins[1][1])
truth={r['read_id']:r['haplotype'] for r in csv.DictReader((ROOT/'review/synthetic/truth_reads_A.tsv').open(),delimiter='\t')}
counts=defaultdict(Counter)
for r in csv.DictReader((ROOT/'review/runs/amplicon_sorter/pipeline_replay/flt3_data/sim_A_validation_read_support.tsv').open(),delimiter='\t'):
    if r['valid_support']=='True': counts[r['ref_alias']][truth[r['read_id']]]+=1
print('Replay validated read truth:',dict(counts))
