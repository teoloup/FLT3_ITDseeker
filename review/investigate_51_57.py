"""Read-only FASTQ investigation; outputs are isolated from pipeline results."""
import sys, os, json, logging
from pathlib import Path
from collections import Counter
import pandas as pd
from Bio import Align

sys.path.insert(0, '/workspace')
from itdseeker.GMM_peaks import fit_gmm_itds
from itdseeker.haplotype_split import split_peaks

logging.basicConfig(level=logging.INFO)
out = Path('/workspace/review/fastq_51_57_analysis')
out.mkdir(exist_ok=True)
def readfq(path):
    reads = {}
    with open(path) as f:
        while True:
            h = f.readline().strip()
            if not h: break
            s = f.readline().strip(); f.readline(); q = f.readline().strip()
            reads[h.split()[0][1:]] = dict(seq=s, qual=q, strand='-' if h.endswith(' rc') else '+')
    return reads

reads = readfq('/inputs/trimmed.fastq')
seqcounts = Counter(r['seq'] for r in reads.values())
wt = seqcounts.most_common(1)[0][0]
aligner = Align.PairwiseAligner()
aligner.mode='global'; aligner.match_score=2; aligner.mismatch_score=-3
aligner.open_gap_score=-8; aligner.extend_gap_score=-0.5
def insertions(s):
    aln=aligner.align(wt,s)[0]
    c=aln.coordinates
    ins=[]
    for i in range(c.shape[1]-1):
        a,b=int(c[0,i]),int(c[0,i+1]); x,y=int(c[1,i]),int(c[1,i+1])
        if a==b and y-x>=12 and 20<a<len(wt)-20:
            ins.append((a,y-x,s[x:y]))
    return ins

events=Counter(); records=[]
for rid,r in reads.items():
    if len(r['seq']) < 365: continue
    ins=insertions(r['seq'])
    for pos,n,s in ins: events[(pos,n,s)]+=1
    records.append(dict(read_id=rid,read_len=len(r['seq']),strand=r['strand'],insertions=json.dumps(ins)))
pd.DataFrame(records).to_csv(out/'long_read_insertions.tsv',sep='\t',index=False)
print('TOP_INSERTIONS',json.dumps(events.most_common(15)),flush=True)
print('INSERTION_LENGTH_COUNTS',Counter({n:sum(count for (p,l,s),count in events.items() if l==n) for n in set(x[1] for x in events)}),flush=True)
for length in (387,393):
    top=Counter({s:n for s,n in seqcounts.items() if len(s)==length}).most_common(3)
    print('TOP_LENGTH',length,[(n,s,insertions(s)) for s,n in top],flush=True)

fit=fit_gmm_itds(reads,min_gmm_fraction=.01,max_itds_detected=6,min_ggmm_peak_distance=10,max_peak_sd=5,assign_width_factor=2,assign_mode='manual',prob_threshold=.85)
fit.comps.to_csv(out/'gmm_components.tsv',sep='\t',index=False)
fit.reads_df.to_csv(out/'gmm_read_assignments.tsv',sep='\t',index=False)
print('GMM',fit.comps.to_string(index=False),flush=True)
res=split_peaks('dada2',comps=fit.comps,reads_df=fit.reads_df,peak_subsets=fit.peak_subsets,wt_amplicon_length=336,threads=4,work_dir=str(out/'haplotypes'),min_child_fraction=.15,min_child_reads=20)
res.comps.to_csv(out/'dada2_components_default.tsv',sep='\t',index=False)
print('DADA2_DEFAULT',res.comps.to_string(index=False),flush=True)
pd.DataFrame([dict(insertion_position=p,insertion_length=l,insertion_sequence=s,read_count=n) for (p,l,s),n in events.most_common()]).to_csv(out/'insertion_counts.tsv',sep='\t',index=False)
