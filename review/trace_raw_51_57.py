from collections import Counter
from pathlib import Path
import json
from Bio import Align
from Bio.Seq import Seq
import pandas as pd

p=Path('/workspace/review/fastq_51_57_analysis')
f=open('/inputs/trimmed.fastq'); trimmed={}; seqs=Counter()
while True:
 h=f.readline().strip()
 if not h: break
 s=f.readline().strip(); f.readline(); f.readline()
 trimmed[h.split()[0][1:]]=s; seqs[s]+=1
wt=seqs.most_common(1)[0][0]
a=Align.PairwiseAligner(); a.mode='global'; a.match_score=2; a.mismatch_score=-3
a.open_gap_score=-8; a.extend_gap_score=-.5; a.end_gap_score=0
f=open('/inputs/region_reads.fastq'); records=[]; total=0
while True:
 h=f.readline().strip()
 if not h: break
 s=f.readline().strip(); f.readline(); f.readline(); total+=1
 if len(s)<375: continue
 rid=h.split()[0][1:]; rc=str(Seq(s).reverse_complement())
 if a.score(wt,rc)>a.score(wt,s): s=rc
 c=a.align(wt,s)[0].coordinates
 for i in range(c.shape[1]-1):
  l,r=int(c[0,i]),int(c[0,i+1]); x,y=int(c[1,i]),int(c[1,i+1])
  if l==r and y-x in (51,57) and 20<l<len(wt)-20:
   records.append(dict(read_id=rid,insertion_length=y-x,insertion_sequence=s[x:y],retained_after_trimming=rid in trimmed))
df=pd.DataFrame(records); df.to_csv(p/'raw_read_support.tsv',sep='\t',index=False)
summary={'raw_total_reads':total,'raw_exact_length_support':df.groupby(['insertion_length','retained_after_trimming']).size().to_string()}
print(json.dumps(summary,indent=2))
(p/'raw_summary.json').write_text(json.dumps(summary,indent=2))
