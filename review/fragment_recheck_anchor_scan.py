"""Measure net insertion size between two exact outer anchors, before alignment."""
import sys,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from fragment_recheck_scan import OUT,WT,START,independent_aligner,parse_insertions

manifest=json.loads((OUT/'manifest.json').read_text())
left=WT[25:40]; right=WT[300:315]; refwindow=WT[25:315]
for entry in manifest['samples']:
    sample=entry['biological_sample']; p=OUT/sample
    if (p/'anchored_all_read_sizes.tsv.gz').exists(): continue
    counts=defaultdict(list); n=0; all_sizes=[]
    for rec in SeqIO.parse(str(p/'fastq/trimmed.fastq'),'fastq'):
        n+=1; seq=str(rec.seq); x=seq.find(left); y=seq.rfind(right)
        if x<0 or y<x+15: continue
        window=seq[x:y+15]; net=len(window)-len(refwindow)
        if not -20<=net<=300: continue
        item=dict(read_id=rec.id,strand='-' if rec.description.endswith(' rc') else '+',mean_q=sum(rec.letter_annotations['phred_quality'])/len(seq))
        counts[window].append(item)
        all_sizes.append(dict(net_insertion_length=net,window_hash=hashlib.sha256(window.encode()).hexdigest()[:12],**item))
    rows=[]; idrows=[]; a=independent_aligner()
    for window,reads in counts.items():
        net=len(window)-len(refwindow)
        if net<3 or len(reads)<3: continue
        allele=WT[:25]+window+WT[315:]; aid=hashlib.sha256(allele.encode()).hexdigest()[:12]
        events=parse_insertions(a,allele)
        main=max(events,key=lambda e:e[1]) if events else (0,0,'',0,0)
        pos,inslen,ins,x,y=main
        rows.append(dict(sample=sample,allele_id=aid,allele=allele,net_insertion_length=net,insertion_length=inslen,insertion_pos=pos,insertion_sequence=ins,genomic_anchor=START+pos-1,n_reads=len(reads),plus=sum(r['strand']=='+' for r in reads),minus=sum(r['strand']=='-' for r in reads),mean_q=sorted(r['mean_q'] for r in reads)[len(reads)//2],alignment_insertion_events=json.dumps([e[:3] for e in events])))
        idrows.extend(dict(allele_id=aid,**r) for r in reads)
    pd.DataFrame(rows).sort_values('n_reads',ascending=False).to_csv(p/'anchored_haplotypes.tsv',sep='\t',index=False)
    pd.DataFrame(idrows).to_csv(p/'anchored_read_ids.tsv.gz',sep='\t',index=False)
    pd.DataFrame(all_sizes).to_csv(p/'anchored_all_read_sizes.tsv.gz',sep='\t',index=False)
    (p/'anchor_stats.json').write_text(json.dumps(dict(total_reads=n,exact_both_anchor_reads=sum(len(v) for v in counts.values()),left_anchor=left,right_anchor=right,WT_window_length=len(refwindow)),indent=2))
    print('ANCHOR_DONE',sample,'anchored',sum(len(v) for v in counts.values()),'of',n,flush=True)
