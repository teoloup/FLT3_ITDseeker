"""Cross-check all-read insertion discovery with the production aligner."""
import sys,json,logging,hashlib
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
sys.path.insert(0,'/workspace')
from itdseeker.Helper_functions import build_default_aligner,find_insertions
from fragment_recheck_scan import OUT,WT,START

def scan(entry):
    sample=entry['biological_sample']; p=OUT/sample
    if (p/'pipeline_allele_groups.tsv').exists(): return sample
    gmm=pd.read_csv(p/'default_gmm_assignments.tsv.gz',sep='\t')
    peaks=dict(zip(gmm.read_id,gmm.gmm_peak_alias)); a=build_default_aligner(); cache={}; events=[]
    for rec in SeqIO.parse(str(p/'fastq/trimmed.fastq'),'fastq'):
        seq=str(rec.seq); qual=rec.letter_annotations['phred_quality']; strand='-' if rec.description.endswith(' rc') else '+'
        if seq not in cache:
            aln=a.align(WT,seq)[0]
            cache[seq]=[(int(pos),int(x),int(n)) for pos,x,n in find_insertions(aln) if 3<=n<=300 and 20<pos<len(WT)-20]
        for pos,x,n in cache[seq]:
            ins=seq[x:x+n]; allele=WT[:pos]+ins+WT[pos:]
            events.append(dict(read_id=rec.id,insertion_pos=pos,genomic_anchor=START+pos-1,insertion_length=n,insertion_sequence=ins,allele=allele,strand=strand,read_len=len(seq),mean_q=float(np.mean(qual)),junction_q=float(np.mean(qual[max(0,x-10):min(len(qual),x+n+10)])),gmm_peak=peaks.get(rec.id,'excluded')))
    df=pd.DataFrame(events); df.to_csv(p/'pipeline_all_read_insertions.tsv.gz',sep='\t',index=False)
    grouped=[]
    cigar=pd.read_csv(p/'bam_cigar_insertions.tsv.gz',sep='\t'); cigar_by=defaultdict(set)
    for row in cigar.itertuples(): cigar_by[(row.insertion_length,row.insertion_sequence)].add(row.read_id)
    for allele,group in df.groupby('allele',sort=False):
        row=group.iloc[0]; pos=int(row.insertion_pos); n=int(row.insertion_length); ins=row.insertion_sequence
        adjacent=[]
        if pos>=n: adjacent.append(sum(x==y for x,y in zip(ins,WT[pos-n:pos]))/n)
        if pos+n<=len(WT): adjacent.append(sum(x==y for x,y in zip(ins,WT[pos:pos+n]))/n)
        grouped.append(dict(sample=sample,allele_id=hashlib.sha256(allele.encode()).hexdigest()[:12],allele=allele,insertion_pos=pos,genomic_anchor=int(row.genomic_anchor),insertion_length=n,insertion_sequence=ins,n_reads=group.read_id.nunique(),plus=int(sum(group.strand=='+')),minus=int(sum(group.strand=='-')),mean_q=float(group.mean_q.median()),junction_q=float(group.junction_q.median()),q20_reads=int(sum(group.mean_q>=20)),gmm_peaks=json.dumps(dict(Counter(group.gmm_peak))),adjacent_dup_identity=max(adjacent or [0]),raw_cigar_exact_shared_ids=len(set(group.read_id)&cigar_by[(n,ins)])))
    pd.DataFrame(grouped).sort_values('n_reads',ascending=False).to_csv(p/'pipeline_allele_groups.tsv',sep='\t',index=False)
    return sample

from collections import defaultdict
if __name__=='__main__':
    logging.basicConfig(level=logging.WARNING)
    manifest=json.loads((OUT/'manifest.json').read_text())
    entries=[e for e in manifest['samples'] if (OUT/e['biological_sample']/'scan_complete.json').exists() and (len(sys.argv)==1 or e['biological_sample'] in sys.argv[1:])]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for f in as_completed([pool.submit(scan,e) for e in entries]): print('PIPELINE_SCAN_DONE',f.result(),flush=True)
