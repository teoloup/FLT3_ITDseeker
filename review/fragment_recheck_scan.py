"""Audit marked fragment peaks independently of GMM/ASV candidate discovery.

Run in the project container with the repository mounted at /workspace.
"""
import os, sys, json, re, logging, hashlib, platform
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
import pysam
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0,'/workspace')
from itdseeker.bam_extractor import FLT3ReadExtractor
from itdseeker.GMM_peaks import fit_gmm_itds

ROOT=Path('/workspace'); OUT=ROOT/'FRAGMENT_DATA/recheck_20261019'
WT=re.search('DEFAULT_REF_WT = """(.*?)"""',(ROOT/'Nano_ITDseeker.py').read_text()).group(1)
START=28033881

def independent_aligner():
    a=Align.PairwiseAligner(); a.mode='global'; a.match_score=2; a.mismatch_score=-3
    a.open_gap_score=-8; a.extend_gap_score=-.5; a.end_gap_score=0
    return a

def parse_insertions(a,seq):
    c=a.align(WT,seq)[0].coordinates; events=[]
    for i in range(c.shape[1]-1):
        p,end=int(c[0,i]),int(c[0,i+1]); x,y=int(c[1,i]),int(c[1,i+1])
        if p==end and 3<=y-x<=300 and 20<p<len(WT)-20:
            events.append((p,y-x,seq[x:y],x,y))
    return events

def bam_evidence(bam_path):
    cigar=[]; meta={}; stats=Counter()
    with pysam.AlignmentFile(str(bam_path),'rb') as bam:
        chrom='chr13' if 'chr13' in bam.references else '13'
        rg=bam.header.to_dict().get('RG',[])
        for r in bam.fetch(chrom,28033300,28034800):
            stats['alignments']+=1
            if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
            stats['primary']+=1
            if r.mapping_quality<20: continue
            stats['primary_mapq20']+=1
            meta[r.query_name]=dict(bam_strand='-' if r.is_reverse else '+',mapq=r.mapping_quality,mean_q=float(np.mean(r.query_qualities)) if r.query_qualities is not None else 0)
            q=0; pos=r.reference_start; sequence=r.query_sequence or ''
            for op,n in r.cigartuples or []:
                if op==1:
                    if 3<=n<=300 and START+19<pos<START+len(WT)-20:
                        cigar.append(dict(read_id=r.query_name,genomic_anchor=pos,insertion_length=n,insertion_sequence=sequence[q:q+n],strand='-' if r.is_reverse else '+',mapq=r.mapping_quality))
                    q+=n
                elif op in (0,7,8): q+=n; pos+=n
                elif op in (2,3): pos+=n
                elif op==4: q+=n
    return cigar,meta,dict(stats),rg

def scan(entry):
    sample=entry['biological_sample']; p=OUT/sample; p.mkdir(exist_ok=True)
    if (p/'scan_complete.json').exists(): return json.loads((p/'scan_complete.json').read_text())
    logging.basicConfig(level=logging.WARNING)
    bam=ROOT/entry['bam']
    cigar,meta,bamstats,rg=bam_evidence(bam)
    pd.DataFrame(cigar).to_csv(p/'bam_cigar_insertions.tsv.gz',sep='\t',index=False)
    extractor=FLT3ReadExtractor(str(bam),'hg38',20,2,320,300,336,str(p/'fastq'))
    existing=p/'fastq/trimmed.fastq'
    if existing.exists() and existing.stat().st_size:
        reads={r.id:dict(seq=str(r.seq),qual=''.join(chr(q+33) for q in r.letter_annotations['phred_quality']),strand='-' if r.description.endswith(' rc') else '+') for r in SeqIO.parse(str(existing),'fastq')}
    else:
        reads,fq=extractor.process_reads()
    total=len(reads); a=independent_aligner(); cache={}; events=[]; lengths=Counter(); gmm_info={}
    # Default GMM replay only; all-read insertion discovery below bypasses it.
    default_reads={rid:r for rid,r in reads.items() if len(r['seq'])>=330}
    try:
        fit=fit_gmm_itds(default_reads,min_gmm_fraction=.01,max_itds_detected=6,min_ggmm_peak_distance=10,max_peak_sd=5,assign_width_factor=2,assign_mode='manual',prob_threshold=.85)
        fit.comps.to_csv(p/'default_gmm_components.tsv',sep='\t',index=False)
        gmm_info=dict(zip(fit.reads_df.read_id,fit.reads_df.gmm_peak_alias))
        fit.reads_df[['read_id','read_len','gmm_peak_alias','strand']].to_csv(p/'default_gmm_assignments.tsv.gz',sep='\t',index=False)
    except Exception as exc: (p/'gmm_error.txt').write_text(str(exc))
    for rid,r in reads.items():
        seq=r['seq']; lengths[len(seq)]+=1
        if seq not in cache: cache[seq]=parse_insertions(a,seq)
        quals=np.frombuffer(r['qual'].encode(),dtype=np.uint8).astype(float)-33
        for pos,n,ins,x,y in cache[seq]:
            allele=WT[:pos]+ins+WT[pos:]
            events.append(dict(read_id=rid,insertion_pos=pos,genomic_anchor=START+pos-1,insertion_length=n,insertion_sequence=ins,allele=allele,strand=r['strand'],read_len=len(seq),mean_q=float(np.mean(quals)),junction_q=float(np.mean(quals[max(0,x-10):min(len(quals),y+10)])),gmm_peak=gmm_info.get(rid,'excluded'),mapq=meta.get(rid,{}).get('mapq',0)))
    df=pd.DataFrame(events)
    df.to_csv(p/'all_read_insertions.tsv.gz',sep='\t',index=False)
    grouped=[]
    for allele,group in df.groupby('allele',sort=False):
        row=group.iloc[0]; n=int(row.insertion_length); pos=int(row.insertion_pos); ins=row.insertion_sequence
        adjacent=[]
        if pos>=n: adjacent.append(sum(x==y for x,y in zip(ins,WT[pos-n:pos]))/n)
        if pos+n<=len(WT): adjacent.append(sum(x==y for x,y in zip(ins,WT[pos:pos+n]))/n)
        cg=[r for r in cigar if r['insertion_length']==n and abs(r['genomic_anchor']-int(row.genomic_anchor))<=n+2]
        cigar_ids={r['read_id'] for r in cg if r['insertion_sequence']==ins}
        grouped.append(dict(sample=sample,allele_id=hashlib.sha256(allele.encode()).hexdigest()[:12],allele=allele,insertion_pos=pos,genomic_anchor=int(row.genomic_anchor),insertion_length=n,insertion_sequence=ins,n_reads=group.read_id.nunique(),plus=int(sum(group.strand=='+')),minus=int(sum(group.strand=='-')),mean_q=float(group.mean_q.median()),junction_q=float(group.junction_q.median()),q20_reads=int(sum(group.mean_q>=20)),gmm_peaks=json.dumps(dict(Counter(group.gmm_peak))),adjacent_dup_identity=max(adjacent or [0]),raw_cigar_exact_shared_ids=len(set(group.read_id)&cigar_ids),raw_cigar_length_nearby=len(cg)))
    gd=pd.DataFrame(grouped).sort_values('n_reads',ascending=False)
    gd.to_csv(p/'allele_groups.tsv',sep='\t',index=False)
    summary=dict(sample=sample,file_sample=entry['file_sample'],bam=entry['bam'],bam_stat=dict(size=bam.stat().st_size,mtime_ns=bam.stat().st_mtime_ns),bam_stats=bamstats,bam_read_groups=rg,trimmed_min320=total,trimmed_min330=len(default_reads),read_length_counts=dict(sorted(lengths.items())),default_gmm_counts=dict(Counter(gmm_info.values())),n_insertion_reads=int(df.read_id.nunique()),n_exact_alleles=len(gd))
    (p/'scan_complete.json').write_text(json.dumps(summary,indent=2))
    return summary

def main():
    manifest=json.loads((OUT/'manifest.json').read_text())
    with ProcessPoolExecutor(max_workers=3) as pool:
        jobs={pool.submit(scan,e):e['biological_sample'] for e in manifest['samples']}
        for future in as_completed(jobs):
            sample=jobs[future]
            try:
                summary=future.result(); print('SCAN_DONE',sample,summary['trimmed_min320'],'reads',summary['n_exact_alleles'],'exact insertion alleles',flush=True)
            except Exception as e:
                import traceback; traceback.print_exc(); print('SCAN_FAILED',sample,str(e),flush=True); raise

if __name__=='__main__': main()
