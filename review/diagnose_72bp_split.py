"""Locate the 72 bp split difference and save read-level evidence."""
from pathlib import Path
from collections import Counter
import importlib.util
import json
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Bio import SeqIO
from Bio.Seq import Seq
import edlib
import pysam
BASE=ROOT/'review/runs/amplicon_sorter'
folder=BASE/'13697_2runs_hg38_RG/sc99__native'
cs=sorted([c for c in json.loads((folder/'summary.json').read_text())['clusters'] if c['length']>400],key=lambda c:c['length'])
short,long=[c['consensus'] for c in cs]
spec=importlib.util.spec_from_file_location('sorter',ROOT/'review/tools/amplicon_sorter/amplicon_sorter.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
inputs=BASE/'inputs/13697_2runs_hg38_RG'
raw={r.id:str(r.seq) for r in SeqIO.parse(inputs/'region_reads.fastq','fastq')}
indexes={r.id:i for i,r in enumerate(SeqIO.parse(inputs/'trimmed.fastq','fastq'))}
with pysam.AlignmentFile(ROOT/'bam_data/test_bam/13697_2runs_hg38_RG.bam','rb') as bam:
    meta={r.query_name:dict(reverse=r.is_reverse,rg=r.get_tag('RG') if r.has_tag('RG') else None)
          for r in bam.fetch('chr13',28033300,28034800) if not r.is_secondary and not r.is_supplementary}
rows=[]
for c in cs:
    reads=list(SeqIO.parse(folder/c['file'],'fastq'))
    starts=Counter();sources=Counter();batches=Counter();strands=Counter()
    for r in reads:
        seq=str(r.seq);starts[seq[:12]]+=1;batches[indexes[r.id]//1000]+=1;strands[str(meta[r.id]['reverse'])]+=1
        if seq.startswith('TACCTTTCAGCA'):
            original=raw[r.id]
            if seq not in original:original=str(Seq(original).reverse_complement())
            pos=original.find(seq);assert pos>=0
            sources['at_raw_start' if pos==0 else 'internal_cutadapt_boundary']+=1
    rows.append(dict(file=c['file'],reads=len(reads),length=c['length'],start_counts=dict(starts),
                     short_prefix_origin=dict(sources),input_batches=dict(batches),bam_reverse_counts=dict(strands)))
result=dict(short_is_exact_suffix=long[3:]==short,removed_prefix=long[:3],
            nw_edit_distance=edlib.align(short,long)['editDistance'],
            tool_HW_similarity=module.distance(short,long,mode='HW'),
            tool_NW_similarity=module.distance(short,long,mode='NW'),clusters=rows)
assert result['short_is_exact_suffix'] and result['tool_HW_similarity']==1
(BASE/'72bp_split_diagnosis.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
