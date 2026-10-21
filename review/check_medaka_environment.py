import json,sys,importlib.util
from pathlib import Path
import pysam
root=Path('/mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker')
print('medaka',importlib.util.find_spec('medaka'),'torch',importlib.util.find_spec('torch'))
for p in sorted((root/'bam_data/test_bam').glob('*.bam')):
    with pysam.AlignmentFile(p,'rb') as b:
        h=b.header.to_dict()
        print(p.name,json.dumps({k:h.get(k) for k in ['RG','PG']})[:6000])
for p in sorted((root/'review/runs/amplicon_sorter/13697_2runs_hg38_RG').glob('*/summary.json')):
    r=json.loads(p.read_text());print(p.parent.name,r.get('completed'),r['seconds'],[(c['n_reads'],c['length']) for c in r['clusters']])
