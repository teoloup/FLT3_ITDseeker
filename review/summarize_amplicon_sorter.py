"""Summarize amplicon_sorter experiments without treating partial runs as negatives."""
from pathlib import Path
from collections import Counter
import hashlib
import csv
import re
import importlib.metadata
import json
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Bio import SeqIO
from Bio.Seq import Seq
import edlib
from itdseeker.Helper_functions import build_default_aligner,find_insertions
from scripts.simulate_itd_data import DEFAULT_REF_WT as WT
BASE=ROOT/'review/runs/amplicon_sorter'
rows=[]
for p in sorted(BASE.rglob('summary.json')):
    row=json.loads(p.read_text())
    folder=p.parent
    remaining=list(folder.rglob('*.group'))
    complete=row.get('completed',row['exit_code']==0 and not remaining and bool(list(folder.rglob('*_consensussequences.fasta'))))
    result={k:row.get(k) for k in ['sample','profile','peak','native','maxreads','length_diff_consensus','replicate','exit_code','seconds','input_reads','assigned_reads','unassigned_reads','duplicate_assignments','unknown_ids']}
    selection=re.search(r'(\d+) out of (\d+) sequences',(folder/'console.log').read_text())
    result['selected_reads']=int(selection.group(1)) if selection else None
    result.update(run=str(folder.relative_to(BASE)),completed=complete,clusters=[])
    for c in row['clusters']:
        seq=c['consensus']
        rc=str(Seq(seq).reverse_complement())
        if edlib.align(rc,WT)['editDistance']<edlib.align(seq,WT)['editDistance']:seq=rc
        insertions=[dict(position=p,length=n) for p,q,n in find_insertions(build_default_aligner().align(WT,seq)[0]) if n>=12] if seq else []
        result['clusters'].append(dict(n_reads=c['n_reads'],length=c['length'],insertions=insertions,
            truth_counts=c['truth_counts'],exact_truth=[k for k,v in c['truth_edit_distances'].items() if v==0],
            ambiguous_bases=sum(base not in 'ACGT' for base in seq)))
    rows.append(result)
(BASE/'comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
tool=ROOT/'review/tools/amplicon_sorter/amplicon_sorter.py'
manifest=dict(upstream='https://github.com/avierstr/amplicon_sorter',
    commit=subprocess.check_output(['git','-C',str(tool.parent),'rev-parse','HEAD'],text=True).strip(),
    script_sha256=hashlib.sha256(tool.read_bytes()).hexdigest(),python=sys.version,
    packages={name:importlib.metadata.version(name) for name in ['edlib','biopython','numpy','pysam']})
(BASE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for r in rows:
    print(r['run'], 'complete='+str(r['completed']), str(r['seconds'])+'s',
          'clusters='+str([(c['n_reads'],c['insertions'],c['truth_counts'],c['exact_truth']) for c in r['clusters']]))

# Independently check reconstructed VCF alleles, without rotation-only matching.
truth_alleles={}
for r in csv.DictReader((ROOT/'review/synthetic/truth_haplotypes_A.tsv').open(),delimiter='\t'):
    if int(r['itd_len']):
        p=int(r['ins_pos_local']);truth_alleles[r['haplotype']]=WT[:p]+r['itd_seq']+WT[p:]
vcf=BASE/'pipeline_replay/sim_A_FLT3_ITD_calls.vcf'
if vcf.exists():
    calls=[]
    for line in vcf.read_text().splitlines():
        if line.startswith('#'):continue
        f=line.split('\t');info=dict(x.split('=',1) for x in f[7].split(';') if '=' in x)
        p=int(f[1])-28033881+1;allele=WT[:p]+f[4][1:]+WT[p:]
        exact=[name for name,seq in truth_alleles.items() if allele==seq]
        assert len(exact)==1, 'Replay call is not an exact truth allele'
        calls.append(dict(haplotype=exact[0],alias=f[2],length=int(info['ITD_LEN']),af=float(info['AF']),dp=int(info['DP'])))
    assert {c['haplotype'] for c in calls}=={'ITD_A','ITD_B','ITD_C'}
    truth_reads={r['read_id']:r['haplotype'] for r in csv.DictReader((ROOT/'review/synthetic/truth_reads_A.tsv').open(),delimiter='\t')}
    support={}
    for r in csv.DictReader((BASE/'pipeline_replay/flt3_data/sim_A_validation_read_support.tsv').open(),delimiter='\t'):
        if r['valid_support']=='True':
            support.setdefault(r['ref_alias'],Counter())[truth_reads[r['read_id']]]+=1
    (BASE/'pipeline_replay_summary.json').write_text(json.dumps(dict(calls=calls,validated_read_truth=support),indent=2)+'\n')
    print('Verified exact replay alleles:',calls)
