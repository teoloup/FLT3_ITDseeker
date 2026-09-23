"""Prepare matched-read consensus comparisons; run only with an explicit Medaka model.

Execute with the Linux validation Python, whose dependencies include Bio/pysam.
Medaka itself stays isolated in /tmp/flt3-medaka-review/venv.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Bio import SeqIO
import edlib
from command_audit import configure_command_log,record_command
from Helper_functions import build_default_aligner,find_insertions
from simulate_itd_data import DEFAULT_REF_WT as WT
BASE=ROOT/'review/runs/medaka'
AS=ROOT/'review/runs/amplicon_sorter'


def msa_draft(sample,length):
    path=ROOT/'review/runs/context_dada2'/sample/'flt3_data'/f'{sample}_itd_consensus_seq.tsv'
    rows=[r for r in csv.DictReader(path.open(),delimiter='\t') if int(r['consensus_len'])==length]
    assert len(rows)==1
    r=rows[0];p=int(r['consensus_ins_pos_ref'])
    return WT[:p]+r['consensus_seq']+WT[p:]


def prepare():
    folder=BASE/'inputs';folder.mkdir(parents=True,exist_ok=True)
    rows=[]
    def add(name,reads,draft,synthetic=False):
        path=folder/(name+'.fasta');path.write_text('>draft\n'+draft+'\n')
        ids=[r.id for r in SeqIO.parse(reads,'fastq')]
        assert ids and len(ids)==len(set(ids))
        rows.append(dict(name=name,reads=str(reads),draft=str(path),n_reads=len(ids),synthetic=synthetic,
                         reads_sha256=hashlib.sha256(reads.read_bytes()).hexdigest(),
                         draft_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    real=json.loads((AS/'13697_2runs_hg38_RG/sc99__ITD_2__native/summary.json').read_text())
    assert real['completed'] and len(real['clusters'])==1 and real['input_reads']==589
    reads=AS/'inputs/13697_2runs_hg38_RG/ITD_2.fastq'
    add('real72_sorter',reads,real['clusters'][0]['consensus'])
    add('real72_msa',reads,msa_draft('13697_2runs_hg38_RG',72))
    original=json.loads((AS/'13697_2runs_hg38_RG/sc99__native/summary.json').read_text())
    short=next(c['consensus'] for c in original['clusters'] if c['length']==405)
    add('real72_short_draft',reads,short)
    pure=json.loads((AS/'sim_no_wt/default__native/summary.json').read_text())
    assert pure['completed'] and len(pure['clusters'])==1
    reads=AS/'inputs/sim_no_wt/trimmed.fastq'
    add('sim60_sorter',reads,pure['clusters'][0]['consensus'],True)
    add('sim60_msa',reads,msa_draft('sim_no_wt',60),True)
    (BASE/'prepared_cases.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--model',help='Matching Medaka model or basecaller-model:consensus; never guessed')
    p.add_argument('--case',action='append',help='Restrict to named cases; repeat this option if needed')
    p.add_argument('--timeout',type=int,default=300)
    p.add_argument('--medaka-bin',type=Path,default=Path('/tmp/flt3-medaka-review/venv/bin'))
    p.add_argument('--tools-bin',type=Path,default=Path('/tmp/flt3-medaka-review/tools/bin'))
    args=p.parse_args()
    if not args.prepare_only and not args.model:p.error('Supply a confirmed --model, or use --prepare-only')
    if args.timeout<1:p.error('timeout must be positive')
    cases=prepare()
    if args.case:
        unknown=set(args.case)-{c['name'] for c in cases}
        if unknown:p.error('Unknown cases: '+str(sorted(unknown)))
        cases=[c for c in cases if c['name'] in args.case]
    if args.prepare_only:
        print(json.dumps(cases,indent=2));return
    env=os.environ.copy()
    env['PATH']=str(args.medaka_bin)+':'+str(args.tools_bin)+':/usr/bin:/bin'
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',SHELLOPTS='xtrace')
    # Bash tracing records commands inside medaka_consensus and its child scripts.
    model_tag=hashlib.sha256(args.model.encode()).hexdigest()[:12]
    for case in cases:
        out=BASE/'runs'/model_tag/case['name']
        if out.exists():raise FileExistsError(f'Refusing to mix results with an existing run: {out}')
        out.mkdir(parents=True)
        configure_command_log(out/'commands.jsonl')
        with tempfile.TemporaryDirectory(prefix='flt3-medaka-benchmark-') as tmp:
            work=Path(tmp)/'out'
            cmd=['/bin/bash','-x',str(args.medaka_bin/'medaka_consensus'),'-i',case['reads'],
                 '-d',case['draft'],'-o',str(work),'-m',args.model,'-t','2','-b','8']
            record_command(cmd,stdin='inherited',stdout=str(out/'console.log'),stderr='merged-with-stdout')
            start=time.monotonic()
            with (out/'console.log').open('w') as log:
                proc=subprocess.Popen(cmd,env=env,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                try:code=proc.wait(timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGTERM)
                    try:proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                    code='timeout'
            if work.exists():shutil.copytree(work,out/'output')
        result=dict(case=case,model=args.model,exit_code=code,seconds=round(time.monotonic()-start,2))
        polished=out/'output/consensus.fasta'
        records=list(SeqIO.parse(polished,'fasta')) if polished.exists() else []
        result['completed']=code==0 and len(records)==1 and bool(records[0].seq)
        if result['completed']:
            seq=str(records[0].seq);draft=str(next(SeqIO.parse(case['draft'],'fasta')).seq)
            a=build_default_aligner().align(WT,seq)[0]
            result.update(sequence=seq,length=len(seq),edit_distance_from_draft=edlib.align(seq,draft)['editDistance'],
                          insertions=[dict(position=p,length=n,sequence=seq[q:q+n]) for p,q,n in find_insertions(a) if n>=12])
            if case['synthetic']:
                truth=WT[:180]+WT[120:180]+WT[180:]
                result.update(exact_synthetic_allele=seq==truth,truth_edit_distance=edlib.align(seq,truth)['editDistance'])
        (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2),flush=True)
        if not result['completed']:raise RuntimeError(f'Medaka did not complete: see {out}')

if __name__=='__main__':main()
