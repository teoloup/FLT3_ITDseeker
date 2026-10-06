import csv,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import simulate_itd_data as sim

def prepare():
    dest=ROOT/'review'/'synthetic'
    dest.mkdir(exist_ok=True)
    sim.SCENARIOS['no_wt']=[sim.Haplotype('ITD60',1.0,ins_pos=180,dup_len=60)]
    for scenario in ['no_wt','A']:
        paths=sim.simulate(scenario,2000,42,str(dest))
        errors=sim.self_check(scenario,paths,2000)
        if errors: raise RuntimeError(errors)

def run(method):
    rows=[]
    tests=sorted((ROOT/'bam_data'/'test_bam').glob('*.bam'))
    tests+=sorted((ROOT/'review'/'synthetic').glob('*.bam'))
    for variant,src in [('baseline',ROOT/'review'/'baseline'),('fixed',ROOT)]:
        for bam in tests:
            sample=bam.stem
            out=ROOT/'review'/'runs'/method/variant/sample
            out.mkdir(parents=True,exist_ok=True)
            cmd=[sys.executable,str(src/'Nano_ITDseeker.py'),'-b',str(bam),'-o',str(out),'-s',sample,'-t','4','--haplotype-method',method,'--min-allele-frequency','0.01']
            start=time.time()
            with (out/'console.log').open('w') as log:
                proc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=src)
            calls=[]
            vcf=out/f'{sample}_FLT3_ITD_calls.vcf'
            if vcf.exists():
                for line in vcf.read_text().splitlines():
                    if not line or line.startswith('#'): continue
                    fields=line.split('\t')
                    info=dict(x.split('=',1) for x in fields[7].split(';') if '=' in x)
                    calls.append(dict(position=fields[1],itd_len=info.get('ITD_LEN'),af=info.get('AF'),dp=info.get('DP')))
            row=dict(method=method,variant=variant,sample=sample,exit_code=proc.returncode,seconds=round(time.time()-start,1),calls=calls)
            rows.append(row)
            (ROOT/'review'/'runs'/f'{method}_summary.json').write_text(json.dumps(rows,indent=2))
            print(json.dumps(row),flush=True)
    if any(r['exit_code'] for r in rows): raise SystemExit(1)

if __name__=='__main__':
    if sys.argv[1]=='prepare': prepare()
    else: run(sys.argv[1])
