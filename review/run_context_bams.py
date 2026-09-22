import json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
inputs=sorted((ROOT/'bam_data/test_bam').glob('*.bam'))
# Keep the original validation set fixed; do not automatically include new inputs.
names={'10808_hg38_RG','11531_all_hg38_RG','13697_2runs_hg38_RG','14219_2runs_hg38_RG','14417_2runs_hg38_RG'}
inputs=[p for p in inputs if p.stem in names]
inputs += [ROOT/'review/synthetic/sim_no_wt.bam',ROOT/'review/synthetic/sim_A.bam']
summary_path=ROOT/"review/runs/context_dada2_summary.json"
rows=json.loads(summary_path.read_text()) if summary_path.exists() else []
completed={r["sample"] for r in rows if r["exit_code"]==0}
for bam in inputs:
    if bam.stem in completed: continue
    sample=bam.stem
    out=ROOT/'review/runs/context_dada2'/sample
    out.mkdir(parents=True,exist_ok=True)
    cmd=[sys.executable,str(ROOT/'Nano_ITDseeker.py'),'-b',str(bam),'-o',str(out),'-s',sample,'-t','4','--haplotype-method','dada2','--min-allele-frequency','0.01']
    start=time.time()
    with (out/'console.log').open('w') as log:
        code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT).returncode
    calls=[]
    vcf=out/f'{sample}_FLT3_ITD_calls.vcf'
    if vcf.exists():
        for line in vcf.read_text().splitlines():
            if not line or line.startswith('#'): continue
            fields=line.split('\t'); info=dict(x.split('=',1) for x in fields[7].split(';') if '=' in x)
            calls.append(dict(position=int(fields[1]),itd_len=int(info['ITD_LEN']),af=float(info['AF']),dp=int(info['DP']),ref=fields[3],alt=fields[4]))
    row=dict(sample=sample,exit_code=code,seconds=round(time.time()-start,1),calls=calls)
    rows.append(row)
    (ROOT/'review/runs/context_dada2_summary.json').write_text(json.dumps(rows,indent=2))
    print(json.dumps(row),flush=True)
assert all(r['exit_code']==0 for r in rows)
