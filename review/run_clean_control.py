import json,sys,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import scripts.simulate_itd_data as sim
sim.ERR_MISMATCH=sim.ERR_INSERTION=sim.ERR_DELETION=0.
sim.SCENARIOS['clean_no_wt']=[sim.Haplotype('ITD60',1.,ins_pos=180,dup_len=60)]
paths=sim.simulate('clean_no_wt',400,42,str(ROOT/'review'/'synthetic'))
assert not sim.self_check('clean_no_wt',paths,400)
rows=[]
for variant,src in [('baseline',ROOT/'review'/'baseline'),('fixed',ROOT)]:
    out=ROOT/'review'/'runs'/'clean_control'/variant
    out.mkdir(parents=True,exist_ok=True)
    cmd=[sys.executable,str(src/'Nano_ITDseeker.py'),'-b',paths['bam'],'-o',str(out),'-s','clean_no_wt','-t','4','--haplotype-method','dada2','--force-number-of-peaks','1']
    with (out/'console.log').open('w') as log:
        code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=src).returncode
    vcf=out/'clean_no_wt_FLT3_ITD_calls.vcf'
    calls=[line for line in vcf.read_text().splitlines() if line and not line.startswith('#')] if vcf.exists() else []
    rows.append(dict(variant=variant,exit_code=code,calls=calls))
print(json.dumps(rows,indent=2))
(ROOT/'review'/'runs'/'clean_control_summary.json').write_text(json.dumps(rows,indent=2))
assert all(r['exit_code']==0 for r in rows)
assert not rows[0]['calls']
assert len(rows[1]['calls'])==1 and 'ITD_LEN=60;' in rows[1]['calls'][0]
