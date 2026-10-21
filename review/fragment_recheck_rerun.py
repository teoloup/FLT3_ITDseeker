"""Run review-only parameter sensitivity tests; preserve original reports."""
import sys,subprocess,json,time
from pathlib import Path
from fragment_recheck_scan import ROOT,OUT
manifest=json.loads((OUT/'manifest.json').read_text())
baseline='--baseline' in sys.argv
samples=[x for x in sys.argv[1:] if not x.startswith('--')] or ['11528','13462']
for entry in manifest['samples']:
    sample=entry['biological_sample']
    if sample not in samples: continue
    p=OUT/sample/('baseline_rerun' if baseline else 'sensitive_rerun'); p.mkdir(exist_ok=True)
    if (p/'run_exit.json').exists(): continue
    cmd=[sys.executable,str(ROOT/'Nano_ITDseeker.py'),'-b',str(ROOT/entry['bam']),'-o',str(p/'results'),'-s',sample+'_review','-g','hg38','-t','3','--min-allele-frequency','0.001','--min-gmm-fraction','0.001','--max-peak-sd','10','--min-subpeak-fraction','0.01','--min-haplotype-reads','10','--max-itds-to-detect','10','--min-read-length','320','--min-itd-size','3','--wt-peak-tolerance','1','--html-report']
    if baseline:
        cmd=[sys.executable,str(ROOT/'Nano_ITDseeker.py'),'-b',str(ROOT/entry['bam']),'-o',str(p/'results'),'-s',sample+'_baseline','-g','hg38','-t','3','--min-allele-frequency','0.005','--html-report']
    (p/'command.json').write_text(json.dumps(cmd,indent=2))
    start=time.time()
    with open(p/'run.log','w') as log:
        res=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    (p/'run_exit.json').write_text(json.dumps(dict(returncode=res.returncode,elapsed_seconds=time.time()-start),indent=2))
    print('SENSITIVE_DONE',sample,'exit',res.returncode,'seconds',round(time.time()-start,1),flush=True)
