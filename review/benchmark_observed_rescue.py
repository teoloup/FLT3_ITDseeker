"""Matched full pipeline comparisons; run inside the project container.

Each worker runs alone in a fresh Python process to measure elapsed time and
Linux ru_maxrss (largest process, not summed simultaneous worker RSS).
"""
import json
import resource
import subprocess
import sys
import time
import tarfile
import tempfile
from pathlib import Path

ROOT = Path('/workspace')
OUT = ROOT / 'review/observed_rescue_benchmark'

if len(sys.argv) > 1 and sys.argv[1] == '--one':
    sample, mode = sys.argv[2:4]
    manifest = json.loads((ROOT/'FRAGMENT_DATA/recheck_20261019/manifest.json').read_text())
    entry = next(r for r in manifest['samples'] if r['biological_sample'] == sample)
    target = OUT/sample/mode
    target.mkdir(parents=True, exist_ok=True)
    source = ROOT
    baseline_temp = None
    if mode == 'original':
        baseline_temp = tempfile.TemporaryDirectory(prefix='flt3-baseline-')
        source = Path(baseline_temp.name)
        # Trusted local git archive of the pre-change tracked pipeline source.
        with tarfile.open(ROOT/'review/observed_rescue_baseline.tar') as archive:
            archive.extractall(source, filter='data')
    command = [sys.executable, str(source/'Nano_ITDseeker.py'), '-b', str(ROOT/entry['bam']),
               '-o',str(target/'results'),'-s',sample,'-t','3',
               '--min-allele-frequency','.005','--html-report']
    if mode == 'disabled':
        command += ['--disable-observed-rescue']
    if mode == 'no_gmm':
        command += ['--max-peak-sd','0']
    start = time.perf_counter()
    with (target/'run.log').open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    summary = dict(sample=sample, mode=mode, command=command, returncode=result.returncode,
                   wall_seconds=time.perf_counter()-start,
                   largest_child_maxrss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    (target/'metrics.json').write_text(json.dumps(summary,indent=2))
    if baseline_temp:
        baseline_temp.cleanup()
    print(json.dumps(summary),flush=True)
    sys.exit(result.returncode)

OUT.mkdir(exist_ok=True)
for sample in sys.argv[1:] or ['11528','13462','11551']:
    for mode in ['original','disabled','enabled']:
        path = OUT/sample/mode/'metrics.json'
        if path.exists() and json.loads(path.read_text())['returncode'] == 0:
            continue
        result = subprocess.run([sys.executable,__file__,'--one',sample,mode])
        if result.returncode:
            raise SystemExit(result.returncode)
