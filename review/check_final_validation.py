"""Check final smoke results and persisted command audits (run after validation)."""
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.simulate_itd_data import DEFAULT_REF_WT as WT

rows = []
for sample in ['sim_A', 'sim_no_wt', '10808_hg38_RG']:
    folder = ROOT / 'review/runs/final_dada2' / sample
    calls = []
    for line in (folder / f'{sample}_FLT3_ITD_calls.vcf').read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
        assert not set(fields[4]) - set('ACGT'), fields[4]
        calls.append(dict(position=int(fields[1]), itd_len=int(info['ITD_LEN']),
                          af=float(info['AF']), dp=int(info['DP']), ref=fields[3], alt=fields[4]))
    audit = [json.loads(line) for line in (folder / f'{sample}_commands.jsonl').read_text().splitlines()]
    assert [Path(r['argv'][0]).name for r in audit[:3]] == ['samtools', 'samtools', 'cutadapt']
    assert audit[0]['argv'][1] == 'view' and audit[1]['argv'][1] == 'fastq'
    r_calls = [r for r in audit if Path(r['argv'][0]).name == 'Rscript']
    assert r_calls and all(r['argv'][-5:] == ['1e-40', '32', '-1.0', '20', '4'] for r in r_calls)
    assert all(r['shell'] is False and r['executable'] for r in audit)
    assert not (folder / 'temp').exists(), 'Audit must survive temporary-file cleanup'
    rows.append(dict(sample=sample, calls=calls, audited_commands=len(audit)))

pure = rows[1]['calls']
assert len(pure) == 1 and pure[0]['itd_len'] == 60 and pure[0]['dp'] == 1536
p = pure[0]['position'] - 28033881 + 1
assert WT[:p] + pure[0]['alt'][1:] + WT[p:] == WT[:180] + WT[120:180] + WT[180:]
assert [c['itd_len'] for c in rows[0]['calls']] == [30]
previous = json.loads((ROOT / 'review/runs/context_dada2_summary.json').read_text())
expected = next(r['calls'] for r in previous if r['sample'] == '10808_hg38_RG')
assert rows[2]['calls'] == expected
out = ROOT / 'review/runs/final_dada2_summary.json'
out.write_text(json.dumps(rows, indent=2) + '\n')
print(json.dumps(rows, indent=2))
