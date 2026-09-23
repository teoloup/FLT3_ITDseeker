"""Verify image smoke results against previously validated calls and audits."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'review/runs/docker_20260923'
expected = {r['sample']: r['calls'] for r in json.loads(
    (ROOT / 'review/runs/final_dada2_summary.json').read_text(encoding='utf-8'))}
expected['14417_2runs_hg38_RG'] = []
results = []
for sample in ['sim_no_wt', '10808_hg38_RG', '14417_2runs_hg38_RG']:
    folder = BASE / sample
    calls = []
    for line in (folder / f'{sample}_FLT3_ITD_calls.vcf').read_text(encoding='utf-8').splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
        calls.append(dict(position=int(fields[1]), itd_len=int(info['ITD_LEN']),
                          af=float(info['AF']), dp=int(info['DP']), ref=fields[3], alt=fields[4]))
    assert calls == expected[sample], (sample, calls, expected[sample])
    commands = [json.loads(line) for line in
                (folder / f'{sample}_commands.jsonl').read_text(encoding='utf-8').splitlines()]
    assert [Path(c['argv'][0]).name for c in commands[:3]] == ['samtools', 'samtools', 'cutadapt']
    if calls:
        assert any(Path(c['argv'][0]).name == 'Rscript' for c in commands)
    assert list(folder.glob('*.html')), sample
    results.append(dict(sample=sample, calls=calls, audited_commands=len(commands)))
(BASE / 'checked_results.json').write_text(json.dumps(results, indent=2) + '\n', encoding='utf-8')
print(json.dumps(results, indent=2))
