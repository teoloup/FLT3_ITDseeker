"""Summarize before/after test BAM calls and check synthetic allele recovery."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from simulate_itd_data import DEFAULT_REF_WT as wt
new=json.loads((ROOT/'review/runs/context_dada2_summary.json').read_text())
old={r['sample']:r for r in json.loads((ROOT/'review/runs/dada2_summary.json').read_text()) if r['variant']=='fixed'}
assert len(new)==7 and all(r['exit_code']==0 for r in new)
lines=['| Sample | Previous calls: bp / AF % | Context calls: bp / AF % |','|---|---|---|']
for row in new:
    def fmt(calls):
        return ', '.join(f"{int(c['itd_len'])} / {float(c['af'])*100:.3f}" for c in sorted(calls,key=lambda c:int(c['itd_len']))) or 'No calls'
    lines.append(f"| {row['sample']} | {fmt(old[row['sample']]['calls'])} | {fmt(row['calls'])} |")
    if row['sample']=='sim_no_wt':
        assert len(row['calls'])==1
        call=row['calls'][0]; pos=call['position']-28033881+1
        assert wt[:pos]+call['alt'][1:]+wt[pos:]==wt[:180]+wt[120:180]+wt[180:]
        lines += []
print('\n'.join(lines))
(ROOT/'review/runs/context_comparison.md').write_text('\n'.join(lines)+'\n')
