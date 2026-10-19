"""Create an explicit workbook-row to physical-BAM mapping without editing input."""
from pathlib import Path
import json, hashlib
import openpyxl

root=Path(__file__).resolve().parents[1]
out=root/'FRAGMENT_DATA'/'recheck_20261019'
out.mkdir(exist_ok=True)
book=root/'FRAGMENT_DATA'/'COMPARISON.xlsx'
wb=openpyxl.load_workbook(book,data_only=True)
headers=list(next(wb.active.values)); rows=[]
for n,values in enumerate(list(wb.active.values)[1:],2):
    if values[0] is None: continue
    row=dict(zip(headers,values)); row['excel_row']=n
    row['biological_sample']=str(int(values[0]))
    row['file_sample']={'11094':'11202','11202':'11094'}.get(row['biological_sample'],row['biological_sample'])
    flagged=str(row.get('recheck','')).strip().upper()=='YES'
    below=values[3] is None and float(values[1])<12
    swap_control=row['biological_sample']=='11202' and values[3] is None
    row['review_reason']='flagged' if flagged else 'additional_below_minimum' if below else 'swap_control' if swap_control else ''
    rows.append(row)
samples={}
for row in rows:
    if not row['review_reason']: continue
    sample=row['biological_sample']; file_sample=row['file_sample']
    candidates=[p for p in (root/'bam_data'/'all_bam').glob(file_sample+'_*.bam') if not p.name.startswith((file_sample+'_1_',file_sample+'_5_'))]
    preferred=[p for p in candidates if any(x in p.stem for x in ('_2runs','_merged','_all_'))]
    selected=preferred or candidates
    if len(selected)!=1: raise ValueError((sample,selected))
    samples.setdefault(sample,dict(biological_sample=sample,file_sample=file_sample,bam=selected[0].relative_to(root).as_posix(),targets=[],known_calls=[]))['targets'].append(row)
for sample,entry in samples.items():
    entry['known_calls']=[row for row in rows if row['biological_sample']==sample and row['Nanopore size (bp)'] is not None]
manifest=dict(workbook=str(book),workbook_sha256=hashlib.sha256(book.read_bytes()).hexdigest(),swap_mapping={'11094':'11202','11202':'11094'},rows=rows,samples=list(samples.values()))
(out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Flagged',sum(r['review_reason']=='flagged' for r in rows),'total targets',sum(len(s['targets']) for s in samples.values()),'BAMs',len(samples))
for s in samples.values(): print(s['biological_sample'],'->',s['bam'],[(r['excel_row'],round(r['Fragment size (bp)'],2),r['review_reason']) for r in s['targets']])
