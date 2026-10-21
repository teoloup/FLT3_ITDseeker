"""Compact benchmark evidence without per-read sequences or large reports."""
import ast
import csv
import json
import re
from collections import Counter
from pathlib import Path

root=Path(__file__).resolve().parents[1]
folder=root/'review/observed_rescue_benchmark'
results=[]
for metrics in sorted(folder.glob('*/*/metrics.json')):
    row=json.loads(metrics.read_text())
    if row['returncode']:
        continue
    directory=metrics.parent
    row['vcf_calls']=[]
    for vcf in (directory/'results').glob('*.vcf'):
        for line in vcf.read_text().splitlines():
            if line.startswith('#'): continue
            fields=line.split('\t')
            info=dict(item.split('=',1) for item in fields[7].split(';') if '=' in item)
            sample=dict(zip(fields[8].split(':'),fields[9].split(':')))
            row['vcf_calls'].append(dict(alias=fields[2],size=int(info['ITD_LEN']),
                af=float(info['AF']),depth=int(info['DP']),
                reads=int(sample['AD'].split(',')[1]),strand=sample['SB']))
    log=(directory/'run.log').read_text()
    discovery=re.search(r'Observed haplotype discovery: (\{.*\})',log)
    if discovery: row['discovery']=ast.literal_eval(discovery.group(1))
    audits=list((directory/'results/flt3_data').glob('*_observed_haplotype_audit.tsv'))
    if audits:
        with audits[0].open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        row['seed_audit']=[{k:r.get(k,'') for k in ['allele_id','net_size','exact_reads',
            'plus','minus','status','peak_alias','source_assignments',
            'failed_consensus_exact_reads','validated_reads','reported','final_status',
            'sequence_alternatives']} for r in rows if r['status'] in
            ('added_for_competition','existing_candidate','candidate_cap')]
        counts=Counter()
        for validation in (directory/'results/flt3_data').glob('*_validation_read_support.tsv'):
            with validation.open() as handle:
                for r in csv.DictReader(handle,delimiter='\t'):
                    if r['valid_support']=='True' and r['support_call']=='ITD-supporting' and r['filter_reason']=='pass':
                        counts[r['ref_alias']]+=1
        row['all_validated_candidate_reads']=dict(counts)
    results.append(row)

def classified_reads(sample, mode):
    paths=list((folder/sample/mode/'results/flt3_data').glob('*_validation_read_support.tsv'))
    if not paths: return {}
    with paths[0].open() as handle:
        return {r['read_id']: r['ref_alias'] if r['valid_support']=='True' and
                r['filter_reason']=='pass' and r['support_call'] in
                ('ITD-supporting','WT-supporting') else 'unclassified_or_filtered'
                for r in csv.DictReader(handle,delimiter='\t')}

for result in results:
    if result['mode']!='enabled': continue
    before=classified_reads(result['sample'],'original')
    after=classified_reads(result['sample'],'enabled')
    flow=Counter((before.get(rid,'absent'),call) for rid,call in after.items())
    result['read_flow_from_original']={f'{a} -> {b}':n for (a,b),n in sorted(flow.items())}
target=root/'review/observed_rescue_results.json'
target.write_text(json.dumps(dict(baseline_commit='080b9a7',threads=3,
    container_image='teoloup/hematology-aml-flt3-itd:2.3',
    container_image_id='sha256:c0f5c73116d5b67e9ea9056db8b6c7fc5e07880780a468b48d723b362a980ecd',
    report_af=.005,min_itd_size=12,memory_measure='Largest waited child process ru_maxrss; not aggregate concurrent RSS',
    runs=results),indent=2))
print(target)
for r in results:
    print(r['sample'],r['mode'],round(r['wall_seconds'],2),
          round(r['largest_child_maxrss_kib']/1024,1),'MiB',
          [(c['size'],c['reads'],c['af']) for c in r['vcf_calls']])
