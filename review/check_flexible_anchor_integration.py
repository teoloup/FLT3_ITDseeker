"""Verify one saved experimental integration against the saved stable result."""
import ast
import csv
import json
import re
from pathlib import Path
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
base=ROOT/'review/observed_rescue_benchmark/13462'
stable=base/'enabled/results/flt3_data'
experiment=base/'flexible_shift10/results/flt3_data'

def classifications(folder):
    with (folder/'13462_validation_read_support.tsv').open() as handle:
        return {r['read_id']:tuple(r[k] for k in ['ref_alias','support_call','valid_support','filter_reason'])
                for r in csv.DictReader(handle,delimiter='\t')}

before,after=classifications(stable),classifications(experiment)
references=lambda folder:{r.id:str(r.seq) for r in SeqIO.parse(str(folder/'13462_validation_refs.fasta'),'fasta')}
matching_references=references(stable)==references(experiment)
with (experiment/'13462_observed_haplotype_audit.tsv').open() as handle:
    audit=[{k:r[k] for k in ['peak_alias','net_size','exact_interior_reads','exact_flank_reads',
                            'normalized_flank_reads','duplication_identity','validated_reads','reported']}
           for r in csv.DictReader(handle,delimiter='\t') if r.get('peak_alias')]

calls=[]
for line in (base/'flexible_shift10/results/13462_FLT3_ITD_calls.vcf').read_text().splitlines():
    if line.startswith('#'): continue
    fields=line.split('\t')
    info=dict(v.split('=',1) for v in fields[7].split(';') if '=' in v)
    calls.append(dict(alias=fields[2],info=info))
    assert 'DUP_IDENTITY' in info and 'FLANK_NORMALIZED_READS' in info
log=(base/'flexible_shift10/run.log').read_text()
stats=ast.literal_eval(re.search(r'Observed haplotype discovery: (\{.*\})',log).group(1))
html=(base/'flexible_shift10/results/13462_itd_report.html').read_text()
assert 'anchor substitutions replaced by WT' in html
assert 'annotation only' in html
metrics=json.loads((base/'flexible_shift10/metrics.json').read_text())
result=dict(sample='13462',metrics=metrics,discovery_stats=stats,
            same_reference_sequences=matching_references,
            read_count=len(after),all_read_classifications_unchanged=before==after,
            seed_audit=audit,reported_vcf_calls=calls,
            html_normalization_and_identity_annotations_present=True)
target=ROOT/'review/flexible_anchor_integration.json'
target.write_text(json.dumps(result,indent=2))
assert matching_references and before==after
assert metrics['returncode']==0
print('Integration passed:',len(after),'unchanged read classifications; identical reference sequences;',
      len(calls),'reported VCF calls; HTML/VCF annotations verified.')
