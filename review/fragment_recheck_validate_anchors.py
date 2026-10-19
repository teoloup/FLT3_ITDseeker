"""Competitive confirmation of complete recurrent haplotypes with measured net size."""
import sys,json,logging
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
from collections import Counter
import pandas as pd
from Bio import SeqIO
sys.path.insert(0,'/workspace')
from itdseeker.Pairwise_aligment_toolkit import align_reads_multi_ref_parallel
from itdseeker.Helper_functions import calculate_allele_frequencies_and_strand_bias
from fragment_recheck_scan import OUT,WT
CLEAN='--clean' in sys.argv
PREFIX='clean_' if CLEAN else 'anchor_'

def run(entry):
    p=OUT/entry['biological_sample']; logging.basicConfig(level=logging.WARNING)
    if (p/(PREFIX+'validation_summary.tsv')).exists(): return entry['biological_sample']
    groups=pd.read_csv(p/'anchored_haplotypes.tsv',sep='\t')
    selected={}; targetids={}
    for t in entry['targets']:
        close=groups[abs(groups.net_insertion_length-t['Fragment size (bp)'])<=2.5]
        top=int(close.n_reads.max()) if len(close) else 0
        close=close[close.n_reads>=max(5,top*.2)].head(3)
        targetids[str(t['excel_row'])]=list(close.allele_id)
        for row in close.to_dict('records'): selected[row['allele_id']]=row
    for known in entry['known_calls']:
        close=groups[abs(groups.net_insertion_length-known['Nanopore size (bp)'])<=.1].head(1)
        for row in close.to_dict('records'): selected[row['allele_id']]=row
    # Strong unsolicited haplotypes avoid allocating their reads to a target.
    if not CLEAN:
        for row in groups[groups.n_reads>=50].head(12).to_dict('records'):
            selected[row['allele_id']]=row
    selected=list(selected.values()); refs={'WT':{'ref_seq_with_itd':WT}}; cons=[]
    if not selected:
        columns=list(groups.columns)+['ref_alias','n_itd_reads','n_total_reads','allele_frequency','plus_reads','minus_reads','input_trimmed_reads','support_fraction_all_trimmed','q20_validated_support']
        pd.DataFrame(columns=columns).to_csv(p/(PREFIX+'validation_summary.tsv'),sep='\t',index=False)
        (p/(PREFIX+'target_candidates.json')).write_text(json.dumps(targetids,indent=2))
        (p/(PREFIX+'validation_counts.json')).write_text(json.dumps(dict(status='No read-supported candidate or known ITD; no artificial variant reference was generated.')))
        return entry['biological_sample']
    for i,row in enumerate(selected,1):
        alias=f'ITD_A{i}'; row['ref_alias']=alias
        refs[alias]={'ref_seq_with_itd':row['allele']}
        cons.append(dict(peak_alias=alias,consensus_ins_pos_ref=row['insertion_pos'],consensus_len=row['insertion_length']))
    cd=pd.DataFrame(selected); cd.to_csv(p/(PREFIX+'candidate_references.tsv'),sep='\t',index=False)
    (p/(PREFIX+'target_candidates.json')).write_text(json.dumps(targetids,indent=2))
    with open(p/(PREFIX+'candidate_references.fasta'),'w') as f:
        for alias,entry2 in refs.items(): f.write('>'+alias+'\n'+entry2['ref_seq_with_itd']+'\n')
    unique=[]; sequence_id={}; records=[]
    for rec in SeqIO.parse(str(p/'fastq/trimmed.fastq'),'fastq'):
        seq=str(rec.seq)
        if seq not in sequence_id:
            uid=f'U{len(sequence_id)}'; sequence_id[seq]=uid
            unique.append(dict(read_id=uid,read_seq=seq,strand='+'))
        records.append(dict(read_id=rec.id,unique_id=sequence_id[seq],strand='-' if rec.description.endswith(' rc') else '+',mean_q=sum(rec.letter_annotations['phred_quality'])/len(seq)))
    val=align_reads_multi_ref_parallel(pd.DataFrame(unique),refs,pd.DataFrame(cons),logger=logging.getLogger('anchor_validation'),threads=2)
    val=val.drop(columns=['strand']).rename(columns={'read_id':'unique_id'})
    expanded=pd.DataFrame(records).merge(val,on='unique_id',validate='many_to_one')
    expanded.to_csv(p/(PREFIX+'competitive_read_support.tsv.gz'),sep='\t',index=False)
    summary=calculate_allele_frequencies_and_strand_bias(expanded,min_allele_frequency=0)
    summary=cd.merge(summary,on='ref_alias',how='left')
    summary['input_trimmed_reads']=len(expanded)
    summary['support_fraction_all_trimmed']=summary.n_itd_reads.fillna(0)/len(expanded)
    summary['q20_validated_support']=summary.ref_alias.map(expanded[(expanded.valid_support)&(expanded.support_call=='ITD-supporting')&(expanded.mean_q>=20)].ref_alias.value_counts()).fillna(0).astype(int)
    summary.to_csv(p/(PREFIX+'validation_summary.tsv'),sep='\t',index=False)
    (p/(PREFIX+'validation_counts.json')).write_text(json.dumps(dict(reads=len(expanded),unique_sequences=len(unique),references=len(refs),filter_reasons=dict(Counter(expanded.filter_reason)),support_calls=dict(Counter(expanded.support_call))),indent=2))
    return entry['biological_sample']

if __name__=='__main__':
    manifest=json.loads((OUT/'manifest.json').read_text())
    with ProcessPoolExecutor(max_workers=3) as pool:
        selected_samples=[x for x in sys.argv[1:] if not x.startswith('--')]
        jobs=[pool.submit(run,e) for e in manifest['samples'] if not selected_samples or e['biological_sample'] in selected_samples]
        for f in as_completed(jobs): print('ANCHOR_VALIDATED',f.result(),flush=True)
