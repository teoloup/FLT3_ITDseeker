"""Validate all read-derived target candidates competitively against WT/known ITDs."""
import sys, json, logging
from pathlib import Path
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from Bio import SeqIO
sys.path.insert(0,'/workspace')
from itdseeker.Pairwise_aligment_toolkit import align_reads_multi_ref_parallel
from itdseeker.Helper_functions import calculate_allele_frequencies_and_strand_bias
from fragment_recheck_scan import ROOT,OUT,WT,START

def validate(entry):
    sample=entry['biological_sample']; p=OUT/sample
    logging.basicConfig(level=logging.WARNING); log=logging.getLogger('recheck')
    if (p/'validation_summary.tsv').exists(): return sample,pd.read_csv(p/'validation_summary.tsv',sep='\t').to_dict('records')
    groups=pd.read_csv(p/'allele_groups.tsv',sep='\t')
    chosen={}; target_matches={}
    for target in entry['targets']:
        size=target['Fragment size (bp)']; window=1.75 if size<12 else 2.5
        close=groups[abs(groups.insertion_length-size)<=window]
        top=int(close.n_reads.max()) if len(close) else 0
        close=close[close.n_reads>=max(5,top*.15)].head(3)
        target_matches[str(target['excel_row'])]=list(close.allele_id)
        for row in close.to_dict('records'): chosen[row['allele_id']]=row
    for known in entry['known_calls']:
        close=groups[abs(groups.insertion_length-known['Nanopore size (bp)'])<=1].head(1)
        for row in close.to_dict('records'): chosen[row['allele_id']]=row
    chosen=list(chosen.values())
    for i,row in enumerate(chosen,1): row['ref_alias']=f'ITD_C{i}'
    cd=pd.DataFrame(chosen); cd.to_csv(p/'candidate_references.tsv',sep='\t',index=False)
    (p/'target_candidates.json').write_text(json.dumps(target_matches,indent=2))
    refs={'WT':{'ref_seq_with_itd':WT}}
    cons=[]
    with open(p/'candidate_references.fasta','w') as f:
        f.write('>WT\n'+WT+'\n')
        for row in chosen:
            refs[row['ref_alias']]={'ref_seq_with_itd':row['allele']}
            cons.append(dict(peak_alias=row['ref_alias'],consensus_ins_pos_ref=row['insertion_pos'],consensus_len=row['insertion_length']))
            f.write('>'+row['ref_alias']+'\n'+row['allele']+'\n')
    records=[]; sequence_id={}; unique=[]
    for rec in SeqIO.parse(str(p/'fastq/trimmed.fastq'),'fastq'):
        seq=str(rec.seq)
        if seq not in sequence_id:
            uid=f'U{len(sequence_id)}'; sequence_id[seq]=uid
            unique.append(dict(read_id=uid,read_seq=seq,strand='+'))
        records.append(dict(read_id=rec.id,unique_id=sequence_id[seq],strand='-' if rec.description.endswith(' rc') else '+',read_len=len(seq),mean_q=sum(rec.letter_annotations['phred_quality'])/len(seq)))
    # Identical sequences have identical competitive alignments. Expand afterward
    # to restore every read ID and its original orientation, without deduplicating AF.
    val=align_reads_multi_ref_parallel(pd.DataFrame(unique),refs,pd.DataFrame(cons),logger=log,threads=2)
    val=val.drop(columns=['strand']).rename(columns={'read_id':'unique_id'})
    expanded=pd.DataFrame(records).merge(val,on='unique_id',validate='many_to_one')
    expanded.to_csv(p/'competitive_read_support.tsv.gz',sep='\t',index=False)
    summary=calculate_allele_frequencies_and_strand_bias(expanded,min_allele_frequency=0)
    summary=cd.merge(summary,on='ref_alias',how='left')
    summary['input_trimmed_reads']=len(expanded)
    summary['support_fraction_all_trimmed']=summary.n_itd_reads.fillna(0)/len(expanded)
    summary['q20_validated_support']=summary.ref_alias.map(expanded[(expanded.valid_support)&(expanded.support_call=='ITD-supporting')&(expanded.mean_q>=20)].ref_alias.value_counts()).fillna(0).astype(int)
    summary.to_csv(p/'validation_summary.tsv',sep='\t',index=False)
    (p/'validation_counts.json').write_text(json.dumps(dict(n_reads=len(expanded),unique_sequences=len(unique),references=len(refs),support_calls=dict(Counter(expanded.support_call)),filter_reasons=dict(Counter(expanded.filter_reason))),indent=2))
    return sample,summary[['insertion_length','n_itd_reads','allele_frequency','plus_reads','minus_reads']].fillna(0).to_dict('records')

def main():
    manifest=json.loads((OUT/'manifest.json').read_text())
    entries=[e for e in manifest['samples'] if len(sys.argv)==1 or e['biological_sample'] in sys.argv[1:]]
    with ProcessPoolExecutor(max_workers=3) as pool:
        jobs=[pool.submit(validate,e) for e in entries]
        for future in as_completed(jobs):
            sample,result=future.result(); print('VALIDATED',sample,json.dumps(result),flush=True)

if __name__=='__main__': main()
