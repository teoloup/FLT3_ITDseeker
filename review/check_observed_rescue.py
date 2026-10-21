"""Replay inexpensive discovery on audit FASTQs without committing read data."""
import json
from pathlib import Path
from Bio import SeqIO
from scripts.simulate_itd_data import DEFAULT_REF_WT as WT
from itdseeker.observed_haplotypes import discover_observed_haplotypes

root=Path('/workspace')
results=[]
for sample in ['11528','13462','11551','10808','12156']:
    fastq=root/'FRAGMENT_DATA/recheck_20261019'/sample/'fastq/trimmed.fastq'
    reads={r.id:dict(seq=str(r.seq),qual=''.join(chr(q+33) for q in r.letter_annotations['phred_quality']),
                    strand='-' if r.description.endswith(' rc') else '+')
           for r in SeqIO.parse(str(fastq),'fastq') if len(r.seq)>=330}
    seeds,audit,stats=discover_observed_haplotypes(reads,WT,min_size=3 if sample=='11551' else 12)
    result=dict(sample=sample,stats=stats,seeds=[{k:s[k] for k in
        ['allele_id','net_size','exact_reads','plus','minus','median_q']} for s in seeds])
    if sample=='11551':
        assert not any(s['net_size']==4 for s in seeds), 'Spurious isolated 4bp seed'
    results.append(result)
target=root/'review/observed_rescue_discovery.json'
target.write_text(json.dumps(results,indent=2))
print(target.read_text())
