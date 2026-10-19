"""Compare stable discovery and experimental flank settings on retained FASTQs.

No GMM fitting, denoising or read alignment is run. Export the stable module with:
git archive --format=tar --output=review/flexible_anchor_baseline.tar fb9cae3 itdseeker/observed_haplotypes.py
"""
import csv
import json
import sys
import tarfile
import types
from collections import Counter
from pathlib import Path

from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.simulate_itd_data import DEFAULT_REF_WT as WT
from itdseeker.observed_haplotypes import discover_observed_haplotypes

stable=types.ModuleType('stable_observed_haplotypes')
with tarfile.open(ROOT/'review/flexible_anchor_baseline.tar') as archive:
    source=archive.extractfile('itdseeker/observed_haplotypes.py').read().decode()
exec(compile(source,'stable_observed_haplotypes.py','exec'),stable.__dict__)

rows=[]
summary=[]
for sample in ['11528','13462','11551','10808','12156']:
    fastq=ROOT/'FRAGMENT_DATA/recheck_20261019'/sample/'fastq/trimmed.fastq'
    reads={r.id:dict(seq=str(r.seq),qual=''.join(chr(q+33) for q in r.letter_annotations['phred_quality']),
                    strand='-' if r.description.endswith(' rc') else '+')
           for r in SeqIO.parse(str(fastq),'fastq') if len(r.seq)>=330}
    legacy_ids=set()
    for name,fn,config in [
        ('stable_15_exact_identity80',stable.discover_observed_haplotypes,{}),
        ('experimental_15_exact',discover_observed_haplotypes,
         dict(anchor_length=15,anchor_shift=0,anchor_max_mismatches=0)),
        ('experimental_20_shift10_mm3',discover_observed_haplotypes,
         dict(anchor_length=20,anchor_shift=10,anchor_max_mismatches=3)),
        ('experimental_20_shift15_mm3',discover_observed_haplotypes,
         dict(anchor_length=20,anchor_shift=15,anchor_max_mismatches=3)),
    ]:
        seeds,audit,stats=fn(reads,WT,min_size=3 if sample=='11551' else 12,**config)
        ids={s['allele_id'] for s in seeds}
        if name.startswith('stable_'): legacy_ids=ids
        seed_rows=[{key:s.get(key) for key in ['allele_id','net_size','exact_reads',
                   'exact_interior_reads','exact_flank_reads','normalized_flank_reads',
                   'left_anchor_substitutions','right_anchor_substitutions',
                   'plus','minus','median_q','duplication_identity']} for s in seeds]
        row=dict(sample=sample,configuration=name,stats=stats,seeds=seed_rows,
                 rejection_counts=dict(Counter(r['status'] for r in audit)),
                 new_seed_ids=sorted(ids-legacy_ids),lost_stable_seed_ids=sorted(legacy_ids-ids))
        rows.append(row)
        brief=dict(sample=sample,configuration=name,total_reads=len(reads),
                   anchored_reads=stats['anchored_reads'],coverage=stats['anchored_reads']/len(reads),
                   ambiguous_pairs=(stats['anchor_outcomes'].get('ambiguous_anchor_pair',0)
                                    if 'anchor_outcomes' in stats else None),
                   normalized_flank_reads=stats.get('normalized_flank_reads',0),
                   seeds=len(seeds),new_seeds=len(ids-legacy_ids),lost_stable_seeds=len(legacy_ids-ids),
                   seconds=stats['seconds'])
        summary.append(brief)
        print(sample,name,'coverage',round(100*brief['coverage'],2),
              'ambiguous',brief['ambiguous_pairs'],'seconds',round(brief['seconds'],3),
              'seeds',[(s['net_size'],s['exact_reads'],round(s['duplication_identity'],3)) for s in seeds],
              'lost',row['lost_stable_seed_ids'],flush=True)

target=ROOT/'review/flexible_anchor_comparison.json'
target.write_text(json.dumps(dict(stable_commit='fb9cae3',
    note='Stable identity >=0.8 gate; experiments annotate identity only. All experiments require a unique ordered anchor pair. 11551 minimum size=3, other samples=12.',
    runs=rows),indent=2))
with (ROOT/'review/flexible_anchor_comparison.tsv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=summary[0].keys(),delimiter='\t')
    writer.writeheader()
    writer.writerows(summary)
