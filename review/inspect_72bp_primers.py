from pathlib import Path
import json,sys
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from Bio import SeqIO
from Bio.Seq import Seq
from GMM_peaks import fit_gmm_itds
folder=ROOT/'review/runs/amplicon_sorter/inputs/13697_2runs_hg38_RG'
reads={r.id:r for r in SeqIO.parse(folder/'trimmed.fastq','fastq')}
fit=fit_gmm_itds({k:str(v.seq) for k,v in reads.items()},min_gmm_fraction=.01,max_itds_detected=6,min_ggmm_peak_distance=10,max_peak_sd=5,assign_width_factor=2,assign_mode='manual',prob_threshold=.85)
print(fit.comps[['peak_alias','mean_bp','effective_read_count']].to_string(index=False))
raw={r.id:str(r.seq) for r in SeqIO.parse(folder/'region_reads.fastq','fastq')}
base=ROOT/'review/runs/amplicon_sorter/13697_2runs_hg38_RG/sc99__native'
for c in json.loads((base/'summary.json').read_text())['clusters']:
    if c['length']<400:continue
    stats=Counter();prefixes=Counter()
    for r in SeqIO.parse(base/c['file'],'fastq'):
        trimmed=str(reads[r.id].seq)
        if not trimmed.startswith('TACCTTTCAGCA'):continue
        original=raw[r.id]
        if trimmed not in original:original=str(Seq(original).reverse_complement())
        pos=original.find(trimmed);assert pos>=0
        stats['n_short_start']+=1;stats['at_raw_start' if pos==0 else 'cut_inside_raw']+=1
        prefixes[original[max(0,pos-10):pos]]+=1
    print(c['file'],stats,prefixes.most_common(10))
