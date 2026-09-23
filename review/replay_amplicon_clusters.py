"""Replay the measured amplicon_sorter clusters through the unchanged pipeline.

Review-only adapter: verifies exact input read membership before using cached
assignments. Does not add a production backend or rerun amplicon_sorter.
"""
from pathlib import Path
import json
import runpy
import sys
from Bio import SeqIO
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import haplotype_split

def replay(**kw):
    folder=ROOT/'review/runs/amplicon_sorter/sim_A/sc995__ITD_2__native'
    summary=json.loads((folder/'summary.json').read_text())
    assert summary['completed'] and summary['exit_code']==0
    expected={r.id for r in SeqIO.parse(ROOT/'review/runs/amplicon_sorter/inputs/sim_A/ITD_2.fastq','fastq')}
    assert set(kw['peak_subsets']['ITD_2'])==expected, 'Cached clustering input no longer matches'
    mapping={}
    for i,cluster in enumerate(summary['clusters']):
        for rec in SeqIO.parse(folder/cluster['file'],'fastq'):
            assert rec.id not in mapping
            mapping[rec.id]=str(i)
    assert set(mapping)==expected
    return haplotype_split.assignments_to_result(
        comps=kw['comps'],reads_df=kw['reads_df'],peak_subsets=kw['peak_subsets'],
        assignments={'ITD_2':mapping},wt_amplicon_length=kw['wt_amplicon_length'],
        wt_peak_tolerance=kw['wt_peak_tolerance'],min_child_fraction=kw['min_child_fraction'],
        min_child_reads=kw['min_child_reads'])

if __name__=='__main__':
    haplotype_split.BACKENDS['amplicon_review']=replay
    out=ROOT/'review/runs/amplicon_sorter/pipeline_replay'
    out.mkdir(parents=True,exist_ok=True)
    sys.argv=[str(ROOT/'Nano_ITDseeker.py'),'-b',str(ROOT/'review/synthetic/sim_A.bam'),
              '-o',str(out),'-s','sim_A','-t','4','--min-allele-frequency','0.01',
              '--haplotype-method','amplicon_review']
    print('Review-only replay of verified cached amplicon_sorter read assignments.',flush=True)
    runpy.run_path(sys.argv[0],run_name='__main__')
