"""Standalone amplicon_sorter evaluation; never changes the production backend."""
import argparse
from collections import Counter
import csv
import json
import os
import re
from pathlib import Path
import signal
import shutil
import tempfile
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Bio import SeqIO
import edlib
from bam_extractor import extract_flt3_reads
from command_audit import configure_command_log, record_command
from simulate_itd_data import DEFAULT_REF_WT as WT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True, choices=['sim_A', 'sim_no_wt', '13697_2runs_hg38_RG', '14417_2runs_hg38_RG'])
    parser.add_argument('--profile', choices=['default', 'sc99', 'sc995'], default='default')
    parser.add_argument('--length-diff-consensus', type=float, default=None)
    parser.add_argument('--replicate', type=int, default=1)
    parser.add_argument('--maxreads', type=int, default=100000)
    parser.add_argument('--native', action='store_true', help='Use native Linux temp storage for upstream intermediate files')
    parser.add_argument('--diagnostic', action='store_true', help='Expose upstream swallowed exceptions; no clustering changes')
    parser.add_argument('--peak', help='Optional first-pass GMM peak alias, for a direct per-peak comparison')
    parser.add_argument('--timeout', type=int, default=600)
    args = parser.parse_args()
    if args.maxreads < 1 or args.timeout < 1: parser.error('maxreads and timeout must be positive')
    base = ROOT / 'review/runs/amplicon_sorter'
    inputs = base / 'inputs' / args.sample
    inputs.mkdir(parents=True, exist_ok=True)
    fq = inputs / 'trimmed.fastq'
    configure_command_log(inputs / 'commands.jsonl')
    if not fq.exists():
        bam = ROOT / ('review/synthetic' if args.sample.startswith('sim_') else 'bam_data/test_bam') / (args.sample + '.bam')
        extract_flt3_reads(str(bam), 'hg38', 4, 20, 330, 300, 336, str(inputs))
    if args.peak:
        from GMM_peaks import fit_gmm_itds
        from haplotype_split import write_peak_fastq
        reads = {r.id: dict(seq=str(r.seq), strand='+',
                 qual=''.join(chr(q+33) for q in r.letter_annotations['phred_quality']))
                 for r in SeqIO.parse(fq, 'fastq')}
        fit = fit_gmm_itds(reads, min_gmm_fraction=.01, max_itds_detected=6,
                          min_ggmm_peak_distance=10, max_peak_sd=5, assign_width_factor=2,
                          assign_mode='manual', prob_threshold=.85)
        if args.peak not in fit.peak_subsets: raise ValueError(args.peak)
        fq = inputs / (args.peak + '.fastq')
        write_peak_fastq(fit.reads_df, fit.peak_subsets[args.peak], str(fq))
        fit.comps.to_csv(inputs / 'gmm_peaks.tsv', sep='\t', index=False)
    label = args.profile + ('__' + args.peak if args.peak else '')
    if args.diagnostic: label += '__diagnostic'
    if args.native: label += '__native'
    if args.maxreads != 100000: label += '__cap' + str(args.maxreads)
    if args.length_diff_consensus is not None: label += '__ldc' + format(args.length_diff_consensus, 'g')
    if args.replicate != 1: label += '__rep' + str(args.replicate)
    out = base / args.sample / label
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'summary.json').exists():
        print((out / 'summary.json').read_text()); return
    configure_command_log(out / 'commands.jsonl')
    work_out = Path(tempfile.mkdtemp(prefix='flt3-amplicon-review-')) if args.native else out
    tool = ROOT / 'review/tools/amplicon_sorter' / ('amplicon_sorter_diagnostic.py' if args.diagnostic else 'amplicon_sorter.py')
    if args.diagnostic:
        source = (tool.parent / 'amplicon_sorter.py').read_text(encoding='utf-8')
        old = '            except Exception: \n                continue'
        if old not in source: raise RuntimeError('Pinned upstream exception handler changed')
        source = source.replace(old, '            except Exception:\n                import traceback\n                traceback.print_exc()\n                raise')
        tool.write_text(source, encoding='utf-8')
    cmd = [sys.executable, '-u', str(tool),
           '-i', str(fq), '--allreads', '--maxreads', str(args.maxreads), '-np', '4',
           '--minlength', '330', '--maxlength', '666', '--save_fastq', '-o', str(work_out)]
    if args.profile != 'default':
        cmd += ['--similar_consensus', '99' if args.profile == 'sc99' else '99.5']
    if args.length_diff_consensus is not None:
        cmd += ['--length_diff_consensus', str(args.length_diff_consensus)]
    record_command(cmd, stdin='inherited', stdout=str(out / 'console.log'), stderr='merged-with-stdout')
    start = time.monotonic()
    with (out / 'console.log').open('w') as handle:
        proc = subprocess.Popen(cmd, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = proc.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGTERM)
            try: proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL); proc.wait()
            code = 'timeout'
    if args.native:
        shutil.copytree(work_out, out / 'native_results')
        shutil.rmtree(work_out)
    remaining_groups = list(out.rglob('*.group'))
    completed = code == 0 and not remaining_groups and bool(list(out.rglob('*_consensussequences.fasta')))
    input_ids = {rec.id for rec in SeqIO.parse(fq, 'fastq')}
    truth_reads, truth_alleles = {}, {'WT': WT}
    if args.sample.startswith('sim_'):
        scenario = args.sample.removeprefix('sim_')
        with (ROOT / f'review/synthetic/truth_reads_{scenario}.tsv').open() as handle:
            truth_reads = {r['read_id']: r['haplotype'] for r in csv.DictReader(handle, delimiter='\t')}
        with (ROOT / f'review/synthetic/truth_haplotypes_{scenario}.tsv').open() as handle:
            for r in csv.DictReader(handle, delimiter='\t'):
                if int(r['itd_len']):
                    p = int(r['ins_pos_local'])
                    truth_alleles[r['haplotype']] = WT[:p] + r['itd_seq'] + WT[p:]
    clusters, assigned, duplicate = [], set(), set()
    for file in sorted(out.rglob('*.fastq')):
        if 'unique' in file.stem: continue
        records = list(SeqIO.parse(file, 'fastq'))
        ids = {r.id for r in records}
        duplicate |= ids & assigned
        assigned |= ids
        fa = file.with_suffix('.fasta')
        consensus = next((str(r.seq) for r in SeqIO.parse(fa, 'fasta') if r.id == 'consensus'), '') if fa.exists() else ''
        distances = {name: edlib.align(consensus, seq, mode='NW')['editDistance'] for name, seq in truth_alleles.items()} if consensus else {}
        clusters.append(dict(file=str(file.relative_to(out)), n_reads=len(records), unique_ids=len(ids),
                             truth_counts=dict(Counter(truth_reads[i] for i in ids if i in truth_reads)),
                             consensus=consensus, length=len(consensus), truth_edit_distances=distances))
    completed = completed and not duplicate and not (assigned-input_ids) and all(
        c['consensus'] and c['n_reads']==c['unique_ids'] for c in clusters)
    selection = re.search(r'(\d+) out of (\d+) sequences', (out / 'console.log').read_text())
    selected = int(selection.group(1)) if selection else None
    row = dict(sample=args.sample, profile=args.profile, peak=args.peak, native=args.native, maxreads=args.maxreads, exit_code=code,
               length_diff_consensus=args.length_diff_consensus, replicate=args.replicate,
               selected_reads=selected,
               unassigned_selected_reads=selected-len(assigned) if selected is not None else None,
               completed=completed, remaining_groups=[str(p.relative_to(out)) for p in remaining_groups],
               seconds=round(time.monotonic()-start, 2), input_reads=len(input_ids),
               assigned_reads=len(assigned), unassigned_reads=len(input_ids-assigned),
               duplicate_assignments=len(duplicate), unknown_ids=len(assigned-input_ids),
               clusters=clusters)
    (out / 'summary.json').write_text(json.dumps(row, indent=2) + '\n')
    print(json.dumps(row, indent=2), flush=True)

if __name__ == '__main__':
    main()
