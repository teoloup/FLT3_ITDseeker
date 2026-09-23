"""Audit completed Medaka runs using the isolated Medaka Python environment."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pysam


def audit(folder):
    result = json.loads((folder / 'summary.json').read_text(encoding='utf-8'))
    if not result['completed']:
        return None
    output = folder / 'output'
    with pysam.AlignmentFile(str(output / 'calls_to_draft.bam'), 'rb') as bam:
        length = bam.lengths[0]
        depth = np.zeros(length, dtype=int)
        mapped = 0
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            mapped += 1
            for position in read.get_reference_positions():
                depth[position] += 1
    datasets = {}
    inferred = set()
    with h5py.File(output / 'consensus_probs.hdf', 'r') as hdf:
        def visit(name, item):
            if not isinstance(item, h5py.Dataset):
                return
            datasets[name] = list(item.shape)
            if name.endswith('/positions'):
                positions = item[()]
                inferred.update(int(p) for p in positions['major'])
        hdf.visititems(visit)
    gaps = (output / 'consensus.fasta.gaps_in_draft_coords.bed').read_text(encoding='utf-8')
    log = (folder / 'console.log').read_text(encoding='utf-8')
    details = dict(
        case=result['case']['name'], mapped_primary_reads=mapped,
        draft_length=length, zero_depth_bases=int((depth == 0).sum()),
        min_depth=int(depth.min()), median_depth=float(np.median(depth)),
        first_ten_depths=depth[:10].tolist(), last_ten_depths=depth[-10:].tolist(),
        inferred_draft_positions=len(inferred),
        all_draft_positions_inferred=set(range(length)).issubset(inferred),
        gap_bed=gaps, hdf_datasets=datasets,
        inference_finished='Finished processing all regions.' in log,
        ambiguous_bases=sum(b not in 'ACGT' for b in result['sequence'].upper()),
        polished_sha256=hashlib.sha256((output / 'consensus.fasta').read_bytes()).hexdigest(),
    )
    assert mapped > 0 and details['inference_finished'] and inferred, details
    (folder / 'coverage_audit.json').write_text(json.dumps(details, indent=2) + '\n', encoding='utf-8')
    return details


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    args = parser.parse_args()
    rows = [audit(p.parent) for p in sorted(args.run_directory.glob('*/summary.json'))]
    print(json.dumps([r for r in rows if r is not None], indent=2))
