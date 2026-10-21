"""Report-only MSAs, each built exclusively from one final allele's support."""
import json
from pathlib import Path
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from .Multiple_seq_aligment_toolkit import (
    dedup_with_counts, select_unique_panel, run_muscle5_on_pairs,
    weighted_consensus_from_msa, plot_consensus_agreement,
)


def plot_final_support_msas(validation, read_sequences, summary, references,
                            output_dir, sample, *, threads=1, max_unique=100):
    """Align at most 100 distinct full reads per call, weighted by recurrence.

    The called reference locates the displayed insertion window only; it does
    not vote in the plotted read consensus. This never changes any variant.
    Cutadapt-trimmed reads already have reference orientation on both strands.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    passed = validation[validation.valid_support.astype(str).str.lower().eq('true')
                        & validation.filter_reason.eq('pass')
                        & validation.support_call.eq('ITD-supporting')]
    result = {}
    for row in summary.itertuples():
        alias = str(row.ref_alias)
        ids = sorted(passed.loc[passed.ref_alias.eq(alias), 'read_id'].astype(str))
        if len(ids) != int(row.n_itd_reads) or len(ids) != len(set(ids)):
            raise ValueError(f'{alias}: final support membership/count mismatch')
        sequences = [str(read_sequences[rid]).upper() for rid in ids]
        counts = dedup_with_counts(sequences)
        panel = select_unique_panel(counts, max_unique, .95)
        represented = sum(n for _, n in panel)
        info = references[alias]
        reference = str(info['ref_seq_with_itd'])
        position = int(info['local_insertion_pos'])
        size = int(info['itd_length'])
        # Locate a short report window per distinct selected read. Alignment here
        # is diagnostic only; final assignments and reported references stay fixed.
        left, right = max(0, position-12), min(len(reference), position+size+12)
        locator = PairwiseAligner(mode='global',match_score=2,mismatch_score=-2,
                                  open_gap_score=-5,extend_gap_score=-.5)
        context_panel = []
        for sequence, weight in panel:
            alignment = locator.align(reference, sequence)[0]
            target_indices, query_indices = alignment.indices
            selected_indices = query_indices[(target_indices >= left) & (target_indices < right) & (query_indices >= 0)]
            if not len(selected_indices):
                raise ValueError(f'{alias}: supporting read has no aligned insertion context')
            context_panel.append((sequence[int(selected_indices.min()):int(selected_indices.max())+1],weight))
        # A zero-weight reference is used to crop common MSA columns, then removed.
        msa, weights = run_muscle5_on_pairs([(reference[left:right], 0)] + context_panel,
                                            prefix='support', threads=threads)
        ref_index = next(i for i,r in enumerate(msa) if r.id == 'support_1')
        ref_columns = [j for j,c in enumerate(str(msa[ref_index].seq)) if c != '-']
        lo = ref_columns[0]
        hi = ref_columns[-1] + 1
        reads_msa = MultipleSeqAlignment([r[lo:hi] for i,r in enumerate(msa) if i != ref_index])
        cc = weighted_consensus_from_msa(reads_msa, weights)
        path = plot_consensus_agreement(
            cc, alias=alias+'_final_support', sample_name=sample, out_dir=str(output_dir),
            resolved=True, ins_start=None, ins_len=size, base_threshold=.7,
            n_reads=represented, panel_used=len(panel),
            title_note=f'FINAL SUPPORT: {represented}/{len(ids)} reads represented; up to {max_unique} distinct sequences; window positions; insertion with flanks')
        selected = {seq for seq,_ in panel}
        result[alias] = dict(total_support=len(ids), represented_reads=represented,
                             distinct_sequences=len(panel), max_unique=max_unique,
                             read_ids=ids, panel_read_ids=[rid for rid in ids if str(read_sequences[rid]).upper() in selected],
                             plot=str(path))
    (output_dir/f'{sample}_final_support_msa_membership.json').write_text(json.dumps(result,indent=2))
    return result
