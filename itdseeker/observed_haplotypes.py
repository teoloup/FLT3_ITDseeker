"""Conservative GMM-independent discovery from recurrent complete observations.

This is candidate discovery, not AF estimation: all reads still compete against
WT and every accepted candidate in the normal validation stage. Exact outer
anchors deliberately trade discovery sensitivity for specificity. Only recurrent
windows are tested as whole insertions; isolated gaps never become candidates.
"""
from collections import defaultdict
from collections import Counter
from hashlib import sha256
from statistics import median
from time import perf_counter

import pandas as pd


def discover_observed_haplotypes(reads, ref_seq, *, min_size=12, max_size=300,
                                 min_reads=20):
    ref = str(ref_seq).upper()
    start, stop, anchor = 25, len(ref) - 21, 15
    left, right = ref[start:start + anchor], ref[stop-anchor:stop]
    windows = defaultdict(list)
    began = perf_counter()
    anchored = 0
    for rid, read in reads.items():
        seq = str(read.get('seq', '')).upper()
        x, y = seq.find(left), seq.rfind(right)
        if x < 0 or y < x + anchor:
            continue
        anchored += 1
        window = seq[x:y + anchor]
        net = len(window) - (stop - start)
        if 0 < net <= max_size:
            windows[window].append((rid, read))

    candidates, audit = [], []
    for window, members in sorted(windows.items(), key=lambda x: (-len(x[1]), x[0])):
        allele = ref[:start] + window + ref[stop:]
        net = len(allele) - len(ref)
        plus = sum(r.get('strand', '+') == '+' for _, r in members)
        qualities = [sum(ord(q)-33 for q in r.get('qual', '')) / len(r['qual'])
                     for _, r in members if r.get('qual')]
        row = dict(allele_id=sha256(allele.encode()).hexdigest()[:12],
                   net_size=net, exact_reads=len(members), plus=plus,
                   minus=len(members)-plus, median_q=median(qualities) if qualities else 0,
                   status='', allele=allele)
        reason = None
        if net < min_size:
            reason = 'below_min_size'
        elif set(window) - set('ACGT'):
            reason = 'unresolved_bases'
        elif len(members) < min_reads:
            reason = 'insufficient_recurrence'
        elif min(plus, len(members)-plus) < 2:
            reason = 'insufficient_strand_support'
        elif row['median_q'] < 20:
            reason = 'low_quality'
        if reason is None:
            # A pure insertion must explain the COMPLETE observed net increase.
            # Reconstructing it must reproduce the full observed allele exactly:
            # no silent deletions, substitutions, or separate short error gaps.
            # Check exact decompositions instead of an alignment's preferred
            # gaps: scoring can favor mismatches even for a valid pure insertion
            # when the inserted copy differs by one junction base.
            valid = [p for p in range(start+anchor,stop-anchor+1)
                     if ref[:p] + allele[p:p+net] + ref[p:] == allele]
            if not valid:
                reason = 'not_single_complete_insertion'
            else:
                def duplication_identity(p):
                    ins = allele[p:p+net]
                    adjacent = [ref[max(0,p-net):p], ref[p:p+net]]
                    return max((sum(a == b for a,b in zip(ins,s))/net
                                for s in adjacent if len(s) == net), default=0)
                p = max(valid, key=lambda pos: (duplication_identity(pos), -pos))
                ins = allele[p:p+net]
                identity = duplication_identity(p)
                row['duplication_identity'] = identity
                if identity < .8:
                    reason = 'weak_tandem_duplication'
                else:
                    row.update(insertion_pos=p, insertion_seq=ins)
                    candidates.append(dict(row, members=members))
        row['status'] = reason or 'eligible'
        audit.append(row)
    stats = dict(total_reads=len(reads), anchored_reads=anchored,
                 unique_positive_windows=len(windows),
                 tested_seeds=sum(r['status'] not in ('below_min_size',
                     'insufficient_recurrence','insufficient_strand_support','low_quality',
                     'unresolved_bases') for r in audit),
                 seconds=perf_counter()-began)
    return candidates, audit, stats


def add_observed_candidates(df_cons, candidates, audit, ref_seq, max_candidates=6,
                            reads_df=None, unresolved_aliases=()):
    """Deduplicate full alleles, retain strongest seeds within the candidate cap."""
    ref = str(ref_seq)
    rows = df_cons.to_dict('records') if df_cons is not None else []
    known = {ref[:int(r['consensus_ins_pos_ref'])] + r['consensus_seq'] +
             ref[int(r['consensus_ins_pos_ref']):]: r['peak_alias'] for r in rows}
    by_id = {r['allele_id']: r for r in audit}
    assignments = (dict(zip(reads_df['read_id'], reads_df['gmm_peak_alias']))
                   if reads_df is not None else {})
    for seed in candidates:
        log = by_id[seed['allele_id']]
        sources = Counter(str(assignments.get(rid)) for rid,_ in seed['members'])
        log['source_assignments'] = ';'.join(f'{k}:{v}' for k,v in sorted(sources.items()))
        log['failed_consensus_exact_reads'] = sum(
            assignments.get(rid) in unresolved_aliases for rid,_ in seed['members'])
        alternatives = [s['allele_id'] for s in candidates if s is not seed and
                        len(s['allele']) == len(seed['allele']) and
                        sum(a != b for a,b in zip(s['allele'],seed['allele'])) <= 1]
        note = ('Observed junction sequence alternatives; distinct biological clones are not established.'
                if alternatives else 'Recurrent observed full haplotype; GMM-independent candidate.')
        if seed['allele'] in known:
            log.update(status='existing_candidate', peak_alias=known[seed['allele']])
            if alternatives:
                for row in rows:
                    if row['peak_alias'] == known[seed['allele']]:
                        row['discovery_note'] = note
                log.update(sequence_alternatives=','.join(alternatives), note=note)
            continue
        if len(rows) >= max_candidates:
            log['status'] = 'candidate_cap'
            continue
        alias = 'ITD_OBS_' + seed['allele_id']
        # Same-size near-identical observations are sequence alternatives, not
        # evidence establishing separate biological clones.
        rows.append(dict(peak_alias=alias, consensus_seq=seed['insertion_seq'],
                         consensus_len=seed['net_size'],
                         consensus_ins_pos_ref=seed['insertion_pos'],
                         median_ins_pos_ref=seed['insertion_pos'],
                         n_total_reads=seed['exact_reads'], n_unique=1,
                         max_minor_fraction=0., discovery_note=note))
        known[seed['allele']] = alias
        log.update(status='added_for_competition', peak_alias=alias,
                   sequence_alternatives=','.join(alternatives), note=note)
    return pd.DataFrame(rows)


def annotate_validation(audit, validation, reported_aliases):
    """Persist final support even when every candidate fails the reporting AF."""
    valid = validation.loc[(validation['valid_support'] == True) &
                           (validation['filter_reason'] == 'pass') &
                           (validation['support_call'] == 'ITD-supporting')]
    counts = valid.groupby('ref_alias')['read_id'].nunique().to_dict()
    for row in audit:
        alias = row.get('peak_alias')
        if alias:
            row['validated_reads'] = int(counts.get(alias,0))
            row['reported'] = alias in reported_aliases
            row['final_status'] = ('reported' if row['reported'] else
                'below_reporting_af' if row['validated_reads'] else 'no_validated_support')
