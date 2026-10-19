"""Conservative GMM-independent discovery from recurrent complete observations.

This is candidate discovery, not AF estimation: all reads still compete against
WT and every accepted candidate in the normal validation stage. Exact outer
anchors locate an exact observed interior. On this experimental branch bounded
anchor substitutions are normalized to WT and explicitly audited. Only recurrent
interiors are tested as whole insertions; isolated gaps never become candidates.
"""
from collections import defaultdict
from collections import Counter
from hashlib import sha256
from statistics import median
from time import perf_counter

import pandas as pd


def anchor_coordinates(ref_length, anchor_length=20, anchor_shift=10):
    """Shift the original [25, len(WT)-21) outer window toward both ends."""
    start, stop = 25-anchor_shift, ref_length-21+anchor_shift
    if anchor_length < 1 or anchor_shift < 0 or start < 0 or stop > ref_length:
        raise ValueError('Anchor length must be positive and shifted anchors must stay inside WT')
    if start + anchor_length >= stop-anchor_length:
        raise ValueError('Anchor windows must leave a nonempty interior')
    return start, stop


def anchor_matches(seq, anchor, max_mismatches):
    """All substitution-only hits, using k+1 exact seeds and bounded Hamming checks.

    With <=k substitutions, at least one of k+1 disjoint anchor pieces is exact.
    No dynamic alignment or third-party approximate-matching dependency is used.
    """
    length=len(anchor)
    if not 0 <= max_mismatches < length:
        raise ValueError('Anchor mismatches must be >=0 and smaller than anchor length')
    starts=set()
    for part in range(max_mismatches+1):
        lo=part*length//(max_mismatches+1)
        hi=(part+1)*length//(max_mismatches+1)
        seed=anchor[lo:hi]
        offset=seq.find(seed)
        while offset >= 0:
            pos=offset-lo
            if 0 <= pos <= len(seq)-length:
                starts.add(pos)
            offset=seq.find(seed,offset+1)
    hits=[]
    for pos in sorted(starts):
        observed=seq[pos:pos+length]
        if set(observed)-set('ACGT'):
            continue
        mismatches=0
        for actual,expected in zip(observed,anchor):
            mismatches += actual != expected
            if mismatches > max_mismatches:
                break
        if mismatches <= max_mismatches:
            hits.append((pos,mismatches))
    return hits


def locate_anchor_pair(seq, left, right, max_mismatches):
    """Require one ordered pair within the mismatch budget; reject alternatives."""
    left_hits=anchor_matches(seq,left,max_mismatches)
    if not left_hits:
        return None,'missing_left_anchor'
    right_hits=anchor_matches(seq,right,max_mismatches)
    if not right_hits:
        return None,'missing_right_anchor'
    pairs=[(x,y,lm,rm) for x,lm in left_hits for y,rm in right_hits if x+len(left) <= y]
    if not pairs:
        return None,'no_ordered_anchor_pair'
    if len(pairs) != 1:
        return None,'ambiguous_anchor_pair'
    return pairs[0],'anchored'


def discover_observed_haplotypes(reads, ref_seq, *, min_size=12, max_size=300,
                                 min_reads=20, anchor_length=20, anchor_shift=10,
                                 anchor_max_mismatches=3):
    ref = str(ref_seq).upper()
    start,stop=anchor_coordinates(len(ref),anchor_length,anchor_shift)
    anchor=anchor_length
    if not 0 <= anchor_max_mismatches < anchor:
        raise ValueError('Anchor mismatches must be >=0 and smaller than anchor length')
    left, right = ref[start:start + anchor], ref[stop-anchor:stop]
    windows = defaultdict(list)
    flank_changes=defaultdict(list)
    outcomes=Counter()
    locations={}
    began = perf_counter()
    anchored = 0
    normalized_anchors = 0
    for rid, read in reads.items():
        seq = str(read.get('seq', '')).upper()
        if seq not in locations:
            locations[seq]=locate_anchor_pair(seq,left,right,anchor_max_mismatches)
        pair,outcome=locations[seq]
        outcomes[outcome]+=1
        if pair is None:
            continue
        x,y,lm,rm=pair
        anchored += 1
        normalized_anchors += lm+rm > 0
        # The locator anchors are normalized; the entire observed interior is
        # retained EXACTLY. Substitutions/indels inside it are never repaired.
        window = left + seq[x+anchor:y] + right
        net = len(window) - (stop - start)
        if 0 < net <= max_size:
            windows[window].append((rid, read))
            flank_changes[window].append((lm,rm))

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
        changes=flank_changes[window]
        row.update(exact_interior_reads=len(members),
                   exact_flank_reads=sum(lm+rm == 0 for lm,rm in changes),
                   normalized_flank_reads=sum(lm+rm > 0 for lm,rm in changes),
                   left_anchor_substitutions=sum(lm for lm,_ in changes),
                   right_anchor_substitutions=sum(rm for _,rm in changes),
                   max_left_anchor_substitutions=max(lm for lm,_ in changes),
                   max_right_anchor_substitutions=max(rm for _,rm in changes),
                   anchor_length=anchor, anchor_shift=anchor_shift,
                   anchor_max_mismatches=anchor_max_mismatches)
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
                # Experimental: identity describes the adjacent copy match;
                # it no longer decides whether an observed insertion survives.
                row.update(insertion_pos=p, insertion_seq=ins)
                candidates.append(dict(row, members=members))
        row['status'] = reason or 'eligible'
        audit.append(row)
    stats = dict(total_reads=len(reads), anchored_reads=anchored,
                 anchor_outcomes=dict(outcomes), unique_sequences_searched=len(locations),
                 anchor_length=anchor,anchor_shift=anchor_shift,
                 anchor_max_mismatches=anchor_max_mismatches,
                 left_anchor_interval=[start,start+anchor],
                 right_anchor_interval=[stop-anchor,stop],
                 left_anchor=left,right_anchor=right,
                 normalized_flank_reads=normalized_anchors,
                 normalized_positive_window_reads=sum(sum(lm+rm>0 for lm,rm in v) for v in flank_changes.values()),
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
                if alternatives else 'Recurrent observed interior haplotype; GMM-independent candidate.')
        note += (f" Experimental flank locator: {seed['normalized_flank_reads']}/{seed['exact_interior_reads']}"
                 " observations have anchor substitutions replaced by WT; these are not called flank variants."
                 f" Adjacent duplication identity {seed['duplication_identity']:.1%} is annotation only.")
        annotations={key:seed[key] for key in ['duplication_identity','exact_interior_reads',
                     'exact_flank_reads','normalized_flank_reads','anchor_length','anchor_shift',
                     'anchor_max_mismatches','left_anchor_substitutions','right_anchor_substitutions']}
        if seed['allele'] in known:
            log.update(status='existing_candidate', peak_alias=known[seed['allele']])
            for row in rows:
                if row['peak_alias'] == known[seed['allele']]:
                    row.update(annotations, discovery_note=note)
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
                         max_minor_fraction=0., discovery_note=note, **annotations))
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
