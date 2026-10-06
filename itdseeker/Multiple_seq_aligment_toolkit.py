import os
os.environ["MPLBACKEND"] = "Agg"       # disable any GUI backend
os.environ["DISPLAY"] = ""             # make sure Tk can't open a window
os.environ["TK_SILENCE_DEPRECATION"] = "1"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging
import time
import pymuscle5
from .Helper_functions import build_default_aligner, find_insertions
from dataclasses import dataclass
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

# WT bases kept either side of the span where reads place the insertion. Every
# allele carries identical WT outside that span, so it adds nothing to the
# consensus; MUSCLE time grows with sequence length, so it is trimmed off and
# these flanks are kept only to anchor the alignment.
CONTEXT_FLANK_BP = 30

# Share of reads the span must reach. A few reads place the insertion far from
# the rest, which would otherwise stretch the window over most of the amplicon.
# At the real insertion those alleles only carry gaps, which never vote on a
# base, so they are left out of the alignment.
CONTEXT_BOUNDARY_SHARE = 0.95

# Share of weighted support disagreeing with a called base at which a consensus
# counts as mixed. Nano_ITDseeker warns at this level; the plot marks columns past it.
MIXED_MINOR_FRACTION = 0.15

BASE_COLOURS = {"A": "#2e8b3e", "C": "#2a63c4", "G": "#c47f00", "T": "#c8323a"}

def dedup_with_counts(seqs: List[str]) -> List[Tuple[str, int]]:
    # Count sequences case-insensitively but preserve first-seen order for
    # deterministic tie-breaking. Returns list of (seq_upper, count).
    counts: Dict[str, int] = {}
    first_seen: Dict[str, int] = {}
    for i, s in enumerate(seqs):
        if not s:
            continue
        su = s.upper()
        counts[su] = counts.get(su, 0) + 1
        if su not in first_seen:
            first_seen[su] = i

    items = [(s, counts[s], first_seen[s]) for s in counts]
    # sort by count desc, then by first_seen asc
    items.sort(key=lambda t: (-t[1], t[2]))
    return [(s, cnt) for s, cnt, _ in items]

def select_unique_panel(seq_counts: List[Tuple[str,int]],
                        max_unique: int,
                        min_weight_coverage: float) -> List[Tuple[str,int]]:
    total = sum(cnt for _, cnt in seq_counts) or 1
    out, acc = [], 0
    for s, cnt in seq_counts:
        if len(out) >= max_unique:
            break
        out.append((s, cnt))
        acc += cnt
        if acc / total >= min_weight_coverage:
            break
    return out

def run_muscle5_on_pairs(panel: List[Tuple[str,int]], prefix: str = "ins", threads: int = 1
    ) -> Tuple[MultipleSeqAlignment, Dict[str, int]]:
    """
    Returns:
      - Biopython MultipleSeqAlignment
      - weights_by_name: {sequence_name -> weight}
    """
    if len(panel) == 0:
        return MultipleSeqAlignment([]), {}
    if len(panel) == 1:
        s, w = panel[0]
        rec = SeqRecord(Seq(s), id=f"{prefix}_1", description="")
        return MultipleSeqAlignment([rec]), {f"{prefix}_1": w}

    pm_sequences, weights_by_name = [], {}
    for i, (s, w) in enumerate(panel, start=1):
        name = f"{prefix}_{i}"
        pm_sequences.append(pymuscle5.Sequence(name.encode(), s.encode()))
        weights_by_name[name] = w

    # The caller divides the requested thread budget across peak workers.
    aligner = pymuscle5.Aligner(threads=max(1, int(threads)))
    msa = aligner.align(pm_sequences)

    bio_records = [
        SeqRecord(Seq(seq.sequence.decode()), id=seq.name.decode(), description="")
        for seq in msa.sequences
    ]
    return MultipleSeqAlignment(bio_records), weights_by_name

def boundary_interval(allele, ref):
    """Every boundary p at which allele == ref[:p] + X + ref[p:] for some X.

    An insertion in a tandem repeat can be written at a run of equivalent
    positions, and the reads' aligner picks any of them. They form a single
    interval: p cannot pass the prefix the allele shares with ref, nor fall
    short of the shared suffix.
    """
    n = min(len(allele), len(ref))
    shared_prefix = 0
    while shared_prefix < n and allele[shared_prefix] == ref[shared_prefix]:
        shared_prefix += 1
    shared_suffix = 0
    while shared_suffix < n and allele[-1 - shared_suffix] == ref[-1 - shared_suffix]:
        shared_suffix += 1
    return max(0, len(ref) - shared_suffix), min(len(ref), shared_prefix)


def context_window(intervals, weights, ref_len, flank=CONTEXT_FLANK_BP,
                   share=CONTEXT_BOUNDARY_SHARE):
    """WT bases to trim from the two ends of each allele: (left, right).

    `intervals` gives each allele's boundary_interval and `weights` its read
    count. Every allele whose interval meets [left, ref_len - right] starts
    with ref[:left] and ends with ref[ref_len - right:], so those ends can be
    dropped. The span starts at the position most reads' intervals share,
    which copies of one ITD do even with sequencing errors in the insertion,
    and widens until it reaches `share` of the reads; a second ITD in the peak
    therefore widens it rather than being dropped. `flank` bases of anchor are
    kept on each side.
    """
    total = sum(weights)

    def reached(lo, hi):
        return sum(w for (a, b), w in zip(intervals, weights) if a <= hi and b >= lo)

    # The most-shared position is always the start of some interval.
    lo = hi = max(sorted({a for a, _ in intervals}), key=lambda p: reached(p, p))
    while reached(lo, hi) < share * total:
        # Widen towards the nearest allele not yet reached.
        _, a, b = min(
            [(a - hi, a, b) for a, b in intervals if a > hi]
            + [(lo - b, a, b) for a, b in intervals if b < lo]
        )
        lo, hi = min(lo, b), max(hi, a)
    left = max(0, lo - flank)
    right = max(0, ref_len - hi - flank)
    return left, right


@dataclass
class ColumnConsensus:
    """Weighted consensus of an MSA, with the per-column detail behind it."""
    consensus: str          # one called base per kept column; never a gap
    max_minor: float        # largest weighted share disagreeing with a call
    columns: np.ndarray     # MSA column behind each consensus base
    agreement: np.ndarray   # weighted share backing each called base
    matrix: np.ndarray      # (n_alleles, aln_len) aligned characters
    weights: np.ndarray     # read count of each aligned allele


def weighted_consensus_from_msa(
    msa, weights_by_name, *,
    base_threshold=0.7,
    min_col_coverage=0.8,
    ambiguous="N",
):
    """
    Compute a weighted consensus sequence from a Multiple Sequence Alignment (MSA).

    Parameters
    ----------
    msa : Bio.Align.MultipleSeqAlignment
        Multiple sequence alignment object.
    weights_by_name : dict
        {record.id: weight} dictionary (e.g. read counts per unique allele).
    base_threshold : float
        Minimum fraction of weighted votes for top base to be accepted.
    min_col_coverage : float
        Minimum weighted non-gap coverage fraction to include a column in consensus.
    ambiguous : str
        Symbol for ambiguous columns (default: 'N').

    Returns
    -------
    ColumnConsensus
        The weighted consensus (no gaps) and its per-column support. Its
        max_minor is what reveals an unbalanced mixture: an N appears only
        once the top base falls under `base_threshold`, so two haplotypes at
        80/20 give a clean consensus, no N, and max_minor near 0.2.
    """
    if len(msa) == 0:
        return ColumnConsensus("", 0.0, np.array([], dtype=int), np.array([]),
                               np.empty((0, 0), dtype="<U1"), np.array([]))

    # --- Prepare alignment matrix ---
    names = [rec.id for rec in msa]
    A = np.array([list(str(rec.seq)) for rec in msa])  # shape: (n_seq, aln_len)
    n, L = A.shape
    W = np.array([weights_by_name.get(name, 1.0) for name in names], dtype=float)
    total_w = W.sum() or 1.0

    # --- Compute per-column weighted coverage ---
    non_gap_w = np.array([(W[A[:, j] != "-"]).sum() for j in range(L)], dtype=float)
    coverage_frac = non_gap_w / total_w
    keep = coverage_frac >= min_col_coverage

    # Fallback: keep at least one column
    if not np.any(keep):
        keep[np.argmax(coverage_frac)] = True

    kept_positions = np.where(keep)[0]

    # --- Weighted consensus computation ---
    consensus_chars, agreement, columns = [], [], []
    for j in kept_positions:
        col = A[:, j]
        mask_col = (col != "-")
        if not mask_col.any():
            continue
        base_weights = {}
        for b, w in zip(col[mask_col], W[mask_col]):
            base_weights[b] = base_weights.get(b, 0.0) + w
        top_base, top_w = max(base_weights.items(), key=lambda kv: kv[1])
        frac = top_w / sum(base_weights.values())
        # An N only appears once the top base drops below base_threshold, so an
        # 80/20 mixture looks clean. Track the agreement itself.
        agreement.append(frac)
        columns.append(j)
        consensus_chars.append(top_base if frac >= base_threshold else ambiguous)

    return ColumnConsensus(
        consensus="".join(consensus_chars),
        max_minor=1.0 - min(agreement) if agreement else 0.0,
        columns=np.array(columns, dtype=int),
        agreement=np.array(agreement),
        matrix=A,
        weights=W,
    )

def plot_consensus_agreement(
    cc, *,
    alias,
    sample_name,
    out_dir,
    resolved,
    ins_start,
    ins_len,
    base_threshold,
    n_reads,
    panel_used,
    pad=10,
    wrap=100,
    max_alleles=8,
    dpi=150,
):
    """Per-base agreement across one ITD's consensus, for the HTML report.

    Shows the inserted bases with `pad` WT bases either side. Bars give the
    share of reads backing each called base, with the call threshold and the
    mixed-consensus level marked. Below them come the consensus and the most
    abundant aligned alleles, printing only where they differ from it, with
    '+' marking bases an allele has that the consensus lacks: a difference
    shared by many reads points to a second haplotype, while scattered single
    differences are sequencing noise. When the insertion cannot be placed, the
    whole aligned window is shown instead.
    """
    n_cols = len(cc.consensus)
    if n_cols == 0:
        return None
    if ins_start is not None:
        lo, hi = max(0, ins_start - pad), min(n_cols, ins_start + ins_len + pad)
    else:
        lo, hi = 0, n_cols
    blocks = [(b0, min(b0 + wrap, hi)) for b0 in range(lo, hi, wrap)]
    order = np.argsort(-cc.weights, kind="stable")[:max_alleles]
    n_rows = 1 + len(order)             # consensus, then alleles
    row_step = 0.2                      # text rows sit below the bars, in bar units
    first_row = -0.28
    bottom = first_row - n_rows * row_step
    mixed_line = 1.0 - MIXED_MINOR_FRACTION
    sample_label = sample_name or "Sample"
    # Columns the consensus dropped (too few reads with a base there) still hold
    # bases for some alleles. Each is drawn as '+' before the next kept column.
    extra_before = np.searchsorted(cc.columns, np.arange(cc.matrix.shape[1]))
    dropped = np.setdiff1d(np.arange(cc.matrix.shape[1]), cc.columns)
    # One block keeps to its own width; several share one scale across rows.
    span = wrap if len(blocks) > 1 else hi - lo
    width = min(14.0, max(9.0, 2.6 + 0.115 * span))

    if resolved:
        what = f"{ins_len} bp insertion" + (" (shaded)" if ins_start is not None else "")
    else:
        what = "consensus not resolved; whole aligned window shown"
    title = [
        f"{sample_label} · {alias}: {what} · {n_reads:,} reads, {panel_used} distinct "
        f"alleles aligned · largest disagreement {cc.max_minor:.1%}",
        "Bars: share of reads backing the called base.",
        "Rows: the most abundant aligned alleles, '·' where they match the consensus, "
        "'+' for bases it lacks.",
    ]
    if width >= 12:
        title = [title[0], f"{title[1]} {title[2]}"]
    title_in = 0.2 + 0.17 * len(title)
    fig_h = title_in + 0.25 + len(blocks) * (1.1 + 0.16 * n_rows)

    fig, axes = plt.subplots(len(blocks), 1, squeeze=False, dpi=dpi, figsize=(width, fig_h))
    for ax, (b0, b1) in zip(axes[:, 0], blocks):
        xs = np.arange(b0, b1)
        agree = cc.agreement[b0:b1]
        ax.bar(xs, agree, width=0.8, zorder=2, color=[
            "#c8323a" if 1.0 - a > MIXED_MINOR_FRACTION else "#9aa3ad" for a in agree
        ])
        if ins_start is not None:
            s0, s1 = max(b0, ins_start), min(b1, ins_start + ins_len)
            if s0 < s1:
                ax.axvspan(s0 - 0.5, s1 - 0.5, color="#f6dcb4", alpha=0.5, zorder=0)
        ax.axhline(base_threshold, color="#555555", ls="--", lw=0.8, zorder=1)
        ax.axhline(mixed_line, color="#c8323a", ls=":", lw=0.8, zorder=1)

        for x in xs:
            base = cc.consensus[x]
            ax.text(x, first_row, base, ha="center", va="center", fontsize=8,
                    family="monospace", fontweight="bold",
                    color=BASE_COLOURS.get(base, "#666666"))
        ax.text(b1 - 0.3, first_row, "  consensus", ha="left", va="center",
                fontsize=7.5, color="#333333")
        for r, idx in enumerate(order, start=1):
            y = first_row - r * row_step
            for x, ch in zip(xs, cc.matrix[idx, cc.columns[b0:b1]]):
                same = ch == cc.consensus[x]
                ax.text(x, y, "·" if same else ch, ha="center", va="center",
                        fontsize=8, family="monospace",
                        color="#b8b8b8" if same else "#c8323a",
                        fontweight="normal" if same else "bold")
            extra = {int(extra_before[j]) for j in dropped if cc.matrix[idx, j] != "-"}
            for x in sorted(extra & set(range(b0, b1 + 1))):
                ax.text(x - 0.5, y, "+", ha="center", va="center", fontsize=8,
                        family="monospace", color="#c8323a", fontweight="bold")
            reads = int(cc.weights[idx])
            ax.text(b1 - 0.3, y, f"  {reads:,} read{'' if reads == 1 else 's'}",
                    ha="left", va="center", fontsize=7.5, color="#333333")

        ax.set_xlim(b0 - 0.6, b0 + span - 0.4)
        ax.set_ylim(bottom, 1.05)
        ax.set_yticks([0, base_threshold, mixed_line, 1.0])
        ax.set_yticklabels(
            ["0", f"{base_threshold:g} call", f"{mixed_line:.2f} mixed", "1"], fontsize=7
        )
        ax.spines["left"].set_bounds(0, 1)
        for side in ("right", "bottom"):
            ax.spines[side].set_visible(False)
        # Ticks count bases into the insertion (1 = its first base); flank
        # positions are 0 and below.
        origin = ins_start if ins_start is not None else 0
        ticks = [x for x in xs if x - origin == 0 or (x - origin + 1 > 0 and (x - origin + 1) % 10 == 0)]
        if ins_start is not None and b0 <= ins_start + ins_len - 1 < b1:
            ticks.append(ins_start + ins_len - 1)
        ax.xaxis.tick_top()
        ax.set_xticks(sorted(set(ticks)))
        ax.set_xticklabels([str(t - origin + 1) for t in sorted(set(ticks))], fontsize=7)

    fig.suptitle("\n".join(title), fontsize=9, x=0.1 / width, ha="left")
    # Margins in inches, so labels keep their room however wide the figure is.
    fig.subplots_adjust(left=0.85 / width, right=1 - 1.05 / width,
                        top=1 - title_in / fig_h, bottom=0.1 / fig_h, hspace=0.35)
    out_path = os.path.join(out_dir, f"{sample_label}_{alias}_MSA_consensus.png")
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    logger.info("[plot_consensus_agreement] Saved: %s", os.path.abspath(out_path))
    return out_path

def _consensus_for_peak(task):
    """Panel selection + MSA + weighted consensus for one peak.

    Module level so it can be pickled to a process pool. Peaks are independent;
    DataFrame-derived values stay in the parent and rows are reassembled in the
    parent's order, so parallel output is identical to serial.
    """
    (alias, seqs, max_unique, min_weight_coverage, base_threshold,
     min_col_coverage, ambiguous, out_dir, sample_name, ref_seq,
     min_itd_size, max_itd_size, msa_threads) = task

    t0 = time.perf_counter()
    seq_counts = dedup_with_counts(seqs)
    panel = select_unique_panel(
        seq_counts,
        max_unique=max_unique,
        min_weight_coverage=min_weight_coverage,
    )
    intervals = [boundary_interval(s, ref_seq.upper()) for s, _ in panel]
    left, right = 0, 0
    if panel:
        left, right = context_window(intervals, [w for _, w in panel], len(ref_seq))
    # Only alleles that can place their insertion inside the window share the
    # trimmed WT ends.
    end = len(ref_seq) - right
    inside = [(s, w) for (s, w), (a, b) in zip(panel, intervals) if a <= end and b >= left]
    if len(inside) < len(panel):
        logger.info(
            "[%s] %d of %d aligned alleles (%d reads) place the insertion outside "
            "the %d bp window and are left out of the MSA.", alias,
            len(panel) - len(inside), len(panel),
            sum(w for _, w in panel) - sum(w for _, w in inside), end - left,
        )
    windowed = [(s[left:len(s) - right], w) for s, w in inside]
    msa, weights_by_name = run_muscle5_on_pairs(windowed, prefix=alias, threads=msa_threads)
    cc = weighted_consensus_from_msa(
        msa,
        weights_by_name,
        base_threshold=base_threshold,
        min_col_coverage=min_col_coverage,
        ambiguous=ambiguous,
    )
    # Every aligned allele carries the trimmed WT ends unchanged, so putting them
    # back gives a full allele in reference coordinates.
    consensus = ""
    if cc.consensus:
        consensus = ref_seq[:left] + cc.consensus + ref_seq[len(ref_seq) - right:]

    # Consensus is an allele in reference context. Recover the insertion and
    # its boundary from the SAME alignment; never combine a rotated payload
    # with an independently computed median read boundary.
    insertion_seq = ""
    insertion_pos = None
    plot_start = None
    if consensus and set(consensus.upper()) - set("ACGT"):
        logger.warning("[%s] Unresolved bases in context consensus; skipping candidate.", alias)
    elif consensus:
        aln = build_default_aligner().align(ref_seq, consensus)[0]
        candidates = [(p, q, n) for p, q, n in find_insertions(aln)
                      if n >= min_itd_size]
        if len(candidates) == 1:
            p, q, n = candidates[0]
            if max_itd_size is None or n <= max_itd_size:
                insertion_seq, insertion_pos = consensus[q:q + n], p
                if 0 <= q - left <= len(cc.consensus) - n:
                    plot_start = q - left
        if not insertion_seq:
            logger.warning("[%s] Context consensus has no single insertion within size bounds; skipping candidate.", alias)

    try:
        plot_consensus_agreement(
            cc,
            alias=alias,
            sample_name=sample_name,
            out_dir=out_dir,
            resolved=bool(insertion_seq),
            ins_start=plot_start,
            ins_len=len(insertion_seq),
            base_threshold=base_threshold,
            n_reads=len(seqs),
            panel_used=len(inside),
        )
    except Exception as e:
        logger.warning("[%s] Could not draw the consensus plot: %s", alias, e)

    return {
        "alias": alias,
        "consensus": insertion_seq,
        "consensus_ins_pos_ref": insertion_pos,
        "allele_consensus_len": len(consensus),
        "n_unique": len(seq_counts),
        "panel_used": len(panel),
        "window_bp": len(ref_seq) - left - right,
        "max_minor_fraction": cc.max_minor,
        "elapsed": time.perf_counter() - t0,
    }


def build_itd_consensus_sequences(
    all_itd_insertions,
    comps,
    *,
    ref_seq,
    min_itd_size=12,
    max_itd_size=None,
    sample_name=None,
    max_unique=200,
    min_weight_coverage=0.98,
    base_threshold=0.7,
    min_col_coverage=0.8,
    ambiguous="N",
    threads=1,
    out_dir,
):
    """
    Build reference-context alleles from each (insertion boundary, payload),
    compute their weighted MSA consensus, and recover a coherent insertion.

    ref_seq is required: isolated payloads can be cyclic rotations of the same
    allele and cannot safely be aligned independently of their boundaries. Only
    the WT span where CONTEXT_BOUNDARY_SHARE of the reads can place their
    insertion, plus CONTEXT_FLANK_BP either side, is sent to MUSCLE; the rest
    is identical in every allele.

    Parameters
    ----------
    all_itd_insertions : pd.DataFrame
        Must contain columns ['peak_alias', 'ins_pos_ref', 'ins_seq'].
    comps : pd.DataFrame
        Must contain ['peak_alias', 'putative_itd_size', 'sd_bp'].
    max_unique : int
        Maximum number of unique insertion sequences sent to MUSCLE.
    min_weight_coverage : float
        Fraction of total reads (by count) to include before truncating unique panel.
    base_threshold : float
        Weighted consensus threshold for keeping the dominant base per column.
    min_col_coverage : float
        Minimum weighted coverage per column to retain in consensus.
    ambiguous : str
        Symbol to use for ambiguous columns (e.g., 'N').
    threads : int
        Worker processes for the per-peak MSA. Peaks are independent, so this
        only changes wall time, never the consensus produced for any peak.
    out_dir : str
        Directory to save results.

    Returns
    -------
    df_cons : pd.DataFrame
        Per-ITD consensus summary.
    """
    # Collect the per-peak work first, keeping the parent's peak order so the
    # assembled table does not depend on which worker finishes first.
    ref_seq = str(ref_seq)
    tasks, meta = [], []
    for _, row in comps.iterrows():
        alias = row["peak_alias"]
        if alias.upper() == "WT":
            continue

        itd_subset = all_itd_insertions.loc[
            all_itd_insertions["peak_alias"] == alias
        ]
        if itd_subset.empty:
            logger.info(f"[build_itd_consensus_sequences] No insertions for {alias}")
            continue

        seqs = []
        for ins in itd_subset.itertuples(index=False):
            if pd.isna(ins.ins_seq) or not ins.ins_seq:
                continue
            pos = int(ins.ins_pos_ref)
            if not 0 <= pos <= len(ref_seq):
                raise ValueError(f"Invalid insertion boundary {pos} for {alias}")
            seqs.append(ref_seq[:pos] + str(ins.ins_seq) + ref_seq[pos:])
        tasks.append((
            alias, seqs, max_unique, min_weight_coverage, base_threshold,
            min_col_coverage, ambiguous, out_dir, sample_name, ref_seq,
            min_itd_size, max_itd_size,
        ))
        meta.append({
            "alias": alias,
            "n_total_reads": len(seqs),
            "raw_median_ins_pos_ref": int(itd_subset["ins_pos_ref"].median()),
            "expected_itd_bp": row["putative_itd_size"],
            "sd_bp": row["sd_bp"],
        })

    n_workers = max(1, min(int(threads), len(tasks)))
    threads_per_peak = max(1, int(threads) // n_workers)
    tasks = [task + (threads_per_peak,) for task in tasks]
    t_all = time.perf_counter()
    if n_workers > 1:
        logger.info(
            "[build_itd_consensus_sequences] Building %d peak consensuses across "
            "%d workers.", len(tasks), n_workers,
        )
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            # map keeps results in submission order, so the table order is fixed
            outputs = list(ex.map(_consensus_for_peak, tasks))
    else:
        outputs = [_consensus_for_peak(t) for t in tasks]

    results = []
    for m, o in zip(meta, outputs):
        if not o["consensus"]:
            continue
        logger.info(
            f"[build_itd_consensus_sequences] {m['alias']}: "
            f"total_reads={m['n_total_reads']}, unique={o['n_unique']}, "
            f"panel_used={o['panel_used']}, window={o['window_bp']} bp of WT"
        )
        results.append({
            "peak_alias": m["alias"],
            "n_total_reads": m["n_total_reads"],
            "n_unique": o["n_unique"],
            "consensus_len": len(o["consensus"]),
            "consensus_seq": o["consensus"],
            "max_minor_fraction": round(float(o["max_minor_fraction"]), 4),
            "consensus_ins_pos_ref": o["consensus_ins_pos_ref"],
            # Retained as a compatibility alias; this is now the consensus anchor.
            "median_ins_pos_ref": o["consensus_ins_pos_ref"],
            "raw_median_ins_pos_ref": m["raw_median_ins_pos_ref"],
            "allele_consensus_len": o["allele_consensus_len"],
            "expected_itd_bp": m["expected_itd_bp"],
            "sd_bp": m["sd_bp"],
        })
        logger.info(
            f"[build_itd_consensus_sequences] {m['alias']}: "
            f"consensus_len={len(o['consensus'])}, elapsed={o['elapsed']:.2f}s"
        )
    logger.info(
        "[build_itd_consensus_sequences] All peak consensuses built in %.2fs.",
        time.perf_counter() - t_all,
    )
    # --- Save summary ---
    df_cons = pd.DataFrame(results)
    consensus_name = f"{sample_name}_itd_consensus_seq.tsv"
    out_tsv = os.path.join(out_dir, consensus_name)
    df_cons.to_csv(out_tsv, sep="\t", index=False)
    logger.info("[build_itd_consensus_sequences] Saved: %s", os.path.abspath(out_tsv))

    return df_cons
