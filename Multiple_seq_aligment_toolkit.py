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
import textwrap
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

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

def run_muscle5_on_pairs(panel: List[Tuple[str,int]], prefix: str = "ins"
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

    aligner = pymuscle5.Aligner()
    msa = aligner.align(pm_sequences)

    bio_records = [
        SeqRecord(Seq(seq.sequence.decode()), id=seq.name.decode(), description="")
        for seq in msa.sequences
    ]
    return MultipleSeqAlignment(bio_records), weights_by_name

def weighted_consensus_from_msa(
    msa, weights_by_name, *,
    base_threshold=0.7,
    min_col_coverage=0.8,
    ambiguous="N",
    create_plot=False,
    out_dir=".",
    sample_name=None,
    peak_alias=None,
    dpi=150,
    seq_track_ypos=0.05,
    fontsize=9,
):
    """
    Compute a weighted consensus sequence from a Multiple Sequence Alignment (MSA)
    and optionally generate a per-peak coverage + consensus plot.

    Parameters
    ----------
    msa : Bio.Align.MultipleSeqAlignment
        Multiple sequence alignment object.
    weights_by_name : dict
        {record.id: weight} dictionary (e.g. read counts per unique insertion).
    base_threshold : float
        Minimum fraction of weighted votes for top base to be accepted.
    min_col_coverage : float
        Minimum weighted non-gap coverage fraction to include a column in consensus.
    ambiguous : str
        Symbol for ambiguous columns (default: 'N').
    create_plot : bool
        If True, saves per-peak MSA coverage + consensus plot.
    out_dir : str
        Directory to save plots.
    sample_name : str or None
        Sample name for plot title and file name.
    peak_alias : str or None
        ITD peak alias (e.g., 'ITD_1') for labeling plots.
    dpi : int
        Plot resolution.
    seq_track_ypos : float
        Vertical offset for consensus letter track.
    fontsize : int
        Font size for consensus bases in the plot.

    Returns
    -------
    (consensus, max_minor_fraction) : tuple[str, float]
        The weighted consensus (no gaps), paired with the largest per-column
        fraction of weighted support that disagreed with the called base. That
        second value is what reveals an unbalanced mixture: an N appears only
        once the top base falls under `base_threshold`, so two haplotypes at
        80/20 give a clean consensus, no N, and max_minor_fraction near 0.2.
    """

    if len(msa) == 0:
        return "", 0.0

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
    consensus_chars = []
    minor_fractions = []
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
        # How much of this column disagrees with the base about to be called.
        # An N only appears once the top base drops below base_threshold, so a
        # lopsided mixture of two haplotypes -- 80/20, say -- yields a clean
        # consensus of the majority and no N at all. Tracking the disagreement
        # itself is what keeps that case visible.
        minor_fractions.append(1.0 - frac)
        consensus_chars.append(top_base if frac >= base_threshold else ambiguous)

    consensus = "".join(consensus_chars).replace("-", "")
    max_minor = max(minor_fractions) if minor_fractions else 0.0

    # --- Optional per-peak plot ---
    if create_plot:
        sample_label = sample_name or "Sample"
        peak_label = peak_alias or "Peak"
        wrapped_consensus = "\n".join(textwrap.wrap(consensus, width=60))
        plot_title = f"Sample: {sample_label} - {peak_label}\nConsensus (length={len(consensus)}):\n{wrapped_consensus}"

        fig, ax = plt.subplots(figsize=(10, 3), dpi=dpi)
        x = np.arange(1, L + 1)
        ax.bar(x, coverage_frac, width=0.9, alpha=0.8, color="gray")
        ax.axhline(min_col_coverage, color="red", linestyle="--", lw=1)
        ax.set_xlim(0.5, L + 0.5)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("MSA column")
        ax.set_ylabel("Weighted coverage (non-gap / total)")
        ax.set_title(plot_title, fontsize=10)

        # Highlight kept columns
        for j in kept_positions:
            ax.axvspan(j + 0.5, j + 1.5, alpha=0.15, color="blue")

        # Annotate consensus bases
        if kept_positions.size > 0:
            ymax = ax.get_ylim()[1]
            y_text = ymax * seq_track_ypos
            for j, base in zip(kept_positions, consensus_chars):
                ax.text(j + 1, y_text, base, ha="center", va="center", fontsize=fontsize, color="navy")

        fig.tight_layout()
        fname = f"{sample_label}_{peak_label}_MSA_consensus.png"
        out_path = os.path.join(out_dir, fname)
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        logger.info("[weighted_consensus_from_msa] Saved: %s", os.path.abspath(out_path))

    return consensus, max_minor

def _consensus_for_peak(task):
    """Panel selection + MSA + weighted consensus for one peak.

    This is lifted verbatim out of build_itd_consensus_sequences so it can be
    pickled to a process pool. Peaks are independent -- each one reads only its
    own insertion sequences and writes only its own plot -- so evaluating them
    concurrently cannot change any individual result. Everything derived from
    the DataFrame (read counts, median insertion position) stays in the parent,
    and rows are reassembled in the parent's original peak order, so the output
    is identical to the serial version rather than merely equivalent.
    """
    (alias, seqs, max_unique, min_weight_coverage, base_threshold,
     min_col_coverage, ambiguous, out_dir, sample_name) = task

    t0 = time.perf_counter()
    seq_counts = dedup_with_counts(seqs)
    panel = select_unique_panel(
        seq_counts,
        max_unique=max_unique,
        min_weight_coverage=min_weight_coverage,
    )
    msa, weights_by_name = run_muscle5_on_pairs(panel, prefix=alias)
    consensus, max_minor = weighted_consensus_from_msa(
        msa,
        weights_by_name,
        base_threshold=base_threshold,
        min_col_coverage=min_col_coverage,
        ambiguous=ambiguous,
        create_plot=True,
        out_dir=out_dir,
        sample_name=sample_name,
        peak_alias=alias,
    )
    return {
        "alias": alias,
        "consensus": consensus,
        "n_unique": len(seq_counts),
        "panel_used": len(panel),
        "max_minor_fraction": max_minor,
        "elapsed": time.perf_counter() - t0,
    }


def build_itd_consensus_sequences(
    all_itd_insertions,
    comps,
    *,
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
    For each ITD peak in comps, deduplicate insertion sequences,
    perform weighted MSA, and extract weighted consensus.

    Parameters
    ----------
    all_itd_insertions : pd.DataFrame
        Must contain columns ['peak_alias', 'ins_seq'].
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

        seqs = itd_subset["ins_seq"].dropna().tolist()
        tasks.append((
            alias, seqs, max_unique, min_weight_coverage, base_threshold,
            min_col_coverage, ambiguous, out_dir, sample_name,
        ))
        meta.append({
            "alias": alias,
            "n_total_reads": len(seqs),
            # median taken here, on the DataFrame, exactly as before
            "median_ins_pos_ref": int(itd_subset["ins_pos_ref"].median()),
            "expected_itd_bp": row["putative_itd_size"],
            "sd_bp": row["sd_bp"],
        })

    n_workers = max(1, min(int(threads), len(tasks)))
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
        logger.info(
            f"[build_itd_consensus_sequences] {m['alias']}: "
            f"total_reads={m['n_total_reads']}, unique={o['n_unique']}, "
            f"panel_used={o['panel_used']}"
        )
        results.append({
            "peak_alias": m["alias"],
            "n_total_reads": m["n_total_reads"],
            "n_unique": o["n_unique"],
            "consensus_len": len(o["consensus"]),
            "consensus_seq": o["consensus"],
            "max_minor_fraction": round(float(o["max_minor_fraction"]), 4),
            "median_ins_pos_ref": m["median_ins_pos_ref"],
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
