import os
os.environ["MPLBACKEND"] = "Agg"       # disable any GUI backend
os.environ["DISPLAY"] = ""             # make sure Tk can't open a window
os.environ["TK_SILENCE_DEPRECATION"] = "1"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import argparse
import logging
import seaborn as sns
import shutil
import pymuscle5
import base64
import textwrap
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple
from venv import logger
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import Align, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.mixture import GaussianMixture
from scipy.stats import fisher_exact
from datetime import datetime

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
    consensus : str
        Weighted consensus sequence (no gaps).
    """

    if len(msa) == 0:
        return ""

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
        consensus_chars.append(top_base if frac >= base_threshold else ambiguous)

    consensus = "".join(consensus_chars).replace("-", "")

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
        print(f"[weighted_consensus_from_msa] Saved: {os.path.abspath(out_path)}")

    return consensus

def build_itd_consensus_sequences(
    all_itd_insertions,
    comps,
    *,
    max_unique=200,
    min_weight_coverage=0.98,
    base_threshold=0.7,
    min_col_coverage=0.8,
    ambiguous="N",
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
    out_dir : str
        Directory to save results.

    Returns
    -------
    df_cons : pd.DataFrame
        Per-ITD consensus summary.
    """
    results = []

    for _, row in comps.iterrows():
        alias = row["peak_alias"]
        if alias.upper() == "WT":
            continue
        mean_len = row["putative_itd_size"]

        itd_subset = all_itd_insertions.loc[
            all_itd_insertions["peak_alias"] == alias
        ]
        if itd_subset.empty:
            print(f"[build_itd_consensus_sequences] No insertions for {alias}")
            continue

        # --- Deduplicate and select coverage panel ---
        seqs = itd_subset["ins_seq"].dropna().tolist()
        seq_counts = dedup_with_counts(seqs)
        panel = select_unique_panel(
            seq_counts,
            max_unique=max_unique,
            min_weight_coverage=min_weight_coverage,
        )

        # --- Run MUSCLE alignment ---
        msa, weights_by_name = run_muscle5_on_pairs(panel, prefix=alias)

        # --- Weighted consensus ---
        consensus = weighted_consensus_from_msa(
            msa,
            weights_by_name,
            base_threshold=base_threshold,
            min_col_coverage=min_col_coverage,
            ambiguous=ambiguous,
            create_plot=True,
            out_dir=out_dir,
            sample_name=sample_name,
            peak_alias=alias
        )

        # Median insertion position for this ITD
        median_pos = int(itd_subset["ins_pos_ref"].median())

        results.append({
            "peak_alias": alias,
            "n_total_reads": len(seqs),
            "n_unique": len(seq_counts),
            "consensus_len": len(consensus),
            "consensus_seq": consensus,
            "median_ins_pos_ref": median_pos,
            "expected_itd_bp": mean_len,
            "sd_bp": row["sd_bp"],
        })

    # --- Save summary ---
    df_cons = pd.DataFrame(results)
    consensus_name = f"{sample_name}_itd_consensus_seq.tsv"
    out_tsv = os.path.join(out_dir, consensus_name)
    df_cons.to_csv(out_tsv, sep="\t", index=False)
    print(f"[build_itd_consensus_sequences] Saved: {os.path.abspath(out_tsv)}")

    return df_cons