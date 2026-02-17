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
from Pairwise_aligment_toolkit import chunk_iterable, process_chunk, build_default_aligner



def extract_itd_insertions_from_subset_parallel(
    reads_df: pd.DataFrame,
    read_ids_subset: list,
    ref_seq: str,
    peak_alias: str,
    comps: pd.DataFrame,
    threads: int = 4,
    itd_sd_factor: float = 1.0
) -> pd.DataFrame:
    """
    Strand-aware alignment for ITD subset.
    Uses read sequences from reads_df instead of passing them directly.

    Parameters
    ----------
    reads_df : pd.DataFrame
        Must contain 'read_id' and 'read_seq' columns.
    read_ids_subset : list
        Read IDs belonging to this ITD peak (subset).
    ref_seq : str
        Reference sequence for alignment.
    peak_alias : str
        ITD alias, e.g. 'ITD_1'.
    comps : pd.DataFrame
        GMM components with mean_bp, sd_bp, putative_itd_size, etc.
    threads : int
        Number of threads (chunks = threads).
    itd_sd_factor : float
        Acceptable deviation from ITD size ± (factor × SD).
    """

    # --- Skip WT ---
    if peak_alias.upper() == "WT":
        print(f"[INFO] Skipping WT subset alignment ({peak_alias})")
        return pd.DataFrame()

    # --- Retrieve expected ITD range ---
    row = comps.loc[comps["peak_alias"] == peak_alias].iloc[0]
    itd_mean = row["putative_itd_size"]
    itd_sd = row["sd_bp"]

    itd_min = itd_mean - itd_sd * itd_sd_factor
    itd_max = itd_mean + itd_sd * itd_sd_factor
    logger.debug(f"Processing {peak_alias}: mean ITD size = {itd_mean:.1f} bp, SD = {itd_sd:.1f} bp. Range: {itd_min:.1f} - {itd_max:.1f} bp")
    # --- Get sequences for this subset ---
    reads_subset_df = reads_df.loc[reads_df["read_id"].isin(read_ids_subset), ["read_id", "read_seq"]]
    reads_list = list(reads_subset_df.itertuples(index=False, name=None))  # [(id, seq), ...]

    # --- Chunk and parallelize ---
    chunks = list(chunk_iterable(reads_list, threads))
    results = []
    # with ThreadPoolExecutor(max_workers=threads) as ex:
    #     futures = {ex.submit(process_chunk, chunk): i for i, chunk in enumerate(chunks)}
    #     for fut in as_completed(futures):
    #         results.extend(fut.result())
    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = {ex.submit(process_chunk, chunk, itd_min, itd_max, ref_seq, peak_alias ): i for i, chunk in enumerate(chunks)}
        for fut in as_completed(futures):
            results.extend(fut.result())

    if not results:
        return pd.DataFrame(columns=[
            "peak_alias", "read_id", "strand", "aln_score", "pct_identity",
            "ins_pos_ref", "ins_len", "ins_seq", "fwd_score", "rev_score"
        ])

    df = pd.DataFrame(results)
    df.sort_values(["ins_pos_ref", "ins_len"], inplace=True, ignore_index=True)
    return df

def plot_itd_size_distribution(all_itd_insertions, out_dir, sample_name=None, bins=50):
    """
    Plot histogram and KDE of putative ITD sizes across all reads/peaks.

    Parameters
    ----------
    all_itd_insertions : pd.DataFrame
        Must contain 'ins_len' and optionally 'peak_alias'.
    out_dir : str
        Directory to save the plot.
    sample_name : str or None
        Optional sample label for title and filename.
    bins : int
        Number of histogram bins.
    """
    if all_itd_insertions.empty:
        print("[plot_itd_size_distribution] No insertions found. Skipping plot.")
        return

    plt.figure(figsize=(8, 5), dpi=150)

    # --- Main histogram + KDE ---
    sns.histplot(
        data=all_itd_insertions,
        x="ins_len",
        bins=bins,
        hue="peak_alias",
        multiple="stack",
        kde=True,
        alpha=0.6,
        edgecolor=None,
    )

    plt.xlabel("Insertion length (bp)")
    plt.ylabel("Read count")
    title = f"ITD insertion length distribution"
    if sample_name:
        title += f" — Sample:{sample_name}"
    plt.title(title)

    # Add mean ± SD lines
    mean_len = all_itd_insertions["ins_len"].mean()
    sd_len = all_itd_insertions["ins_len"].std()
    plt.axvline(mean_len, color="black", linestyle="--", lw=1.2, label=f"Mean = {mean_len:.1f}")
    plt.axvspan(mean_len - sd_len, mean_len + sd_len, alpha=0.15, color="gray", label=f"±SD = {sd_len:.1f}")

    plt.legend()
    plt.tight_layout()

    fname = os.path.join(out_dir, f"{sample_name or 'sample'}_itd_size_distribution.png")
    plt.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"[plot_itd_size_distribution] Saved: {os.path.abspath(fname)}")

def build_itd_reference_per_peak(
    df_cons,
    ref_seq: str,
    genome_build: str,
    sample_name: str,
    amplicon_coords: dict,
    exon_boundaries: list,
    out_dir: str,
    exon_labels: list,
    flank_bp: int = 100,

):
    """
    Build ITD reference sequences per peak, map insertion to genomic coordinates,
    and produce per-ITD plots showing its genomic position within FLT3.

    Parameters
    ----------
    df_cons : pd.DataFrame
        Consensus dataframe with at least:
        ['peak_alias', 'consensus', 'median_insertion_pos']
    ref_seq : str
        Reference amplicon sequence (covering the FLT3 JM domain region).
    genome_build : str
        'hg19' or 'hg38'.
    amplicon_coords : dict
        {build: {'chr': 'chr13', 'start': int, 'end': int, 'strand': +1}}
        Coordinates are UCSC-style 0-based (BED-like).
    exon_boundaries : list[(start, end)]
        Exon coordinates (UCSC 0-based) for plotting context.
    out_dir : str
        Directory for output plots.
    exon_labels : list[str], optional
        Names/numbers of the exons for labeling.
    flank_bp : int
        Bases around the insertion site to show in plots.

    Returns
    -------
    itd_ref_dict : dict
        {peak_alias: {'alias', 'itd_seq', 'ref_seq_with_itd',
                      'local_insertion_pos', 'genomic_insertion_pos', 'chr'}}
    """

    amp = amplicon_coords[genome_build]
    itd_ref_dict = {}

    for _, row in df_cons.iterrows():
        alias = row["peak_alias"]
        itd_seq = row["consensus_seq"]
        ins_pos_local = int(row["median_ins_pos_ref"])

        # --- Genomic coordinate conversion (UCSC 0-based → human-readable 1-based)
        genomic_ins_pos = amp["start"] + ins_pos_local  # still 0-based
        chr_name = amp["chr"]

        # --- Build ITD-inserted reference
        itd_ref_seq = ref_seq[:ins_pos_local] + itd_seq + ref_seq[ins_pos_local:]
        itd_len = len(itd_seq)

        itd_ref_dict[alias] = {
            "alias": alias,
            "itd_seq": itd_seq,
            "ref_seq_with_itd": itd_ref_seq,
            "local_insertion_pos": ins_pos_local,
            "genomic_insertion_pos": genomic_ins_pos,
            "chr": chr_name,
            "itd_length": itd_len,
        }

        # --- Visualization per ITD ---
        plot_itd_vs_ref_with_genome(
            itd_ref_seq=itd_ref_seq,
            ref_seq=ref_seq,
            sample_name=sample_name,
            peak_alias=alias,
            out_dir=out_dir,
            genomic_info=amp,
            exon_boundaries=exon_boundaries,
            genome_build=genome_build,
            insertion_genomic_pos=genomic_ins_pos,
            exon_labels=exon_labels,
            flank_bp=flank_bp,
        )

    return itd_ref_dict

def plot_itd_vs_ref_with_genome(
    itd_ref_seq: str,
    ref_seq: str,
    sample_name: str,
    peak_alias: str,
    out_dir: str,
    genomic_info: dict,
    exon_boundaries: list,
    genome_build: str,
    insertion_genomic_pos: int,
    exon_labels: list,
    flank_bp: int = 100,
    dpi: int = 150,
):
    """
    Visualize the ITD insertion within the FLT3 genomic context.

    Parameters
    ----------
    itd_ref_seq : str
        Reference sequence with the ITD inserted.
    ref_seq : str
        Original amplicon reference sequence.
    sample_name : str
        Sample name for labeling.
    peak_alias : str
        ITD alias (e.g., "ITD_1").
    out_dir : str
        Directory to save the figure.
    genomic_info : dict
        {"chr": "chr13", "start": int, "end": int, "strand": +1}
        Coordinates of the amplicon (0-based UCSC convention).
    exon_boundaries : list[(start, end)]
        List of exon intervals (0-based).
    genome_build : str
        "hg19" or "hg38".
    insertion_genomic_pos : int
        0-based genomic coordinate of the ITD insertion site.
    exon_labels : list[str], optional
        Labels for exons (e.g., ["Ex1", "Ex2", ...]).
    flank_bp : int
        Number of base pairs upstream/downstream to display.
    dpi : int
        Plot resolution.
    """
    chr_name = genomic_info["chr"]
    amp_start = genomic_info["start"]
    amp_end = genomic_info["end"]

    itd_len = len(itd_ref_seq) - len(ref_seq)

    # Genomic position for display (1-based UCSC)
    ins_pos_display = insertion_genomic_pos + 1

    fig, ax = plt.subplots(figsize=(15, 2.5), dpi=dpi)

    # Determine plot range
    left = insertion_genomic_pos - (flank_bp)
    right = insertion_genomic_pos + flank_bp
    # ----- Plot only exons 10 and 11 -----
    for i, (start, end) in enumerate(exon_boundaries):
        if i not in [9, 10]:  # Skip all exons except 10 and 11
            continue
        rect = mpatches.Rectangle(
            (start, 0.4),
            end - start,
            0.2,
            color="green",
            alpha=0.4,
            lw=0
        )
        ax.add_patch(rect)

        # Draw exon labels, but clip them if they are too close to the right edge
        if exon_labels and i < len(exon_labels):
            xpos = min((start + end) / 2, right - 100)  # keep label visible
            ax.text(xpos, 0.65, exon_labels[i],
                    ha="center", va="bottom", fontsize=7, color="darkgreen",
                    clip_on=True)

    # ----- Highlight amplicon -----
    ax.add_patch(mpatches.Rectangle(
        (amp_start, 0.3),
        amp_end - amp_start,
        0.4,
        color="orange",
        alpha=0.3,
        lw=0
    ))

    # ----- Mark insertion site -----
    ax.axvline(insertion_genomic_pos, color="red", lw=2, linestyle="--", label="ITD insertion")

    # Annotate ITD
    ax.text(insertion_genomic_pos, 0.05,
            f"ITD ({itd_len} bp)",
            ha="center", va="bottom", color="red", fontsize=8)

    # ----- Axis setup -----

    ax.set_xlim(left, right)

    # Use readable genomic coordinates with commas
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))

    ax.set_ylim(0, 1)
    ax.set_xlabel(f"Genomic position on {chr_name} ({genome_build})")
    ax.set_yticks([])

    # Improved title placement
    ax.set_title(
        f"{sample_name} – {peak_alias}\nInsertion at {chr_name}:{ins_pos_display:,} (ITD {itd_len} bp)",
        pad=35, fontsize=11
    )

    # ----- Legend -----
    handles = [
        mpatches.Patch(color="green", alpha=0.4, label="Exons"),
        mpatches.Patch(color="orange", alpha=0.3, label="Amplicon"),
        mpatches.Patch(color="red", alpha=0.4, label="ITD insertion"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=False)

    # ----- Save -----
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{sample_name}_{peak_alias}_itd_ref_{genome_build}.png")

    # Balanced layout adjustments
    fig.subplots_adjust(top=0.82, bottom=0.22)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)

    print(f"[plot_itd_vs_ref_with_genome] Saved plot: {os.path.abspath(out_path)}")

    aligner = build_default_aligner()
    alignment = aligner.align(ref_seq, itd_ref_seq)[0]
    aln_score = alignment.score

    aln_out_path = os.path.join(
        out_dir,
        f"{sample_name}_pairwisealignment_{peak_alias}.txt"
    )

    # Write a readable alignment block and metadata
    with open(aln_out_path, "w") as f:
        f.write(f"# Sample: {sample_name}\n")
        f.write(f"# Peak alias: {peak_alias}\n")
        f.write(f"# Genome build: {genome_build}\n")
        f.write(f"# Insertion position: {chr_name}:{ins_pos_display:,}\n")
        f.write(f"# ITD length: {itd_len} bp\n")
        f.write(f"# Alignment score: {aln_score}\n")
        f.write("#" + "-" * 70 + "\n\n")
        f.write(str(alignment))
        f.write("\n")

def make_validation_refs(ref_seq, itd_ref_dict, sample_name, out_dir):
    """
    Build a multi-reference FASTA (WT + all ITDs) for competitive alignment validation.

    Parameters
    ----------
    ref_seq : str
        Wild-type reference amplicon sequence.
    itd_ref_dict : dict
        Dictionary of ITD reference info as produced by build_itd_reference_per_peak.
    sample_name : str
        Sample identifier (used for file naming).
    out_dir : str
        Output directory for saving the FASTA.

    Returns
    -------
    ref_dict : dict
        {"WT": <ref_seq>, "ITD_1": <itd_seq_with_insertion>, ...}
    fasta_path : str
        Absolute path to saved multi-reference FASTA file.
    """
    os.makedirs(out_dir, exist_ok=True)
    fasta_path = os.path.join(out_dir, f"{sample_name}_validation_refs.fasta")

    # Prepare records
    records = []
    ref_dict = {}

    # 1. WT entry
    wt_record = SeqRecord(Seq(ref_seq), id="WT", description="Wild-type amplicon reference")
    records.append(wt_record)
    ref_dict["WT"] = {
        "ref_seq_with_itd": ref_seq,
        "chr": itd_ref_dict[next(iter(itd_ref_dict))]["chr"] if itd_ref_dict else None,
        "genomic_insertion_pos": None,
        "itd_length": 0,
    }
    # --- 2. ITD references ---
    for alias, entry in itd_ref_dict.items():
        itd_seq = entry["ref_seq_with_itd"]
        rec = SeqRecord(
            Seq(itd_seq),
            id=alias,
            description=f"ITD {alias}, len={len(itd_seq)} bp, ins@{entry['genomic_insertion_pos']}",
        )
        records.append(rec)

        ref_dict[alias] = {
            "ref_seq_with_itd": itd_seq,
            "chr": entry.get("chr"),
            "genomic_insertion_pos": entry.get("genomic_insertion_pos"),
            "itd_length": entry.get("itd_length", len(itd_seq)),
        }

    # Save multi-fasta
    SeqIO.write(records, fasta_path, "fasta")

    logger.info(f"Saved multi-reference FASTA: {os.path.abspath(fasta_path)}")
    logger.info(f"References included: {list(ref_dict.keys())}")

    return ref_dict, fasta_path

def prepare_validation_reads(reads_df):
    """
    Prepare high-confidence reads for ITD validation alignment.

    Parameters
    ----------
    reads_df : pd.DataFrame
        Output from fit_gmm_itds(), must include:
        'read_id', 'read_seq', 'read_len', 'gmm_peak_alias', 'is_ambiguous'.

    Returns
    -------
    filtered_df : pd.DataFrame
        Cleaned DataFrame containing only reads to use for validation.
    """

    # --- Sanity checks ---
    required_cols = {"read_id", "read_seq", "gmm_peak_alias", "is_ambiguous"}
    missing = required_cols - set(reads_df.columns)
    if missing:
        raise ValueError(f"reads_df missing required columns: {missing}")

    # --- Filter reads ---
    filtered_df = reads_df.loc[
        (~reads_df["is_ambiguous"]) & (reads_df["gmm_peak_alias"].notnull())
    ].copy()

    # Drop empty sequences (if any)
    filtered_df = filtered_df[filtered_df["read_seq"].str.len() > 0]

    # Log summary
    n_total = len(reads_df)
    n_used = len(filtered_df)
    logger.info(f"[prepare_validation_reads] Using {n_used}/{n_total} reads ({100*n_used/n_total:.1f}%) for validation.")

    return filtered_df

def calculate_allele_frequencies_and_strand_bias(validation_results_df, min_allele_frequency: float):

    df_valid = validation_results_df.query(
        "valid_support == True and filter_reason == 'pass'"
    ).copy()

    # Count ITD-supporting and WT-supporting reads
    itd_counts = (
        df_valid[df_valid["support_call"] == "ITD-supporting"]
        .groupby("ref_alias")
        .agg(n_itd_reads=("read_id", "nunique"))
        .reset_index() 
    )

    itd_counts["n_total_reads"] = len(df_valid)
    itd_counts["allele_frequency"] = itd_counts["n_itd_reads"] / itd_counts["n_total_reads"]

    # Calculate strand bias, make sure ref_alias is stored correctly

    strand_counts = (
        df_valid.groupby(["ref_alias", "strand"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    # Ensure both "+" and "-" columns exist
    if "+" not in strand_counts.columns:
        strand_counts["+"] = 0
    if "-" not in strand_counts.columns:
        strand_counts["-"] = 0

    strand_stats = []
    for _, row in strand_counts.iterrows():
        alias = str(row["ref_alias"])  # ✅ correct alias name
        plus = int(row["+"])
        minus = int(row["-"])
        total = plus + minus
        if total == 0:
            logger.warning(f"No reads found for {alias}, skipping strand bias calculation.")
            continue

        ratio = plus / (total + 1e-8)
        # Fisher’s exact test
        wt_plus = df_valid.query("ref_alias == 'WT' and strand == '+'").shape[0]
        wt_minus = df_valid.query("ref_alias == 'WT' and strand == '-'").shape[0]
        oddsratio, pval = fisher_exact([[plus, minus], [wt_plus, wt_minus]])

        strand_stats.append({
            "ref_alias": alias,
            "plus_reads": plus,
            "minus_reads": minus,
            "strand_ratio": ratio,
            "odds_ratio": oddsratio,
            "fisher_p": pval
        })

    strand_stats = pd.DataFrame(strand_stats)

    logger.info("Applying allele frequency filter to detected ITDs...")
    logger.info(f"ITD counts before applying min allele frequency of {min_allele_frequency}:")
    logger.info(itd_counts)

    itd_counts = itd_counts[itd_counts["allele_frequency"] >= min_allele_frequency].copy()

    strand_stats["ref_alias"] = strand_stats["ref_alias"].astype(str)
    itd_counts = itd_counts.reset_index()
    itd_counts["ref_alias"] = itd_counts["ref_alias"].astype(str)

    logger.info(f"Detected ITDs after applying min allele frequency of {min_allele_frequency}:")
    logger.info(itd_counts)

    summary_df = itd_counts.merge(strand_stats, on="ref_alias", how="left")
    summary_df = summary_df.sort_values("allele_frequency", ascending=False)

    logger.debug("printing summary_df")
    logger.debug(type(summary_df))
    logger.debug(summary_df)

    return summary_df
