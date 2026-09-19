import os
os.environ["MPLBACKEND"] = "Agg"       # disable any GUI backend
os.environ["DISPLAY"] = ""             # make sure Tk can't open a window
os.environ["TK_SILENCE_DEPRECATION"] = "1"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import logging
import time
import seaborn as sns
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import fisher_exact
logger = logging.getLogger(__name__)

def build_default_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -5
    aligner.target_open_gap_score = -15
    aligner.target_extend_gap_score = -0.05
    aligner.query_open_gap_score = -200
    aligner.query_extend_gap_score = -1.5
    aligner.query_end_gap_score = 0
    aligner.target_end_gap_score = 0
    aligner.mode = "global"
    return aligner

def insertions_from_aligned_blocks(aln_blocks):
    """Detect insertions (in query) as large jumps between aligned blocks.

    Blocks returned by Bio.Align are gapless by construction: inside a single
    block the target and query spans are always the same length. All indel
    evidence therefore lives in the jumps *between* consecutive blocks, never
    inside one -- so any check that subtracts the two spans within a block is
    identically zero and silently passes everything.

    Parameters
    ----------
    aln_blocks : array-like or None
        An ``Alignment.aligned`` value, shape (2, n_blocks, 2), or an empty
        sequence when no alignment was produced.

    Returns
    -------
    list[tuple[int, int, int]]
        (target_end, query_end, insertion_length) for every inter-block gap
        where the query advanced further than the target.
    """
    if aln_blocks is None or len(aln_blocks) != 2:
        return []

    t_blocks, q_blocks = aln_blocks[0], aln_blocks[1]
    insertions = []

    for i in range(len(t_blocks) - 1):
        t_end = int(t_blocks[i][1])
        t_next = int(t_blocks[i + 1][0])
        q_end = int(q_blocks[i][1])
        q_next = int(q_blocks[i + 1][0])

        # If query advanced more than target between blocks -> insertion
        ins_len = (q_next - q_end) - (t_next - t_end)
        if ins_len > 0:
            insertions.append((t_end, q_end, ins_len))
    return insertions

def find_insertions_between_blocks(aln):
    """Detect insertions (in query) as large jumps between aligned blocks."""
    return insertions_from_aligned_blocks(aln.aligned)

def chunk_iterable(data, n):
    """Split list into n roughly equal chunks."""
    k = max(1, len(data) // n)
    for i in range(0, len(data), k):
        yield data[i:i + k]

def percent_identity(aln):
    ref = aln.target
    qry = aln.query
    ident = 0
    total = 0
    for (rs, re), (qs, qe) in zip(*aln.aligned):
        for r_i, q_i in zip(range(rs, re), range(qs, qe)):
            total += 1
            if ref[r_i] == qry[q_i]:
                ident += 1
    return ident / total if total else 0.0

def find_insertions(aln):
    # Insertions only ever appear between aligned blocks (see
    # insertions_from_aligned_blocks); there is no within-block case to fall
    # back to.
    return find_insertions_between_blocks(aln)

def _seq_in_reference_orientation(read_seq: str, read_strand: str) -> str:
    """
    Use stored sequence as-is for alignment.
    With cutadapt --rc, sequences marked with 'rc' are already reverse-complemented
    in output when that orientation gives the better adapter match.
    """
    return read_seq

def process_chunk(chunk, itd_min, itd_max, ref_seq, peak_alias):
    aligner = build_default_aligner()
    out_rows = []
    for read_id, read_seq, read_strand in chunk:
        try:
            aligned_seq = _seq_in_reference_orientation(read_seq, read_strand)
            aln = aligner.align(ref_seq, aligned_seq)[0]
            pid = percent_identity(aln)
            ins_regions = find_insertions(aln)

            for ts, qs, ins_len in ins_regions:
                if itd_min <= ins_len <= itd_max:
                    ins_seq = aligned_seq[qs: qs + ins_len]
                    out_rows.append({
                        "peak_alias": peak_alias,
                        "read_id": read_id,
                        "strand": read_strand,
                        "aln_score": aln.score,
                        "pct_identity": pid,
                        "ins_pos_ref": ts,
                        "ins_len": ins_len,
                        "ins_seq": ins_seq,
                        "fwd_score": aln.score,
                        "rev_score": np.nan
                    })
        except Exception as e:
            logger.warning(f"[WARN] Alignment failed for {read_id}: {e}")
    return out_rows

def compute_adjusted_score(aln, alpha=0.5):
    """Compute hybrid PID + gap-penalized score from a Biopython alignment."""
    if aln is None:
        return 0.0

    try:
        # --- Percent identity (0-1) ---
        pid = percent_identity(aln)

        # --- Find insertions (query gaps between blocks) ---
        insertions = find_insertions_between_blocks(aln)
        total_ins = sum(ins_len for _, _, ins_len in insertions)

        # --- Compute aligned length from target blocks ---
        t_blocks, _ = aln.aligned
        aln_len = sum(e - s for s, e in t_blocks) or 1
        ref_len = len(aln.target)
        coverage = aln_len / ref_len

        # --- Gap ratio and adjusted score ---
        gap_ratio = total_ins / aln_len
        adjusted = (pid - alpha * gap_ratio) * coverage
        if total_ins > 0:
            logger.debug(f"Insertion(s) detected: total={total_ins}bp, pid={pid:.3f}, adjusted={adjusted:.3f}")

        return max(adjusted, 0.0)

    except Exception as e:
        logger.error(f"[compute_adjusted_score] Failed on alignment: {e}")
        return 0.0

def softmax(x, beta):
    """Stable softmax with temperature scaling."""
    x = np.array(x, dtype=float)
    x = x - np.max(x)  # numerical stability
    exp_x = np.exp(beta * x)
    return exp_x / (np.sum(exp_x) + 1e-8)

def classify_read_support(row, metric_used, z_thresh=1.0, delta_thresh=0.05):
    """
    Classify reads as WT-supporting, ITD-supporting, or Ambiguous.
    Assumes ref_alias is always either "WT" or "ITD_N" (e.g., ITD_1, ITD_2).
    """
    alias = row["ref_alias"]
    val = row["metric_value"]

    # --- Z-score metric ---
    if metric_used == "z_score":
        if abs(val) >= z_thresh:
            return "ITD-supporting" if "ITD" in alias else "WT-supporting"
        else:
            return "Ambiguous"

    # --- Delta-style metrics (delta_score or prob_delta) ---
    elif metric_used in ("delta_score", "prob_delta"):
        delta = row.get("delta", 0)
        if delta >= delta_thresh:
            return "ITD-supporting" if alias.startswith("ITD") else "WT-supporting"
        elif delta <= -delta_thresh:
            return "WT-supporting" if alias.startswith("ITD") else "ITD-supporting"
        else:
            return "Ambiguous"

    # --- Default case ---
    return "Ambiguous"

def validate_itd_supporting_reads(
    df_best,
    df_cons,
    *,
    gap_window,
    max_gap_bp,
    min_pid,
    max_wt_gap_bp=10,
):
    """
    Validate both WT and ITD reads by inspecting precomputed alignments.
    """
    validated = []
    reasons = []

    for _, row in df_best.iterrows():
        alias = row["ref_alias"]
        aln_blocks = row.get("aligned_blocks")
        pid = row.get("pct_identity", 0)

        # --- Handle WT reads ---
        if alias == "WT":
            if aln_blocks is None or len(aln_blocks) != 2:
                validated.append(False)
                reasons.append("missing_alignment_blocks")
                continue

            # Total inserted bases relative to WT, measured between blocks.
            total_gap_len = sum(
                ins_len
                for _t_end, _q_end, ins_len in insertions_from_aligned_blocks(aln_blocks)
            )

            if pid < min_pid:
                validated.append(False)
                reasons.append(f"low_pid({pid:.2f})")
            elif total_gap_len > max_wt_gap_bp:
                validated.append(False)
                reasons.append(f"wt_with_gaps({total_gap_len} bp)")
            else:
                validated.append(True)
                reasons.append("pass")
            continue

        # --- Handle ITD reads ---
        elif alias.startswith("ITD"):
            ins_row = df_cons.loc[df_cons["peak_alias"] == alias]
            if ins_row.empty:
                validated.append(False)
                reasons.append("missing_insertion_position")
                continue

            ins_pos = int(ins_row["median_ins_pos_ref"].iloc[0])
            itd_len = int(ins_row["consensus_len"].iloc[0]) if "consensus_len" in ins_row else 0

            if aln_blocks is None or len(aln_blocks) != 2:
                validated.append(False)
                reasons.append("missing_alignment_blocks")
                continue

            # Residual insertions in the read relative to the ITD reference.
            # The read is aligned to wt[:ins_pos] + itd + wt[ins_pos:], so in
            # ITD-reference coordinates the duplicated segment spans
            # [ins_pos, ins_pos + itd_len); anything within gap_window of that
            # span counts as sitting on a breakpoint.
            lo = ins_pos - gap_window
            hi = ins_pos + itd_len + gap_window
            gaps_near = sum(
                ins_len
                for t_end, _q_end, ins_len in insertions_from_aligned_blocks(aln_blocks)
                if lo <= t_end <= hi
            )

            if pid < min_pid:
                validated.append(False)
                reasons.append(f"low_pid({pid:.2f})")
            elif gaps_near == 0:
                validated.append(True)
                reasons.append("pass")
            elif gaps_near > max_gap_bp:
                validated.append(False)
                reasons.append(f"too_many_gaps({gaps_near})")
            else:
                validated.append(True)
                reasons.append("pass")
            continue

        # --- Handle unexpected aliases ---
        else:
            validated.append(False)
            reasons.append(f"unknown_alias:{alias}")
            continue

    # ---- Finalize ----
    df_best["valid_support"] = validated
    df_best["filter_reason"] = reasons

    total = len(df_best)
    passed = df_best["valid_support"].sum()
    failed = total - passed

    logger.debug(f"[validate_itd_supporting_reads] {passed}/{total} reads passed validation.")
    logger.debug(f"Failed reads: {failed}")
    logger.debug("Reason breakdown:")
    logger.debug(df_best["filter_reason"].value_counts().to_string())

    itd_summary = (
        df_best[df_best["ref_alias"].str.startswith("ITD")]
        .groupby(["ref_alias", "valid_support"])
        .size()
        .unstack(fill_value=0)
    )
    if not itd_summary.empty:
        logger.debug("\nPer-ITD validation summary:")
        logger.debug(itd_summary.to_string())

    return df_best

def count_gaps_near(aln, insertion_pos, window=15):
    """Count total gap bases in the reference +/-window bp around the insertion site."""
    ref = aln.target
    qry = aln.query
    total_gap_bp = 0
    ref_pos = -1
    for r, q in zip(ref, qry):
        if r != "-":
            ref_pos += 1
        if (insertion_pos - window) <= ref_pos <= (insertion_pos + window):
            if q == "-":
                total_gap_bp += 1
    return total_gap_bp



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
        Must contain 'read_id', 'read_seq', and 'strand' columns.
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
    reads_subset_df = reads_df.loc[reads_df["read_id"].isin(read_ids_subset), ["read_id", "read_seq", "strand"]]
    reads_list = list(reads_subset_df.itertuples(index=False, name=None))  # [(id, seq), ...]

    # --- Chunk and parallelize ---
    if threads < 1:
        threads = 1
    chunks = list(chunk_iterable(reads_list, threads))
    n_reads = len(reads_list)
    n_batches = len(chunks)
    est_comparisons = n_reads
    logger.info(
        f"[extract_itd_insertions_from_subset_parallel] Starting peak={peak_alias}: "
        f"reads={n_reads}, batches={n_batches}, workers={threads}, "
        f"estimated_comparisons={est_comparisons}"
    )

    t0 = time.perf_counter()
    results = []
    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = {ex.submit(process_chunk, chunk, itd_min, itd_max, ref_seq, peak_alias ): i for i, chunk in enumerate(chunks)}
        for fut in as_completed(futures):
            results.extend(fut.result())
    elapsed_sec = time.perf_counter() - t0
    reads_per_sec = (n_reads / elapsed_sec) if elapsed_sec > 0 else 0.0
    comps_per_sec = (est_comparisons / elapsed_sec) if elapsed_sec > 0 else 0.0
    logger.info(
        f"[extract_itd_insertions_from_subset_parallel] Finished peak={peak_alias} in "
        f"{elapsed_sec:.2f}s ({reads_per_sec:.1f} reads/s, {comps_per_sec:.1f} est comparisons/s)"
    )

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
    # Add per-peak mean +/- SD overlays (not pooled across peaks)
    peak_aliases = [p for p in all_itd_insertions["peak_alias"].dropna().unique()]
    colors = sns.color_palette("husl", len(peak_aliases)) if peak_aliases else []
    color_map = {alias: colors[i] for i, alias in enumerate(peak_aliases)}

    for alias in peak_aliases:
        peak_vals = all_itd_insertions.loc[all_itd_insertions["peak_alias"] == alias, "ins_len"]
        if peak_vals.empty:
            continue

        mean_len = peak_vals.mean()
        sd_len = peak_vals.std()
        c = color_map[alias]

        plt.axvline(
            mean_len,
            color=c,
            linestyle="--",
            lw=1.5,
            alpha=0.95,
            label=f"{alias} mean={mean_len:.1f} (n={len(peak_vals)})"
        )
        if pd.notna(sd_len) and sd_len > 0:
            plt.axvspan(
                mean_len - sd_len,
                mean_len + sd_len,
                alpha=0.12,
                color=c,
                label=f"{alias} +/-SD={sd_len:.1f}"
            )

    # De-duplicate legend entries from histogram + overlays
    handles, labels = plt.gca().get_legend_handles_labels()
    dedup = {}
    for h, l in zip(handles, labels):
        if l not in dedup:
            dedup[l] = h
    plt.legend(dedup.values(), dedup.keys())
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
        Amplicon coordinates are 1-based inclusive genomic positions.
    exon_boundaries : list[(start, end)]
        Exon coordinates for plotting context.
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

        # Convert local insertion boundary to a 1-based genomic anchor coordinate.
        if ins_pos_local < 1:
            logger.warning(
                f"[build_itd_reference_per_peak] {alias}: local insertion position {ins_pos_local} "
                f"is <1; clamping to 1 for genomic anchoring."
            )
            ins_pos_local = 1
        elif ins_pos_local > len(ref_seq):
            logger.warning(
                f"[build_itd_reference_per_peak] {alias}: local insertion position {ins_pos_local} "
                f"exceeds reference length {len(ref_seq)}; clamping to {len(ref_seq)}."
            )
            ins_pos_local = len(ref_seq)
        genomic_ins_pos = amp["start"] + ins_pos_local - 1
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
        Coordinates of the amplicon in 1-based genomic positions.
    exon_boundaries : list[(start, end)]
        List of exon intervals in genomic positions.
    genome_build : str
        "hg19" or "hg38".
    insertion_genomic_pos : int
        1-based genomic anchor coordinate of the ITD insertion site.
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

    # Already a 1-based genomic anchor coordinate
    ins_pos_display = insertion_genomic_pos

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
        "local_insertion_pos": None,
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
            "local_insertion_pos": entry.get("local_insertion_pos"),
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
        'read_id', 'read_seq', 'strand', 'read_len', 'gmm_peak_alias', 'is_ambiguous'.

    Returns
    -------
    filtered_df : pd.DataFrame
        Cleaned DataFrame containing only reads to use for validation.
    """

    # --- Sanity checks ---
    required_cols = {"read_id", "read_seq", "strand", "gmm_peak_alias", "is_ambiguous"}
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

    # Keep only explicit support calls for downstream summary metrics.
    df_support = df_valid[df_valid["support_call"].isin(["ITD-supporting", "WT-supporting"])].copy()
    df_itd_support = df_support[df_support["support_call"] == "ITD-supporting"].copy()
    df_wt_support = df_support[df_support["support_call"] == "WT-supporting"].copy()

    # Count ITD-supporting and WT-supporting reads
    itd_counts = (
        df_itd_support
        .groupby("ref_alias")
        .agg(n_itd_reads=("read_id", "nunique"))
        .reset_index() 
    )

    # Depth is total validated/classified reads used in support calls.
    itd_counts["n_total_reads"] = len(df_support)
    itd_counts["allele_frequency"] = itd_counts["n_itd_reads"] / itd_counts["n_total_reads"]

    # Strand support counts are ITD-supporting reads only.

    strand_counts = (
        df_itd_support.groupby(["ref_alias", "strand"])
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
    wt_plus = int((df_wt_support["strand"] == "+").sum())
    wt_minus = int((df_wt_support["strand"] == "-").sum())

    for _, row in strand_counts.iterrows():
        alias = str(row["ref_alias"])
        plus = int(row["+"])
        minus = int(row["-"])
        total = plus + minus
        if total == 0:
            logger.warning(f"No reads found for {alias}, skipping strand bias calculation.")
            continue

        ratio = plus / (total + 1e-8)
        # Fisher's exact test: ITD strand split vs WT-supporting strand split.
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
    if strand_stats.empty:
        strand_stats = pd.DataFrame(
            columns=["ref_alias", "plus_reads", "minus_reads", "strand_ratio", "odds_ratio", "fisher_p"]
        )

    logger.info("Applying allele frequency filter to detected ITDs...")
    logger.info(f"ITD counts before applying min allele frequency of {min_allele_frequency}:")
    logger.info(itd_counts)

    itd_counts = itd_counts[itd_counts["allele_frequency"] >= min_allele_frequency].copy()
    if itd_counts.empty:
        itd_counts = pd.DataFrame(columns=["ref_alias", "n_itd_reads", "n_total_reads", "allele_frequency"])

    strand_stats["ref_alias"] = strand_stats["ref_alias"].astype(str)
    itd_counts = itd_counts.reset_index(drop=True)
    itd_counts["ref_alias"] = itd_counts["ref_alias"].astype(str)

    logger.info(f"Detected ITDs after applying min allele frequency of {min_allele_frequency}:")
    logger.info(itd_counts)

    summary_df = itd_counts.merge(strand_stats, on="ref_alias", how="left")
    summary_df = summary_df.sort_values("allele_frequency", ascending=False)

    logger.debug("printing summary_df")
    logger.debug(type(summary_df))
    logger.debug(summary_df)

    return summary_df

