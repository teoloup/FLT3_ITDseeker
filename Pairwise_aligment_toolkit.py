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

def build_validation_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.target_open_gap_score = -15
    aligner.target_extend_gap_score = -0.01
    aligner.query_open_gap_score = -200
    aligner.query_extend_gap_score = -2
    aligner.query_end_gap_score = -200
    aligner.target_end_gap_score = -200
    aligner.mode = "global"
    return aligner

def find_insertions_between_blocks(aln):
    """Detect insertions (in query) as large jumps between aligned blocks."""
    t_blocks, q_blocks = aln.aligned
    insertions = []

    for i in range(len(t_blocks) - 1):
        t_end = t_blocks[i][1]
        t_next = t_blocks[i + 1][0]
        q_end = q_blocks[i][1]
        q_next = q_blocks[i + 1][0]

        # If query advanced more than target between blocks → insertion
        if (q_next - q_end) > (t_next - t_end):
            ins_len = (q_next - q_end) - (t_next - t_end)
            insertions.append((t_end, q_end, ins_len))
    return insertions

def chunk_iterable(data, n):
    """Split list into n roughly equal chunks."""
    k = max(1, len(data) // n)
    for i in range(0, len(data), k):
        yield data[i:i + k]

def aligned_length(aln):
    # sum of lengths of aligned blocks on either sequence (they're equal)
    return sum(e - s for s, e in aln.aligned[0])

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
    ins = find_insertions_between_blocks(aln)
    if not ins:
        # fallback to within-block method (very short indels)
        t_blocks, q_blocks = aln.aligned
        for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
            if (qe - qs) > (te - ts):
                ins.append((ts, qs, qe - qs - (te - ts)))
    return ins

def best_orientation_alignment(
    ref_seq: str,
    read_seq_str: str,
    aligner: Align.PairwiseAligner
):
    """
    Align both forward and reverse-complement orientations of a read
    against a reference and return the best alignment.

    Returns
    -------
    aln : Bio.Align.Alignment
        Best alignment object
    strand : str
        '+' if forward, '-' if reverse complement
    used_read : str
        The sequence used for alignment (in correct orientation)
    score_f : float
        Forward alignment score
    score_r : float
        Reverse alignment score
    """
    # --- Forward alignment ---
    aln_f = aligner.align(ref_seq, read_seq_str)[0]
    score_f = aln_f.score

    # --- Reverse complement alignment ---
    rc = str(Seq(read_seq_str).reverse_complement())
    aln_r = aligner.align(ref_seq, rc)[0]
    score_r = aln_r.score

    # --- Select best ---
    if score_f >= score_r:
        return aln_f, "+", read_seq_str, score_f, score_r
    else:
        return aln_r, "-", rc, score_f, score_r

def process_chunk(chunk, itd_min, itd_max, ref_seq, peak_alias):
        aligner = build_default_aligner()
        out_rows = []
        for read_id, read_seq in chunk:
            try:
                aln, strand, used_seq, f_score, r_score = best_orientation_alignment(ref_seq, read_seq, aligner)
                pid = percent_identity(aln)
                ins_regions = find_insertions(aln)

                for ts, qs, ins_len in ins_regions:
                    if itd_min <= ins_len <= itd_max:
                        ins_seq = used_seq[qs: qs + ins_len]
                        out_rows.append({
                            "peak_alias": peak_alias,
                            "read_id": read_id,
                            "strand": strand,
                            "aln_score": aln.score,
                            "pct_identity": pid,
                            "ins_pos_ref": ts,
                            "ins_len": ins_len,
                            "ins_seq": ins_seq,
                            "fwd_score": f_score,
                            "rev_score": r_score
                        })
            except Exception as e:
                logger.warning(f"[WARN] Alignment failed for {read_id}: {e}")
            # logger.debug(f"Processed read {read_id}: found {len(ins_regions)} insertions. Insertions in ITD range: {[ins for ins in ins_regions if itd_min <= ins[2] <= itd_max]}")
            # logger.debug(f"Alignment details for {read_id}: strand={strand}, fwd_score={f_score}, rev_score={r_score}, aln_score={aln.score}, pct_identity={pid:.2%}")
            # logger.debug(f"Alignment:\n{aln}")
        return out_rows

def _align_reads_chunk(reads_chunk, ref_dict, aligner, alpha=1):
    """
    Align a chunk of reads to all reference sequences (WT + ITDs),
    computing both raw and adjusted scores and selecting best strand.
    """
    results = []
    logger.debug(f"[_align_reads_chunk] Starting chunk of {len(reads_chunk)} reads...")

    for rid, read_seq in reads_chunk:
        if not read_seq or not isinstance(read_seq, str):
            logger.warning(f"[_align_reads_chunk] Skipping empty or invalid read {rid}")
            continue

        strand_data = {}

        for strand, seq in [("+", read_seq), ("-", str(Seq(read_seq).reverse_complement()))]:
            strand_alignments = []
            strand_scores = []

            for alias, ref_entry in ref_dict.items():
                ref_seq = ref_entry.get("ref_seq_with_itd", "")
                if not ref_seq:
                    logger.error(f"[_align_reads_chunk] Empty reference for {alias}")
                    aln_obj = None
                    score = 0.0
                else:
                    alns = aligner.align(ref_seq, seq)
                    if len(alns) == 0:
                        aln_obj = None
                        score = 0.0
                    else:
                        aln_obj = alns[0]
                        score = aln_obj.score

                strand_alignments.append((alias, aln_obj, score))
                strand_scores.append(score)

            # Choose max score instead of mean
            strand_score = np.max(strand_scores) if strand_scores else 0.0
            strand_data[strand] = (strand_score, strand_alignments)

        # Strand decision
        if abs(strand_data["+"][0] - strand_data["-"][0]) < 5:
            logger.debug(
                f"Read {rid} strand ambiguous: +={strand_data['+'][0]:.1f}, -={strand_data['-'][0]:.1f}"
            )

        best_strand, (_, best_alignments) = max(strand_data.items(), key=lambda kv: kv[1][0])

        # Record results
        for alias, aln, raw_score in best_alignments:
            pid = percent_identity(aln) * 100 if aln else 0.0
            adjusted_score = compute_adjusted_score(aln, alpha)
            aligned_blocks = aln.aligned if aln else []
            results.append({
                "read_id": rid,
                "ref_alias": alias,
                "strand": best_strand,
                "score": raw_score,
                "adjusted_score": adjusted_score,
                "pct_identity": pid,
                "aligned_blocks": aligned_blocks,
            })

        if len(results) % 1000 == 0:
            logger.debug(f"[_align_reads_chunk] Processed {len(results)} alignments so far...")

    logger.info(f"[_align_reads_chunk] Completed chunk ({len(reads_chunk)} reads → {len(results)} alignments).")

    # Diagnostics
    if results and logger.isEnabledFor(logging.DEBUG):
        scores = [r["score"] for r in results]
        adj_scores = [r["adjusted_score"] for r in results]
        logger.debug(
            f"[_align_reads_chunk] Raw score stats: min={np.min(scores):.2f}, "
            f"max={np.max(scores):.2f}, mean={np.mean(scores):.2f}"
        )
        logger.debug(
            f"[_align_reads_chunk] Adjusted score stats: min={np.min(adj_scores):.2f}, "
            f"max={np.max(adj_scores):.2f}, mean={np.mean(adj_scores):.2f}"
        )

    return results

def align_reads_multi_ref_parallel(reads_df, ref_dict, df_cons, logger, threads=8, beta=500.0):
    """
    Parallel alignment of reads to WT + ITD references with softmax-normalized adjusted scores.
    Supports adaptive metric selection for 2 vs. ≥3 references.
    """
    aligner = build_validation_aligner()
    reads_list = list(zip(reads_df["read_id"], reads_df["read_seq"]))
    chunk_size = int(np.ceil(len(reads_list) / threads))

    logger.info(
        f"Running parallel alignment on {len(reads_list)} reads across {len(ref_dict)} references "
        f"using {threads} threads."
    )

    # --- Parallel alignment ---
    all_results = []
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futures = [
            ex.submit(_align_reads_chunk, reads_list[i:i + chunk_size], ref_dict, aligner)
            for i in range(0, len(reads_list), chunk_size)
        ]
        for f in as_completed(futures):
            res = f.result()
            if res:
                all_results.extend(res)

    df_results = pd.DataFrame(all_results)
    logger.info(f"[align_reads_multi_ref_parallel] Collected {len(df_results)} total alignments.")

    if df_results.empty:
        logger.error("No alignments produced! Check read sequences or reference dictionary.")
        return pd.DataFrame()

    # --- Use adjusted_score instead of raw score ---
    df_results["score_used"] = df_results["adjusted_score"].fillna(0.0)

    # --- Softmax normalization per read ---
    def softmax_transform(scores):
        return softmax(scores, beta=beta)

    df_results["prob_score"] = df_results.groupby("read_id")["score_used"].transform(softmax_transform)

    n_refs = df_results["ref_alias"].nunique()
    logger.info(f"[align_reads_multi_ref_parallel] Detected {n_refs} reference sequences.")

    # --- Metric selection ---
    if n_refs == 2:
        # ==========================
        # CASE A — Binary WT vs ITD
        # ==========================
        logger.info("Using probability delta metric (WT vs ITD).")

        df_results["metric_value"] = df_results["prob_score"]
        df_results["metric_type"] = "prob_delta"
        metric_used = "prob_delta"

        # Only need top 2 per read (WT and ITD)
        df_results = df_results.sort_values(["read_id", "metric_value"], ascending=[True, False])

        top_hits = (
            df_results
            .groupby("read_id", group_keys=False)
            .head(2)
            .reset_index(drop=True)
        )

        df_best = top_hits.groupby("read_id").nth(0).reset_index()
        df_second = top_hits.groupby("read_id").nth(1).reset_index()

        df_best = df_best.merge(
            df_second[["read_id", "ref_alias", "metric_value"]],
            on="read_id",
            suffixes=("", "_second"),
            how="left"
        )

    else:
        # ==========================
        # CASE B — Multi-ITD (>2 refs)
        # ==========================
        logger.info("Using z-score metric for multi-ITD comparison.")

        df_results["metric_value"] = (
            df_results.groupby(["read_id", "strand"])["prob_score"]
            .transform(lambda x: (x - x.mean()) / (x.std(ddof=0) + 1e-8))
        )
        df_results["metric_type"] = "z_score"
        metric_used = "z_score"

        # Keep all refs for meaningful z-scores
        df_best = (
            df_results.loc[df_results.groupby("read_id")["metric_value"].idxmax()]
            .reset_index(drop=True)
        )

        # Find the 2nd best for delta comparison
        second_best = (
            df_results
            .sort_values(["read_id", "metric_value"], ascending=[True, False])
            .groupby("read_id", group_keys=False)
            .nth(1)
            .reset_index()
        )

        df_best = df_best.merge(
            second_best[["read_id", "ref_alias", "metric_value"]],
            on="read_id",
            suffixes=("", "_second"),
            how="left"
        )

    # --- Compute delta between best and 2nd best ---
    df_best["second_best_ref"] = df_best["ref_alias_second"]
    df_best["second_best_score"] = df_best["metric_value_second"]
    df_best["delta"] = df_best["metric_value"] - df_best["metric_value_second"]
    df_best.drop(columns=["ref_alias_second", "metric_value_second"], inplace=True, errors="ignore")

    # --- Classification ---
    df_best["support_call"] = df_best.apply(
        lambda r: classify_read_support(r, metric_used, z_thresh=1, delta_thresh=0.05),
        axis=1
    )

    logger.info(f"Metric used: {metric_used}, refs={n_refs}")
    logger.info(f"Z-score range: {df_results['metric_value'].min():.3f}–{df_results['metric_value'].max():.3f}")


    # --- Validate ITD-supporting reads ---
    logger.info("[align_reads_multi_ref_parallel] Validating ITD-supporting reads by gap inspection.")
    df_best = validate_itd_supporting_reads(
        df_best,
        df_cons,
        gap_window=15,
        max_gap_bp=5,
        min_pid=0.9,
    )

    # --- Summary ---
    logger.info(f"[align_reads_multi_ref_parallel] Validation complete using {metric_used}.")
    logger.debug(
        f"[align_reads_multi_ref_parallel] Support summary:\n"
        f"{df_best['support_call'].value_counts(dropna=False).to_string()}"
    )

    return df_best

def compute_adjusted_score(aln, alpha=0.5):
    """Compute hybrid PID + gap-penalized score from a Biopython alignment."""
    if aln is None:
        return 0.0

    try:
        # --- Percent identity (0–1) ---
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
        # The z-score itself already encodes "confidence" — how extreme this ref is vs all others.
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
    min_pid
):
    """
    Validate both WT and ITD reads by inspecting precomputed alignments.

    WT validation:
      - Must have high PID
      - No large insertions (total gap length > 10 bp → fail)

    ITD validation:
      - Must have PID ≥ min_pid
      - Must show an insertion near expected breakpoint (±gap_window)
      - Total gap bases near insertion ≤ max_gap_bp
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

            # Compute total gaps (insertions relative to ref)
            total_gap_len = sum(
                max(0, (q_end - q_start) - (t_end - t_start))
                for (t_start, t_end), (q_start, q_end) in zip(*aln_blocks)
            )

            if pid < min_pid:
                validated.append(False)
                reasons.append(f"low_pid({pid:.2f})")
            elif total_gap_len > 10:
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

            ins_pos = ins_row["median_ins_pos_ref"].iloc[0]

            if aln_blocks is None or len(aln_blocks) != 2:
                validated.append(False)
                reasons.append("missing_alignment_blocks")
                continue

            # Check for insertions near ITD breakpoint
            gaps_near = 0
            for (t_start, t_end), (q_start, q_end) in zip(*aln_blocks):
                if (q_end - q_start) > (t_end - t_start):
                    gap_pos = (t_start + t_end) / 2
                    if abs(gap_pos - ins_pos) <= gap_window:
                        gaps_near += (q_end - q_start) - (t_end - t_start)

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
    """Count total gap bases in the reference ±window bp around the insertion site."""
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
