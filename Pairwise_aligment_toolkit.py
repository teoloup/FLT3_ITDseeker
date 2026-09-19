import os
os.environ["MPLBACKEND"] = "Agg"       # disable any GUI backend
os.environ["DISPLAY"] = ""             # make sure Tk can't open a window
os.environ["TK_SILENCE_DEPRECATION"] = "1"
import numpy as np
import pandas as pd
import logging
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import Align
from Helper_functions import (
    percent_identity,
    span_identity,
    compute_adjusted_score,
    softmax,
    classify_read_support,
    validate_itd_supporting_reads,
)

logger = logging.getLogger(__name__)


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

def _align_reads_chunk(reads_chunk, ref_dict, alpha=2.0):
    """
    Align a chunk of reads to all reference sequences (WT + ITDs),
    computing both raw and adjusted scores on the stored read orientation.
    """
    aligner = build_validation_aligner()
    results = []
    logger.debug(f"[_align_reads_chunk] Starting chunk of {len(reads_chunk)} reads...")

    for rid, read_seq, read_strand in reads_chunk:
        if not read_seq or not isinstance(read_seq, str):
            logger.warning(f"[_align_reads_chunk] Skipping empty or invalid read {rid}")
            continue

        best_alignments = []
        for alias, ref_entry in ref_dict.items():
            ref_seq = ref_entry.get("ref_seq_with_itd", "")
            if not ref_seq:
                logger.error(f"[_align_reads_chunk] Empty reference for {alias}")
                aln_obj = None
                score = 0.0
            else:
                alns = aligner.align(ref_seq, read_seq)
                if len(alns) == 0:
                    aln_obj = None
                    score = 0.0
                else:
                    aln_obj = alns[0]
                    score = aln_obj.score
            best_alignments.append((alias, aln_obj, score))

        # Record results
        for alias, aln, raw_score in best_alignments:
            # Two different identities, for two different jobs:
            #   pct_identity (block-based) answers "is this read good enough?" and
            #     feeds the min_pid quality filter;
            #   span identity answers "which reference does it fit best?" and feeds
            #     the competitive score. Block identity cannot do the second job --
            #     it reads 1.0 against every reference for a clean read.
            pid = percent_identity(aln) if aln else 0.0
            span_pid = span_identity(aln) if aln else 0.0
            adjusted_score = compute_adjusted_score(aln, alpha, pid=span_pid)
            aligned_blocks = aln.aligned if aln else []
            results.append({
                "read_id": rid,
                "ref_alias": alias,
                "strand": read_strand,
                "score": raw_score,
                "adjusted_score": adjusted_score,
                "pct_identity": pid,
                "span_identity": span_pid,
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
    reads_list = list(zip(reads_df["read_id"], reads_df["read_seq"], reads_df["strand"]))
    if threads < 1:
        threads = 1
    # Several chunks per worker rather than exactly one: with one batch each,
    # a single slow worker stalls the whole pass and nothing can be stolen
    # from it. Smaller batches let the pool rebalance; total IPC is unchanged.
    chunks_per_worker = 4
    chunk_size = max(1, int(np.ceil(len(reads_list) / (threads * chunks_per_worker))))
    batches = [reads_list[i:i + chunk_size] for i in range(0, len(reads_list), chunk_size)]
    n_batches = len(batches)
    logger.debug(
        "[align_reads_multi_ref_parallel] %d batches of <=%d reads across %d workers",
        n_batches, chunk_size, threads,
    )
    n_reads = len(reads_list)
    n_refs = len(ref_dict)
    est_comparisons = n_reads * n_refs

    logger.info(
        "[align_reads_multi_ref_parallel] Starting pairwise alignment: "
        f"reads={n_reads}, refs={n_refs}, batches={n_batches}, workers={threads}, "
        f"estimated_comparisons={est_comparisons}"
    )

    # --- Parallel alignment ---
    t0 = time.perf_counter()
    all_results = []
    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [
            ex.submit(_align_reads_chunk, batch, ref_dict)
            for batch in batches
        ]
        for f in as_completed(futures):
            res = f.result()
            if res:
                all_results.extend(res)
    elapsed_sec = time.perf_counter() - t0

    reads_per_sec = (n_reads / elapsed_sec) if elapsed_sec > 0 else 0.0
    comps_per_sec = (est_comparisons / elapsed_sec) if elapsed_sec > 0 else 0.0
    logger.info(
        "[align_reads_multi_ref_parallel] Finished pairwise alignment in "
        f"{elapsed_sec:.2f}s ({reads_per_sec:.1f} reads/s, {comps_per_sec:.1f} est comparisons/s)"
    )

    df_results = pd.DataFrame(all_results)
    logger.info(f"[align_reads_multi_ref_parallel] Collected {len(df_results)} total alignments.")

    # Workers complete in nondeterministic order, so row order here varies
    # between runs. Pin it before any ranking, otherwise exact score ties are
    # broken by whichever worker happened to finish first.
    df_results = df_results.sort_values(
        ["read_id", "ref_alias"], kind="mergesort"
    ).reset_index(drop=True)

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

    # A read supports whichever reference scores highest. Ties break towards WT
    # (a coin flip should not become a variant call), then alphabetically, so the
    # outcome never depends on which worker finished first.
    df_results["_wt_first"] = (df_results["ref_alias"] == "WT").astype(int)

    # --- Metric selection ---
    if n_refs == 2:
        # ==========================
        # CASE A - Binary WT vs ITD
        # ==========================
        logger.info("Using probability delta metric (WT vs ITD).")

        df_results["metric_value"] = df_results["prob_score"]
        df_results["metric_type"] = "prob_delta"
        metric_used = "prob_delta"

    else:
        # ==========================
        # CASE B - Multi-ITD (>2 refs)
        # ==========================
        logger.info("Using z-score metric for multi-ITD comparison.")

        df_results["metric_value"] = (
            df_results.groupby("read_id")["prob_score"]
            .transform(lambda x: (x - x.mean()) / (x.std(ddof=0) + 1e-8))
        )
        df_results["metric_type"] = "z_score"
        metric_used = "z_score"

    # --- Rank references per read, best first ---
    ranked = df_results.sort_values(
        ["read_id", "metric_value", "_wt_first", "ref_alias"],
        ascending=[True, False, False, True],
        kind="mergesort",
    )
    ranked["_rank"] = ranked.groupby("read_id").cumcount()

    # cumcount rather than GroupBy.nth: nth changed between pandas 1.x and 2.x
    # (index vs filter semantics); this behaves the same on both.
    df_best = ranked[ranked["_rank"] == 0].reset_index(drop=True)
    df_second = ranked[ranked["_rank"] == 1].reset_index(drop=True)

    df_best = df_best.merge(
        df_second[["read_id", "ref_alias", "metric_value"]],
        on="read_id",
        suffixes=("", "_second"),
        how="left",
    )
    df_best = df_best.drop(columns=["_wt_first", "_rank"], errors="ignore")

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

    # aligned_blocks holds raw numpy index arrays used only by the validation
    # step above. Left in place it is written verbatim into
    # *_validation_read_support.tsv as array reprs (tens of MB of unreadable
    # text on a deep run), so drop it once validation has consumed it.
    df_best = df_best.drop(columns=["aligned_blocks"], errors="ignore")

    # --- Summary ---
    logger.info(f"[align_reads_multi_ref_parallel] Validation complete using {metric_used}.")
    logger.debug(
        f"[align_reads_multi_ref_parallel] Support summary:\n"
        f"{df_best['support_call'].value_counts(dropna=False).to_string()}"
    )

    return df_best

