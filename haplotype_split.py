"""
Per-peak haplotype splitting.

`fit_gmm_itds` groups reads by length. Two ITDs of the same length land in the
same peak and no length-based method can separate them, so the per-peak MSA
mixes haplotypes and the consensus fills with `N`. This module replaces that
second pass with pluggable per-peak *sequence* clustering.

Every backend has the same contract as `refine_peak_substructure_once`: it takes
`(comps, reads_df, peak_subsets)` and returns a `PeakRefineResult` with the same
three fields, so callers downstream are unaffected.

Backends receive the full primer-trimmed reads assigned to one peak, with their
base qualities. Reads within a peak already share a length and are already
trimmed, so each tool sees the kind of input it was built for. Peaks are
clustered independently.
"""

import logging
import os
import shutil
import subprocess
import tempfile
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from GMM_peaks import PeakRefineResult, refine_peak_substructure_once

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def write_peak_fastq(reads_df: pd.DataFrame, read_ids: List[str], path: str) -> int:
    """Write one peak's reads to FASTQ. Returns the number written.

    Reads without stored qualities get a flat placeholder so the file stays
    valid FASTQ; a backend that genuinely needs real qualities should check
    `reads_df["read_qual"]` itself rather than trusting the placeholder.
    """
    subset = reads_df.loc[reads_df["read_id"].isin(read_ids)]
    n = 0
    with open(path, "w", newline="\n") as fh:
        for row in subset.itertuples(index=False):
            seq = row.read_seq
            if not seq:
                continue
            qual = getattr(row, "read_qual", "") or ""
            if len(qual) != len(seq):
                qual = "I" * len(seq)  # Q40 placeholder
            fh.write(f"@{row.read_id}\n{seq}\n+\n{qual}\n")
            n += 1
    return n


def _require_tool(binary: str, backend: str) -> str:
    """Resolve an external tool or fail loudly. Never silently fall back."""
    found = shutil.which(binary)
    if not found:
        raise RuntimeError(
            f"Haplotype backend '{backend}' needs '{binary}' on PATH and it was not "
            f"found. Install it (see the Dockerfile) or choose a different "
            f"--haplotype-method."
        )
    return found


def assignments_to_result(
    *,
    comps: pd.DataFrame,
    reads_df: pd.DataFrame,
    peak_subsets: Dict[str, List[str]],
    assignments: Dict[str, Dict[str, str]],
    wt_amplicon_length: float,
    min_child_fraction: float,
    min_child_reads: int,
) -> PeakRefineResult:
    """Turn per-peak read->cluster labels into the standard refinement result.

    `assignments` maps peak_alias -> {read_id: cluster_label}. A peak absent from
    it is left untouched.

    Guardrails mirror the existing length-based refinement: a cluster must hold
    at least `min_child_fraction` of the parent's reads and at least
    `min_child_reads` reads to become a haplotype. If fewer than two clusters
    survive, the peak is left unsplit.

    Reads in sub-threshold clusters are reassigned to the largest surviving
    cluster rather than dropped. Dropping them would shrink the allele-frequency
    denominator and inflate every other call; they are most likely error-driven
    splits off the dominant haplotype, so folding them back is also the
    conservative choice for ITD count.
    """
    reads_out = reads_df.copy()
    subsets_out = {k: list(v) for k, v in peak_subsets.items()}

    len_by_read = dict(zip(reads_df["read_id"], reads_df["read_len"]))

    wt_rows = comps.loc[comps["peak_alias"].astype(str).str.upper() == "WT"]
    if wt_rows.empty:
        wt_mean = float(
            comps.loc[(comps["mean_bp"] - wt_amplicon_length).abs().idxmin(), "mean_bp"]
        )
    else:
        wt_mean = float(wt_rows.iloc[0]["mean_bp"])

    split_children: Dict[str, List[Dict]] = {}

    for parent_alias, read_to_cluster in assignments.items():
        parent_reads = subsets_out.get(parent_alias, [])
        if not parent_reads or not read_to_cluster:
            continue

        clusters: Dict[str, List[str]] = {}
        for rid in parent_reads:
            lbl = read_to_cluster.get(rid)
            if lbl is None:
                continue
            clusters.setdefault(str(lbl), []).append(rid)
        if not clusters:
            continue

        n_parent = sum(len(v) for v in clusters.values())
        ordered = sorted(clusters.items(), key=lambda kv: (-len(kv[1]), kv[0]))
        keep = [
            (lbl, rids) for lbl, rids in ordered
            if len(rids) >= min_child_reads and len(rids) / n_parent >= min_child_fraction
        ]

        if len(keep) < 2:
            logger.info(
                "[haplotype_split] %s: %d cluster(s) from %d reads, %d passed the "
                "size guardrails -- leaving the peak unsplit.",
                parent_alias, len(clusters), n_parent, len(keep),
            )
            continue

        # fold sub-threshold clusters into the largest survivor
        kept_labels = {lbl for lbl, _ in keep}
        folded = 0
        biggest = keep[0][0]
        merged: Dict[str, List[str]] = {lbl: list(rids) for lbl, rids in keep}
        for lbl, rids in ordered:
            if lbl not in kept_labels:
                merged[biggest].extend(rids)
                folded += len(rids)
        if folded:
            logger.info(
                "[haplotype_split] %s: folded %d reads from sub-threshold clusters "
                "into %s_H1.", parent_alias, folded, parent_alias,
            )

        children = []
        for i, (lbl, rids) in enumerate(
            sorted(merged.items(), key=lambda kv: (-len(kv[1]), kv[0])), start=1
        ):
            lengths = np.array([len_by_read[r] for r in rids if r in len_by_read],
                               dtype=float)
            child_alias = f"{parent_alias}_H{i}"
            children.append({
                "alias": child_alias,
                "read_ids": rids,
                "mean_bp": float(lengths.mean()) if lengths.size else float("nan"),
                "sd_bp": float(lengths.std(ddof=0)) if lengths.size else 0.0,
                "cluster_label": lbl,
            })
            logger.info(
                "[haplotype_split] %s -> %s (n=%d, mean=%.1f bp, sd=%.2f bp)",
                parent_alias, child_alias, len(rids),
                children[-1]["mean_bp"], children[-1]["sd_bp"],
            )

        split_children[parent_alias] = children
        subsets_out.pop(parent_alias, None)
        for ch in children:
            subsets_out[ch["alias"]] = ch["read_ids"]
            reads_out.loc[
                reads_out["read_id"].isin(ch["read_ids"]), "gmm_peak_alias"
            ] = ch["alias"]

    if not split_children:
        return PeakRefineResult(comps=comps, reads_df=reads_out, peak_subsets=subsets_out)

    eff_counts = {alias: len(ids) for alias, ids in subsets_out.items()}
    total_eff = sum(eff_counts.values()) or 1

    rows = []
    for _, row in comps.iterrows():
        alias = str(row["peak_alias"])
        if alias in split_children:
            for ch in split_children[alias]:
                count = int(eff_counts.get(ch["alias"], 0))
                frac = count / total_eff
                rows.append({
                    "mean_bp": ch["mean_bp"],
                    "sd_bp": ch["sd_bp"],
                    "fraction": frac,
                    "read_count": count,
                    "effective_read_count": count,
                    "effective_allele_freq": frac,
                    "putative_itd_size": ch["mean_bp"] - wt_mean,
                    "is_wt": False,
                    "peak_alias": ch["alias"],
                    "parent_peak_alias": alias,
                    "refinement_level": 1,
                    "is_refined_child": True,
                })
        else:
            count = int(eff_counts.get(alias, 0))
            frac = count / total_eff
            rows.append({
                "mean_bp": float(row["mean_bp"]),
                "sd_bp": float(row["sd_bp"]),
                "fraction": frac,
                "read_count": count,
                "effective_read_count": count,
                "effective_allele_freq": frac,
                "putative_itd_size": float(row["mean_bp"]) - wt_mean,
                "is_wt": bool(alias.upper() == "WT"),
                "peak_alias": alias,
                "parent_peak_alias": alias,
                "refinement_level": int(row.get("refinement_level", 0)),
                "is_refined_child": bool(row.get("is_refined_child", False)),
            })

    comps_out = (
        pd.DataFrame(rows).sort_values("fraction", ascending=False).reset_index(drop=True)
    )
    logger.info(
        "[haplotype_split] %d -> %d peaks after haplotype splitting.",
        len(comps), len(comps_out),
    )
    return PeakRefineResult(comps=comps_out, reads_df=reads_out, peak_subsets=subsets_out)


def _target_peaks(comps: pd.DataFrame, cluster_wt_peak: bool) -> List[str]:
    aliases = [str(a) for a in comps["peak_alias"]]
    if cluster_wt_peak:
        return aliases
    return [a for a in aliases if a.upper() != "WT"]


# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------

def _backend_none(**kw) -> PeakRefineResult:
    """No splitting at all: the first-pass GMM peaks are the final peaks."""
    logger.info("[haplotype_split] method=none -- keeping first-pass GMM peaks.")
    return PeakRefineResult(
        comps=kw["comps"], reads_df=kw["reads_df"], peak_subsets=kw["peak_subsets"]
    )


def _backend_gmm2pass(**kw) -> PeakRefineResult:
    """Baseline: the existing one-level length-based refinement, unchanged."""
    opts = kw.get("gmm_kwargs") or {}
    return refine_peak_substructure_once(
        comps=kw["comps"],
        reads_df=kw["reads_df"],
        peak_subsets=kw["peak_subsets"],
        wt_amplicon_length=kw["wt_amplicon_length"],
        **opts,
    )


def _backend_isonclust(**kw) -> PeakRefineResult:
    """Cluster each peak's reads with isONclust.

    isONclust is quality-aware and takes FASTQ directly, writing
    `final_clusters.tsv` as (cluster_id, read_id). Its `--ont` preset (k=13,
    w=20) is tuned to separate genes; ITD haplotypes differ far less than that,
    so k/w are exposed for tuning.
    """
    binary = _require_tool("isONclust", "isonclust")
    comps, reads_df = kw["comps"], kw["reads_df"]
    peak_subsets = kw["peak_subsets"]
    opts = kw.get("tool_kwargs") or {}
    k = int(opts.get("k", 13))
    w = int(opts.get("w", 20))
    # isONclust defaults aligned_threshold to 0.4 and mapped_threshold to 0.7,
    # which are right for deciding whether two reads came from the same *gene*.
    # ITD haplotypes differ over a few dozen bases of an otherwise identical
    # read, so at those defaults every haplotype in a peak joins one cluster.
    # Both are raised here and exposed for tuning. 0.80 was picked by sweeping
    # against simulated data with known haplotypes: 0.4 and 0.6 return a single
    # cluster (no separation at all), 0.8 gives ARI 0.68 at 89% cluster purity,
    # and 0.95+ fragments into hundreds of clusters (ARI 0.25 and falling).
    aligned_threshold = float(opts.get("aligned_threshold", 0.80))
    mapped_threshold = float(opts.get("mapped_threshold", 0.90))
    threads = int(kw.get("threads", 1))
    work_dir = kw["work_dir"]

    assignments: Dict[str, Dict[str, str]] = {}
    for alias in _target_peaks(comps, kw.get("cluster_wt_peak", False)):
        read_ids = peak_subsets.get(alias, [])
        if len(read_ids) < kw["min_child_reads"] * 2:
            continue

        peak_dir = os.path.join(work_dir, f"isonclust_{alias}")
        os.makedirs(peak_dir, exist_ok=True)
        fq = os.path.join(peak_dir, "reads.fastq")
        n = write_peak_fastq(reads_df, read_ids, fq)
        logger.info(
            "[haplotype_split] isONclust on %s (%d reads, k=%d w=%d "
            "aligned_threshold=%.2f mapped_threshold=%.2f)",
            alias, n, k, w, aligned_threshold, mapped_threshold,
        )

        out_dir = os.path.join(peak_dir, "out")
        cmd = [binary, "--fastq", fq, "--outfolder", out_dir,
               "--k", str(k), "--w", str(w), "--t", str(threads),
               "--aligned_threshold", str(aligned_threshold),
               "--mapped_threshold", str(mapped_threshold)]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if res.returncode != 0:
            raise RuntimeError(
                f"isONclust failed on peak {alias}: {res.stderr.decode()[-2000:]}"
            )

        tsv = os.path.join(out_dir, "final_clusters.tsv")
        if not os.path.exists(tsv):
            raise RuntimeError(
                f"isONclust produced no final_clusters.tsv for peak {alias} "
                f"(looked in {out_dir})"
            )
        mapping: Dict[str, str] = {}
        with open(tsv) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    mapping[parts[1]] = parts[0]
        assignments[alias] = mapping

    return assignments_to_result(
        comps=comps, reads_df=reads_df, peak_subsets=peak_subsets,
        assignments=assignments,
        wt_amplicon_length=kw["wt_amplicon_length"],
        min_child_fraction=kw["min_child_fraction"],
        min_child_reads=kw["min_child_reads"],
    )


BACKENDS: Dict[str, Callable[..., PeakRefineResult]] = {
    "none": _backend_none,
    "gmm2pass": _backend_gmm2pass,
    "isonclust": _backend_isonclust,
}


def split_peaks(
    method: str,
    *,
    comps: pd.DataFrame,
    reads_df: pd.DataFrame,
    peak_subsets: Dict[str, List[str]],
    wt_amplicon_length: float,
    threads: int = 1,
    work_dir: Optional[str] = None,
    cluster_wt_peak: bool = False,
    min_child_fraction: float = 0.15,
    min_child_reads: int = 20,
    gmm_kwargs: Optional[dict] = None,
    tool_kwargs: Optional[dict] = None,
) -> PeakRefineResult:
    """Split GMM peaks into haplotypes using the named method."""
    if method not in BACKENDS:
        raise ValueError(
            f"Unknown haplotype method {method!r}; available: {sorted(BACKENDS)}"
        )
    if comps.empty or reads_df.empty:
        return PeakRefineResult(comps=comps, reads_df=reads_df, peak_subsets=peak_subsets)

    created_tmp = False
    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="haplotype_split_")
        created_tmp = True
    else:
        os.makedirs(work_dir, exist_ok=True)

    try:
        return BACKENDS[method](
            comps=comps,
            reads_df=reads_df,
            peak_subsets=peak_subsets,
            wt_amplicon_length=wt_amplicon_length,
            threads=threads,
            work_dir=work_dir,
            cluster_wt_peak=cluster_wt_peak,
            min_child_fraction=min_child_fraction,
            min_child_reads=min_child_reads,
            gmm_kwargs=gmm_kwargs,
            tool_kwargs=tool_kwargs,
        )
    finally:
        if created_tmp:
            shutil.rmtree(work_dir, ignore_errors=True)
