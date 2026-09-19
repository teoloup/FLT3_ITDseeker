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
from typing import Callable, Dict, List, Optional, Tuple

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


def _backend_dada2(**kw) -> PeakRefineResult:
    """Denoise each peak's reads into ASVs with DADA2, via an Rscript wrapper.

    DADA2 is R-only, so this shells out to dada2_cluster.R next to this module.
    Its error model is substitution-oriented; the wrapper uses the settings
    DADA2 documents for long indel-prone reads (PacBio CCS), which is the
    closest supported regime to ONT.
    """
    # Bioconductor lags new R releases, so the system Rscript is often too new
    # for dada2. ITDSEEKER_RSCRIPT points at an interpreter that has it (for
    # example a bioconda env, which pins a compatible R alongside the package).
    rscript = os.environ.get("ITDSEEKER_RSCRIPT") or _require_tool("Rscript", "dada2")
    if not os.path.exists(rscript):
        raise RuntimeError(
            f"ITDSEEKER_RSCRIPT points at {rscript}, which does not exist."
        )
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "dada2_cluster.R")
    if not os.path.exists(script):
        raise RuntimeError(f"dada2 backend needs {script}, which is missing.")

    comps, reads_df = kw["comps"], kw["reads_df"]
    peak_subsets = kw["peak_subsets"]
    opts = kw.get("tool_kwargs") or {}
    omega_a = opts.get("dada2_omega_a", 1e-40)
    band_size = int(opts.get("dada2_band_size", 32))
    hp_penalty = opts.get("dada2_homopolymer_gap_penalty", -1)
    work_dir = kw["work_dir"]

    assignments: Dict[str, Dict[str, str]] = {}
    for alias in _target_peaks(comps, kw.get("cluster_wt_peak", False)):
        read_ids = peak_subsets.get(alias, [])
        if len(read_ids) < kw["min_child_reads"] * 2:
            continue

        peak_dir = os.path.join(work_dir, f"dada2_{alias}")
        os.makedirs(peak_dir, exist_ok=True)
        fq = os.path.join(peak_dir, "reads.fastq")
        n = write_peak_fastq(reads_df, read_ids, fq)
        tsv = os.path.join(peak_dir, "clusters.tsv")
        fa = os.path.join(peak_dir, "asvs.fasta")
        logger.info(
            "[haplotype_split] DADA2 on %s (%d reads, OMEGA_A=%g band=%d hp_gap=%s)",
            alias, n, omega_a, band_size, hp_penalty,
        )

        cmd = [rscript, script, fq, tsv, fa, str(omega_a), str(band_size),
               str(hp_penalty), str(kw["min_child_reads"])]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if res.returncode != 0:
            raise RuntimeError(
                f"dada2_cluster.R failed on peak {alias}: {res.stderr.decode()[-2000:]}"
            )
        for line in res.stderr.decode().splitlines():
            if "[dada2_cluster]" in line:
                logger.info("  %s", line.strip())
        if not os.path.exists(tsv):
            raise RuntimeError(f"dada2_cluster.R produced no cluster table for {alias}")

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


def _pad_reads_to_equal_length(reads_df: pd.DataFrame, read_ids: List[str],
                               path: str, length_pct: float = 5.0
                               ) -> Tuple[int, int, int]:
    """Write a FASTQ with every read the same length, for AmpliCI.

    AmpliCI requires equal-length reads with no ambiguous bases. Reads are
    trimmed from the left (which keeps the amplicon start, and therefore the
    insertion, in frame) and anything shorter than the target is dropped rather
    than padded -- padding would invent bases the error model would then try to
    explain.

    The target is a low percentile of the peak's length distribution rather than
    the mode: trimming to the mode discarded roughly half the reads on real
    data, because ONT deletion errors put a long tail below it. `length_pct`
    trades a few trimmed bases for keeping almost every read.

    Returns (written, dropped_short, dropped_ambiguous).
    """
    subset = reads_df.loc[reads_df["read_id"].isin(read_ids)]
    lengths = [len(s) for s in subset["read_seq"] if s]
    if not lengths:
        return 0, 0, 0
    target_len = int(np.percentile(lengths, length_pct))

    written = short = ambig = 0
    with open(path, "w", newline="\n") as fh:
        for row in subset.itertuples(index=False):
            seq = row.read_seq or ""
            if len(seq) < target_len:
                short += 1
                continue
            seq = seq[:target_len]
            if set(seq) - set("ACGT"):
                ambig += 1
                continue
            qual = (getattr(row, "read_qual", "") or "")[:target_len]
            if len(qual) != target_len:
                qual = "I" * target_len
            fh.write(f"@{row.read_id}\n{seq}\n+\n{qual}\n")
            written += 1
    return written, short, ambig


def _backend_amplici(**kw) -> PeakRefineResult:
    """Denoise each peak's reads into haplotypes with AmpliCI.

    AmpliCI models substitution *and* indel errors and estimates them from the
    sample, which is the closest fit of the three tools to what ONT data needs
    and to the case of two near-identical haplotypes. Its constraints are the
    awkward part: equal-length reads, no ambiguous bases, and an indel rate
    parameter whose default (6e-5) is an Illumina figure. The measured ONT indel
    rate on this data is ~0.7%, so the default here is raised to match; leaving
    it at the Illumina value makes AmpliCI read every ONT indel as a distinct
    haplotype.
    """
    binary = _require_tool("run_AmpliCI", "amplici")
    comps, reads_df = kw["comps"], kw["reads_df"]
    peak_subsets = kw["peak_subsets"]
    opts = kw.get("tool_kwargs") or {}
    indel_rate = float(opts.get("amplici_indel_rate", 0.007))
    abundance = float(opts.get("amplici_abundance", 2.0))
    # AmpliCI refuses to assign a read whose log-likelihood under the winning
    # haplotype is below --log_likelihood, default -100. That is an Illumina
    # figure: on ONT reads the per-read log-likelihoods run to -2000 and below,
    # so the default left 85% of reads unassigned (NA) in testing, which would
    # gut the allele-frequency denominator.
    log_likelihood = float(opts.get("amplici_log_likelihood", -100000.0))
    length_pct = float(opts.get("amplici_length_percentile", 5.0))
    work_dir = kw["work_dir"]

    assignments: Dict[str, Dict[str, str]] = {}
    for alias in _target_peaks(comps, kw.get("cluster_wt_peak", False)):
        read_ids = peak_subsets.get(alias, [])
        if len(read_ids) < kw["min_child_reads"] * 2:
            continue

        peak_dir = os.path.join(work_dir, f"amplici_{alias}")
        os.makedirs(peak_dir, exist_ok=True)
        fq = os.path.join(peak_dir, "reads.fastq")
        n, n_short, n_ambig = _pad_reads_to_equal_length(
            reads_df, read_ids, fq, length_pct=length_pct
        )
        if n < kw["min_child_reads"] * 2:
            logger.info(
                "[haplotype_split] AmpliCI skipping %s: only %d of %d reads survived "
                "the equal-length requirement.", alias, n, len(read_ids),
            )
            continue
        logger.info(
            "[haplotype_split] AmpliCI on %s (%d reads kept, %d too short, %d "
            "ambiguous; indel_rate=%g abundance=%g ll_floor=%g)",
            alias, n, n_short, n_ambig, indel_rate, abundance, log_likelihood,
        )

        base = os.path.join(peak_dir, "out")
        cmd = [binary, "--fastq", fq, "--outfile", base,
               "--abundance", str(abundance), "--indel", str(indel_rate),
               "--log_likelihood", str(log_likelihood)]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # AmpliCI returns a non-zero exit code even on a successful run and logs
        # everything, including ordinary progress, to stderr at INFO level. The
        # only reliable success signal is that it wrote its result file.
        out_file = base + ".out" if os.path.exists(base + ".out") else base
        if not os.path.exists(out_file):
            raise RuntimeError(
                f"run_AmpliCI produced no result file for peak {alias} "
                f"(exit {res.returncode}): {res.stderr.decode()[-2000:]}"
            )

        mapping = _parse_amplici_assignments(base, fq)
        if mapping:
            assignments[alias] = mapping
        else:
            logger.warning(
                "[haplotype_split] AmpliCI produced no read assignments for %s", alias
            )

    return assignments_to_result(
        comps=comps, reads_df=reads_df, peak_subsets=peak_subsets,
        assignments=assignments,
        wt_amplicon_length=kw["wt_amplicon_length"],
        min_child_fraction=kw["min_child_fraction"],
        min_child_reads=kw["min_child_reads"],
    )


def _parse_amplici_assignments(base: str, fastq_path: str) -> Dict[str, str]:
    """Read AmpliCI's per-read haplotype assignments out of its .out file.

    The file is a key/value text report; the assignment vector is one integer
    per input read, in input order, under a header naming it. Read ids come back
    from the FASTQ because AmpliCI reports positions, not names.
    """
    out_file = base if os.path.exists(base) else base + ".out"
    if not os.path.exists(out_file):
        return {}

    ids: List[str] = []
    with open(fastq_path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                ids.append(line[1:].strip().split()[0])

    text = open(out_file).read()
    values: List[str] = []
    capture = False
    for line in text.splitlines():
        stripped = line.strip()
        low = stripped.lower()
        if low.startswith("assignments:") or low.startswith("cluster assignments:"):
            capture = True
            rest = stripped.split(":", 1)[1].strip()
            if rest:
                values.extend(rest.split())
            continue
        if capture:
            if not stripped or ":" in stripped:
                break
            values.extend(stripped.split())

    if len(values) < len(ids):
        return {}
    return {
        rid: f"H{val}" for rid, val in zip(ids, values[:len(ids)])
        if val.lstrip("-").isdigit() and int(val) >= 0
    }


BACKENDS: Dict[str, Callable[..., PeakRefineResult]] = {
    "none": _backend_none,
    "gmm2pass": _backend_gmm2pass,
    "isonclust": _backend_isonclust,
    "dada2": _backend_dada2,
    "amplici": _backend_amplici,
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
