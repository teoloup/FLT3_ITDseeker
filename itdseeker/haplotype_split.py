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

After insertion extraction, split_by_insertion_length makes a second pass on
the extracted insertion lengths, which separate ITDs a few bp apart that both
the length GMM and DADA2 can merge.
"""

from .command_audit import record_command

import logging
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .GMM_peaks import PeakRefineResult

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
    wt_peak_tolerance: float = 5.0,
) -> PeakRefineResult:
    """Turn per-peak read->cluster labels into the standard refinement result.

    `assignments` maps peak_alias -> {read_id: cluster_label}. A peak absent from
    it is left untouched.

    Guardrails mirror the existing length-based refinement: a cluster must hold
    at least `min_child_fraction` of the parent's reads and at least
    `min_child_reads` reads to become a haplotype. If fewer than two clusters
    survive, the peak is left unsplit.

    Reads in sub-threshold clusters, and parent reads the backend left without
    a cluster, fold into the largest survivor rather than being dropped, which
    would shrink the AF denominator and inflate every other call.
    """
    reads_out = reads_df.copy()
    subsets_out = {k: list(v) for k, v in peak_subsets.items()}

    len_by_read = dict(zip(reads_df["read_id"], reads_df["read_len"]))

    wt_rows = comps.loc[comps["peak_alias"].astype(str).str.upper() == "WT"]
    if wt_rows.empty:
        wt_mean = float(wt_amplicon_length)
    else:
        wt_mean = float(wt_rows.iloc[0]["mean_bp"])

    split_children: Dict[str, List[Dict]] = {}
    used_aliases = set(subsets_out)

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

        # fold sub-threshold clusters, and reads given no cluster, into the
        # largest survivor
        kept_labels = {lbl for lbl, _ in keep}
        folded = 0
        biggest = keep[0][0]
        merged: Dict[str, List[str]] = {lbl: list(rids) for lbl, rids in keep}
        for lbl, rids in ordered:
            if lbl not in kept_labels:
                merged[biggest].extend(rids)
                folded += len(rids)
        unassigned = [rid for rid in parent_reads if read_to_cluster.get(rid) is None]
        merged[biggest].extend(unassigned)
        if folded or unassigned:
            logger.info(
                "[haplotype_split] %s: folded %d reads from sub-threshold clusters "
                "and %d reads without a cluster into %s_H1.",
                parent_alias, folded, len(unassigned), parent_alias,
            )

        children = []
        for i, (lbl, rids) in enumerate(
            sorted(merged.items(), key=lambda kv: (-len(kv[1]), kv[0])), start=1
        ):
            lengths = np.array([len_by_read[r] for r in rids if r in len_by_read],
                               dtype=float)
            prefix = "ITD_WT" if parent_alias.upper() == "WT" else parent_alias
            base_alias = f"{prefix}_H{i}"
            child_alias = base_alias
            suffix = 2
            while child_alias in used_aliases:
                child_alias = f"{base_alias}_{suffix}"
                suffix += 1
            used_aliases.add(child_alias)
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

        if parent_alias.upper() == "WT":
            # WT-like length is a provisional label. All other children must use
            # variant aliases understood by extraction and competitive validation.
            closest = min(children, key=lambda ch: abs(ch["mean_bp"] - wt_amplicon_length))
            if abs(closest["mean_bp"] - wt_amplicon_length) <= wt_peak_tolerance:
                closest["alias"] = "WT"
                wt_mean = closest["mean_bp"]
            else:
                wt_mean = float(wt_amplicon_length)

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
                    "is_wt": ch["alias"] == "WT",
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
               str(hp_penalty), str(kw["min_child_reads"]),
               str(max(1, int(kw.get("threads", 1))))]
        record_command(cmd, stdin="inherited", stdout="captured", stderr="captured")
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
        wt_peak_tolerance=kw.get("wt_peak_tolerance", 5.0),
    )


BACKENDS: Dict[str, Callable[..., PeakRefineResult]] = {
    "none": _backend_none,
    "dada2": _backend_dada2,
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
    wt_peak_tolerance: float = 5.0,
    min_child_fraction: float = 0.15,
    min_child_reads: int = 20,
    tool_kwargs: Optional[dict] = None,
) -> PeakRefineResult:
    """Split GMM peaks into haplotypes using the named method."""
    if method not in BACKENDS:
        raise ValueError(
            f"Unknown haplotype method {method!r}; available: {sorted(BACKENDS)} "
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
            wt_peak_tolerance=wt_peak_tolerance,
            min_child_fraction=min_child_fraction,
            min_child_reads=min_child_reads,
            tool_kwargs=tool_kwargs,
        )
    finally:
        if created_tmp:
            shutil.rmtree(work_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Split by insertion length, after insertion extraction
# ---------------------------------------------------------------------------

# Nanopore indels blur one ITD's insertion length by a base or two either way,
# so two lengths count as separate ITDs only this far apart...
MIN_INSERTION_LENGTH_GAP_BP = 3
# ...and only if the read count between them falls to this share of the smaller
# group, so a shoulder of one ITD is not mistaken for a second.
MAX_INSERTION_VALLEY_FRACTION = 0.5


def insertion_length_modes(
    lengths,
    *,
    min_reads: int,
    min_fraction: float,
    min_gap: int = MIN_INSERTION_LENGTH_GAP_BP,
    max_valley: float = MAX_INSERTION_VALLEY_FRACTION,
) -> List[int]:
    """Insertion lengths that each mark a separate ITD within one peak.

    A candidate is a local maximum of the per-length read count holding at
    least `min_reads` reads and `min_fraction` of all of them. Taken largest
    first, a candidate is kept only if it lies `min_gap` bp or more from every
    length already kept and the count between the two drops to `max_valley`
    of the smaller one.
    """
    counts = Counter(int(x) for x in lengths)
    total = sum(counts.values())
    if not total:
        return []
    candidates = sorted(
        (length for length, n in counts.items()
         if n >= min_reads and n / total >= min_fraction
         and n >= counts.get(length - 1, 0) and n > counts.get(length + 1, 0)),
        key=lambda length: (-counts[length], length),
    )
    kept: List[int] = []
    for length in candidates:
        if all(
            abs(length - other) >= min_gap
            and min(counts.get(x, 0) for x in range(min(length, other) + 1, max(length, other)))
            <= max_valley * min(counts[length], counts[other])
            for other in kept
        ):
            kept.append(length)
    return sorted(kept)


def split_by_insertion_length(
    *,
    comps: pd.DataFrame,
    reads_df: pd.DataFrame,
    peak_subsets: Dict[str, List[str]],
    insertions_df: pd.DataFrame,
    wt_amplicon_length: float,
    min_child_fraction: float,
    min_child_reads: int,
    wt_peak_tolerance: float = 5.0,
) -> Tuple[PeakRefineResult, pd.DataFrame]:
    """Split peaks whose reads carry clearly different insertion lengths.

    Two ITDs a few bp apart give read lengths too close for the length GMM to
    separate, and DADA2 can merge them as well: it scores whole reads, where
    the extra bases are one indel among each read's own sequencing errors. The
    extracted insertion lengths are far sharper, since only the inserted
    segment contributes errors, so a second ITD stands out as its own group.

    Every read of a split peak is assigned: reads with an insertion go to the
    nearest length group, the rest to the largest. The usual size guardrails
    apply. Returns the refined peaks and the insertion table relabelled with
    the new aliases.
    """
    assignments: Dict[str, Dict[str, str]] = {}
    for alias, group in insertions_df.groupby("peak_alias", sort=False):
        if str(alias).upper() == "WT" or alias not in peak_subsets:
            continue
        lengths = group["ins_len"].astype(int)
        modes = insertion_length_modes(
            lengths, min_reads=min_child_reads, min_fraction=min_child_fraction
        )
        if len(modes) < 2:
            continue
        counts = Counter(lengths)
        largest = max(modes, key=lambda m: (counts[m], -m))

        def nearest(n):
            # ties go to the group with more reads
            return min(modes, key=lambda m: (abs(m - n), -counts[m]))

        by_read = {rid: f"{nearest(n)}bp" for rid, n in zip(group["read_id"], lengths)}
        assignments[alias] = {
            rid: by_read.get(rid, f"{largest}bp") for rid in peak_subsets[alias]
        }
        logger.info(
            "[insertion_split] %s: separate insertion lengths %s; splitting the peak.",
            alias, ", ".join(f"{m} bp ({counts[m]} reads)" for m in modes),
        )

    if not assignments:
        return PeakRefineResult(comps=comps, reads_df=reads_df,
                                peak_subsets=peak_subsets), insertions_df

    result = assignments_to_result(
        comps=comps,
        reads_df=reads_df,
        peak_subsets=peak_subsets,
        assignments=assignments,
        wt_amplicon_length=wt_amplicon_length,
        min_child_fraction=min_child_fraction,
        min_child_reads=min_child_reads,
        wt_peak_tolerance=wt_peak_tolerance,
    )
    alias_of = dict(zip(result.reads_df["read_id"], result.reads_df["gmm_peak_alias"]))
    relabelled = insertions_df.copy()
    moved = relabelled["peak_alias"].isin(list(assignments))
    relabelled.loc[moved, "peak_alias"] = relabelled.loc[moved, "read_id"].map(alias_of)
    return result, relabelled
