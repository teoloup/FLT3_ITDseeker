import os
import re
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging
import seaborn as sns
from collections import Counter
from typing import NamedTuple, Dict, List
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


class GMMFitResult(NamedTuple):
    gmm: GaussianMixture
    comps: pd.DataFrame
    reads_df: pd.DataFrame
    peak_subsets: Dict[str, List[str]]


class PeakRefineResult(NamedTuple):
    comps: pd.DataFrame
    reads_df: pd.DataFrame
    peak_subsets: Dict[str, List[str]]

def fit_gmm_itds(
    reads_dict,
    *,
    min_gmm_fraction,
    max_itds_detected,
    min_ggmm_peak_distance,
    max_peak_sd,
    assign_width_factor,
    assign_mode,        # "manual", "predict_proba", or "hybrid"
    prob_threshold,
    wt_amplicon_length=336,
    wt_peak_tolerance=5.0,
    force_k=None,
    reg=1e-3,
    seed=42
):
    """
    Fit a Gaussian Mixture Model to read lengths and assign reads to mixture components.

    Returns
    -------
    gmm : GaussianMixture
        Fitted model.
    comps : pd.DataFrame
        Components table with model + effective counts.
    df : pd.DataFrame
        Read-level assignments.
    subsets : dict
        {component_index: list of read IDs assigned to that peak}.
    """

    # --- Step 1: Prepare data ---
    data = []
    for rid, entry in reads_dict.items():
        if isinstance(entry, dict):
            seq = entry.get("seq", "")
            strand = entry.get("strand", "+")
            # Phred+33 string; may be absent for sequence-only inputs.
            qual = entry.get("qual", "")
        else:
            # Backward compatibility with older {read_id: seq} format.
            seq = entry
            strand = "+"
            qual = ""
        data.append((rid, seq, qual, strand, len(seq)))
    df = pd.DataFrame(
        data, columns=["read_id", "read_seq", "read_qual", "strand", "read_len"]
    )
    X = df["read_len"].to_numpy(dtype=float).reshape(-1, 1)
    n = len(X)

    # --- Step 2: Fit GMM ---
    if force_k is not None:
        gmm = GaussianMixture(
            n_components=force_k,
            covariance_type="full",
            reg_covar=reg,
            n_init=5,
            random_state=seed
        ).fit(X)
    else:
        best = None
        # One component is the WT peak, so k must reach max_itds_detected + 1
        # for max_itds_detected ITDs to be representable.
        max_components = max_itds_detected + 1
        for k in range(1, max_components + 1):
            g = GaussianMixture(
                n_components=k,
                covariance_type="full",
                reg_covar=reg,
                n_init=5,
                random_state=seed
            ).fit(X)
            bic = g.bic(X)
            if best is None or bic < best[0]:
                best = (bic, g)
        gmm = best[1]

    # --- Step 3: Build component table ---
    means = gmm.means_.ravel()
    covs = gmm.covariances_
    sds = np.sqrt(covs.reshape(-1)) if covs.ndim == 3 else np.sqrt(covs.ravel())
    weights = gmm.weights_.ravel()

    comps_raw = pd.DataFrame({
        "orig_id": np.arange(len(means)),
        "mean_bp": means,
        "sd_bp": sds,
        "fraction": weights,
        "read_count": (weights * n).round().astype(int)
    }).sort_values("mean_bp").reset_index(drop=True)

    # Drop tiny components
    len_before = len(comps_raw)
    comps = comps_raw[comps_raw["fraction"] >= min_gmm_fraction].copy().reset_index(drop=True)
    comps = comps[comps["sd_bp"] < max_peak_sd].copy().reset_index(drop=True)
    if comps.empty:
        raise RuntimeError("No GMM components passed filtering criteria.")
    if len(comps) < len(comps_raw):
        logging.warning(f"Dropped {len_before - len(comps)} GMM components below min fraction "
                        f"or above max SD; {len(comps)} remain.")
    # --- Step 3b: Merge close peaks ---
    merged, skip = [], set()
    merge_map = {}  # map original IDs to merged IDs
    for i, row in comps.iterrows():
        if i in skip:
            continue
        mu, sd = row["mean_bp"], row["sd_bp"]
        group = [i]
        for j in range(i + 1, len(comps)):
            if abs(comps.loc[j, "mean_bp"] - mu) < min_ggmm_peak_distance:
                group.append(j)
                skip.add(j)
        if len(group) == 1:
            merged.append(row.to_dict())
            for g in group:
                merge_map[int(comps.loc[g, "orig_id"])] = len(merged) - 1
        else:
            logger.info(
                "Merging close GMM peaks: %s (means: %s, SDs: %s)",
                group,
                [round(v, 1) for v in comps.loc[group, "mean_bp"]],
                [round(v, 2) for v in comps.loc[group, "sd_bp"]],
            )
            total_frac = comps.loc[group, "fraction"].sum()
            weights_norm = comps.loc[group, "fraction"] / total_frac
            merged_mu = (comps.loc[group, "mean_bp"] * weights_norm).sum()
            # Law of total variance: a mixture's spread is the weighted mean of the
            # child variances PLUS the spread of the child means about the merged
            # mean. Averaging the child SDs ignores the second term and returns a
            # merged peak far narrower than the population it represents, which then
            # causes manual/hybrid assignment (mean +/- factor*SD) to discard most of
            # the reads that motivated the merge in the first place.
            within_var = (comps.loc[group, "sd_bp"] ** 2 * weights_norm).sum()
            between_var = (((comps.loc[group, "mean_bp"] - merged_mu) ** 2) * weights_norm).sum()
            merged_sd = float(np.sqrt(within_var + between_var))
            merged_row = {
                "mean_bp": merged_mu,
                "sd_bp": merged_sd,
                "fraction": total_frac,
                "read_count": int(round(total_frac * n))
            }
            logger.info(
                "Merged peak: mean=%.1f bp, sd=%.2f bp (weighted mean of child SDs "
                "would have been %.2f bp), fraction=%.4f",
                merged_mu, merged_sd,
                float((comps.loc[group, "sd_bp"] * weights_norm).sum()), total_frac,
            )
            merged.append(merged_row)
            for g in group:
                merge_map[int(comps.loc[g, "orig_id"])] = len(merged) - 1
    logger.info(f"Merged {len(comps) - len(merged)} close GMM peaks; {len(merged)} final peaks.")
    comps = pd.DataFrame(merged).sort_values("mean_bp").reset_index(drop=True)

    # --- Step 4: Assign reads ---
    assignments, confidences = [], []
    probs = gmm.predict_proba(X) if assign_mode != "manual" else None

    for i, L in enumerate(df["read_len"]):
        if assign_mode == "manual":
            matches = [
                idx for idx, r in comps.iterrows()
                if abs(L - r["mean_bp"]) <= assign_width_factor * r["sd_bp"]
            ]
            if len(matches) == 1:
                assignments.append(matches[0]); confidences.append(1.0)
            elif len(matches) > 1:
                assignments.append("ambig"); confidences.append(0.0)
            else:
                assignments.append(None); confidences.append(0.0)

        elif assign_mode == "predict_proba":
            p_vec = probs[i]
            best_orig = int(np.argmax(p_vec))
            conf = float(np.max(p_vec))
            if conf >= prob_threshold:
                best_merged = merge_map.get(best_orig, None)
                assignments.append(best_merged)
                confidences.append(conf)
            else:
                assignments.append("ambig"); confidences.append(conf)

        elif assign_mode == "hybrid":
            p_vec = probs[i]
            best_orig = int(np.argmax(p_vec))
            conf = float(np.max(p_vec))
            best_merged = merge_map.get(best_orig, None)
            if best_merged is None:
                assignments.append("ambig"); confidences.append(conf)
                continue
            mu, sd = comps.loc[best_merged, ["mean_bp", "sd_bp"]]
            in_window = abs(L - mu) <= assign_width_factor * sd
            if conf >= prob_threshold and in_window:
                assignments.append(best_merged)
                confidences.append(conf)
            else:
                assignments.append("ambig"); confidences.append(conf)
        else:
            raise ValueError(f"Unknown assign_mode: {assign_mode}")

    df["gmm_peak"] = assignments
    df["confidence"] = confidences
    df["is_ambiguous"] = df["gmm_peak"].eq("ambig")

    # --- Step 5: Filter and compute effective per-peak support ---
    # Keep only confidently assigned, non-ambiguous reads
    df = df.loc[df["gmm_peak"].notnull() & ~df["is_ambiguous"]].reset_index(drop=True)

    # Count how many reads belong to each valid peak index
    eff_counts = Counter(df["gmm_peak"])
    total_eff = sum(eff_counts.values()) or 1

    logger.info(f"Kept {len(df)} reads after filtering unassigned/ambiguous "
            f"({len(df)/n:.1%} of total)")


    comps["effective_read_count"] = [eff_counts[i] for i in comps.index]
    comps["effective_allele_freq"] = comps["effective_read_count"] / total_eff

    # A sample need not contain WT. Do not relabel an arbitrary ITD as WT.
    closest = (comps["mean_bp"] - wt_amplicon_length).abs().idxmin()
    wt_peak_id = (closest if abs(comps.loc[closest, "mean_bp"] - wt_amplicon_length)
                  <= wt_peak_tolerance else None)
    wt_mean = (float(comps.loc[wt_peak_id, "mean_bp"]) if wt_peak_id is not None
               else float(wt_amplicon_length))
    comps["putative_itd_size"] = comps["mean_bp"] - wt_mean
    comps["is_wt"] = comps.index == wt_peak_id
    if wt_peak_id is None:
        logger.info("No WT peak within %.1f bp of %d; using reference length as baseline.",
                    wt_peak_tolerance, wt_amplicon_length)
    else:
        logger.info("WT peak: mean=%.1f bp (expected=%d bp)", wt_mean, wt_amplicon_length)

    #Assign aliases before sorting(store in a new column, not the DataFrame index)
    comps["peak_alias"] = [
        "WT" if i == wt_peak_id else f"ITD_{i}"
        for i in comps.index
    ]

    # Build alias_map on the original indices (which match df["gmm_peak"])
    alias_map = {i: alias for i, alias in zip(comps.index, comps["peak_alias"])}

    # safely sort comps for reporting/plotting
    comps = comps.sort_values("fraction", ascending=False).reset_index(drop=True)

    # Apply alias_map to per-read assignments
    df["gmm_peak_alias"] = df["gmm_peak"].map(alias_map)

    # Log summary info
    logger.info("GMM fitting complete.")
    logger.info("Method used to assign reads: %s", assign_mode)
    logger.info("Effective counts: %s", eff_counts)
    logger.info("Sum of per-peak effective counts: %d", sum(comps["effective_read_count"]))
    logger.info("Sum of effective AF: %f", comps["effective_allele_freq"].sum())

    # --- Step 6: Build subsets of reads ---
    subsets = {
        i: df.loc[df["gmm_peak_alias"] == i, "read_id"].tolist()
        for i in comps["peak_alias"].unique()
    }
    logger.debug(f"df columns: {df.columns.tolist()}")
    return GMMFitResult(gmm=gmm, comps=comps, reads_df=df, peak_subsets=subsets)

def refine_peak_substructure_once(
    comps,
    reads_df,
    peak_subsets,
    *,
    min_reads_for_refinement=150,
    min_child_fraction=0.15,
    min_subpeak_distance=3.0,
    max_subpeak_sd=5.0,
    min_bic_gain_for_split=10.0,
    wt_amplicon_length=336,
    reg=1e-3,
    seed=42,
):
    """
    One-level local refinement:
    for each non-WT peak, test k=1 vs k=2 on that peak's reads and split once
    if evidence supports two close but distinct subpeaks.

    Returns
    -------
    comps_refined : pd.DataFrame
        Updated components table (WT + unsplit peaks + split child peaks).
    reads_df_refined : pd.DataFrame
        reads_df with updated gmm_peak_alias for split peaks.
    peak_subsets_refined : dict
        Updated mapping {peak_alias: [read_ids]}.
    """
    if comps.empty or reads_df.empty:
        return PeakRefineResult(comps=comps, reads_df=reads_df, peak_subsets=peak_subsets)

    reads_df_refined = reads_df.copy()
    peak_subsets_refined = {k: list(v) for k, v in peak_subsets.items()}

    wt_rows = comps.loc[comps["peak_alias"].str.upper() == "WT"]
    if wt_rows.empty:
        wt_mean = float(wt_amplicon_length)
    else:
        wt_mean = float(wt_rows.iloc[0]["mean_bp"])

    split_models = {}
    split_alias_map = {}

    for _, row in comps.iterrows():
        parent_alias = str(row["peak_alias"])
        if parent_alias.upper() == "WT":
            continue

        read_ids = peak_subsets_refined.get(parent_alias, [])
        n_parent = len(read_ids)
        if n_parent < min_reads_for_refinement:
            logger.debug(
                f"[refine_peak_substructure_once] Skip {parent_alias}: n={n_parent} < min_reads_for_refinement={min_reads_for_refinement}"
            )
            continue

        X = reads_df_refined.loc[
            reads_df_refined["read_id"].isin(read_ids), "read_len"
        ].to_numpy(dtype=float).reshape(-1, 1)
        if X.shape[0] < min_reads_for_refinement:
            continue

        g1 = GaussianMixture(
            n_components=1,
            covariance_type="full",
            reg_covar=reg,
            n_init=5,
            random_state=seed,
        ).fit(X)
        g2 = GaussianMixture(
            n_components=2,
            covariance_type="full",
            reg_covar=reg,
            n_init=5,
            random_state=seed,
        ).fit(X)

        bic_gain = g1.bic(X) - g2.bic(X)
        means = g2.means_.ravel()
        covs = g2.covariances_
        sds = np.sqrt(covs.reshape(-1)) if covs.ndim == 3 else np.sqrt(covs.ravel())
        weights = g2.weights_.ravel()

        order = np.argsort(means)
        means = means[order]
        sds = sds[order]
        weights = weights[order]
        mean_delta = abs(means[1] - means[0])

        labels_raw = g2.predict(X)
        labels = np.array([0 if x == order[0] else 1 for x in labels_raw], dtype=int)
        child_counts = np.array([(labels == i).sum() for i in [0, 1]], dtype=int)
        child_fracs = child_counts / max(int(X.shape[0]), 1)

        child1_ok = child_fracs[0] >= min_child_fraction
        child2_ok = child_fracs[1] >= min_child_fraction
        sd_ok = bool(np.all(sds <= max_subpeak_sd))
        dist_ok = mean_delta >= min_subpeak_distance
        bic_ok = bic_gain >= min_bic_gain_for_split

        if not (bic_ok and dist_ok and sd_ok and child1_ok and child2_ok):
            logger.debug(
                f"[refine_peak_substructure_once] No split for {parent_alias}: "
                f"bic_gain={bic_gain:.2f} (>= {min_bic_gain_for_split}), "
                f"delta={mean_delta:.2f} (>= {min_subpeak_distance}), "
                f"sds=({sds[0]:.2f},{sds[1]:.2f}) (<= {max_subpeak_sd}), "
                f"child_fracs=({child_fracs[0]:.3f},{child_fracs[1]:.3f}) (>= {min_child_fraction})"
            )
            continue

        child_aliases = [f"{parent_alias}_S1", f"{parent_alias}_S2"]
        split_alias_map[parent_alias] = child_aliases
        split_models[parent_alias] = {
            "means": means,
            "sds": sds,
            "labels": labels,
            "child_counts": child_counts,
        }

        local_ids = reads_df_refined.loc[
            reads_df_refined["read_id"].isin(read_ids), "read_id"
        ].to_numpy()
        child_ids_0 = local_ids[labels == 0].tolist()
        child_ids_1 = local_ids[labels == 1].tolist()

        peak_subsets_refined.pop(parent_alias, None)
        peak_subsets_refined[child_aliases[0]] = child_ids_0
        peak_subsets_refined[child_aliases[1]] = child_ids_1

        reads_df_refined.loc[
            reads_df_refined["read_id"].isin(child_ids_0), "gmm_peak_alias"
        ] = child_aliases[0]
        reads_df_refined.loc[
            reads_df_refined["read_id"].isin(child_ids_1), "gmm_peak_alias"
        ] = child_aliases[1]

        logger.info(
            f"[refine_peak_substructure_once] Split {parent_alias} -> "
            f"{child_aliases[0]}(n={len(child_ids_0)}, mean={means[0]:.2f}, sd={sds[0]:.2f}) and "
            f"{child_aliases[1]}(n={len(child_ids_1)}, mean={means[1]:.2f}, sd={sds[1]:.2f}); "
            f"bic_gain={bic_gain:.2f}, delta={mean_delta:.2f}"
        )

    if not split_models:
        return PeakRefineResult(comps=comps, reads_df=reads_df_refined, peak_subsets=peak_subsets_refined)

    eff_counts = {
        alias: len(ids)
        for alias, ids in peak_subsets_refined.items()
        if alias in set(reads_df_refined["gmm_peak_alias"])
    }
    total_eff = sum(eff_counts.values()) or 1

    refined_rows = []
    for _, row in comps.iterrows():
        alias = str(row["peak_alias"])
        if alias in split_models:
            child_aliases = split_alias_map[alias]
            model = split_models[alias]
            for i, child_alias in enumerate(child_aliases):
                count = int(eff_counts.get(child_alias, 0))
                frac = count / total_eff
                refined_rows.append({
                    "mean_bp": float(model["means"][i]),
                    "sd_bp": float(model["sds"][i]),
                    "fraction": frac,
                    "read_count": count,
                    "effective_read_count": count,
                    "effective_allele_freq": frac,
                    "putative_itd_size": float(model["means"][i] - wt_mean),
                    "is_wt": False,
                    "peak_alias": child_alias,
                    "parent_peak_alias": alias,
                    "refinement_level": 1,
                    "is_refined_child": True,
                })
        else:
            count = int(eff_counts.get(alias, 0))
            frac = count / total_eff
            refined_rows.append({
                "mean_bp": float(row["mean_bp"]),
                "sd_bp": float(row["sd_bp"]),
                "fraction": frac,
                "read_count": count,
                "effective_read_count": count,
                "effective_allele_freq": frac,
                "putative_itd_size": float(row["mean_bp"] - wt_mean),
                "is_wt": bool(str(alias).upper() == "WT"),
                "peak_alias": alias,
                "parent_peak_alias": alias,
                "refinement_level": int(row.get("refinement_level", 0)),
                "is_refined_child": bool(row.get("is_refined_child", False)),
            })

    comps_refined = pd.DataFrame(refined_rows).sort_values("fraction", ascending=False).reset_index(drop=True)
    logger.info(
        f"[refine_peak_substructure_once] Refinement complete: {len(comps)} -> {len(comps_refined)} peaks."
    )
    return PeakRefineResult(
        comps=comps_refined,
        reads_df=reads_df_refined,
        peak_subsets=peak_subsets_refined,
    )

def _root_alias(alias):
    """The GMM peak a split child came from: ITD_1_H2_H1 -> ITD_1."""
    return re.sub(r"(_H\d+(_\d+)?)+$", "", str(alias))


def _shades(base, k):
    """k shades of one colour, darkest first, for the children of one peak."""
    base = np.array(mcolors.to_rgb(base))
    return [tuple(base + (1 - base) * (0.35 * j / max(k - 1, 1))) for j in range(k)]


def plot_gmm_itds(
    reads_df,
    comps,
    *,
    assign_mode=None,
    bins=None,
    out_prefix="gmm_plot",
    title,
    reported_aliases=None,
    validated=None,
    all_read_lengths=None,
    dpi=150
):
    """
    Plot read lengths with the fitted peaks, and a zoomed panel of the ITD region.

    The top panel shows every read; the lower one zooms in on the ITD peaks,
    which the WT peak otherwise flattens. Split peaks share one colour in
    shades, with a note saying how they were split. Each ITD is labelled with
    its reported AF when `validated` is given, so the plot and the VCF agree.

    Parameters
    ----------
    reads_df : pd.DataFrame
        Reads assigned to a peak, from fit_gmm_itds() or a split; must include
        'read_len'.
    comps : pd.DataFrame
        Peaks with ['mean_bp','sd_bp','fraction','effective_allele_freq',
        'peak_alias']; split children also carry 'is_refined_child' and
        'split_by'.
    assign_mode : str
        'manual', 'predict_proba', 'hybrid' (for title).
    bins : int or None
        Histogram bins; None uses 1 bp bins, since read lengths are integers.
    out_prefix : str
        File prefix for saved figure.
    title : str
        Plot title.
    reported_aliases : set[str] or None
        Aliases that survived validation and the allele-frequency filter. Peaks
        outside this set are drawn dimmed and dotted and labelled "not reported",
        so the plot cannot be read as claiming more ITDs than the VCF contains.
        None draws every peak as reported, which is correct before validation has
        run.
    validated : dict or None
        {alias: {"af": float, "reads": int, "itd_len": int}} for reported ITDs.
    all_read_lengths : array-like or None
        Lengths of every trimmed read, drawn behind the peak-assigned reads so
        the reads outside every peak stay visible.
    dpi : int
        Image resolution.
    """
    mpl_logger = logging.getLogger('matplotlib')
    prev_level = mpl_logger.getEffectiveLevel()
    mpl_logger.setLevel(logging.INFO)


    # --- Prepare data ---
    if reads_df.empty or "read_len" not in reads_df.columns:
        logger.warning("[plot_gmm_itds] No assigned reads available for plotting. Skipping GMM plot.")
        mpl_logger.setLevel(prev_level)
        return

    x = reads_df["read_len"].to_numpy(dtype=float)
    if x.size == 0:
        logger.warning("[plot_gmm_itds] Read-length array is empty. Skipping GMM plot.")
        mpl_logger.setLevel(prev_level)
        return
    x_all = (np.asarray(all_read_lengths, dtype=float)
             if all_read_lengths is not None and len(all_read_lengths) else x)

    n = len(x)
    x_min = float(min(x.min(), x_all.min()))
    x_max = float(max(x.max(), x_all.max()))
    if bins is None:
        edges = np.arange(np.floor(x_min) - 0.5, np.ceil(x_max) + 1.5, 1.0)
    else:
        edges = np.linspace(x_min, x_max if x_max > x_min else x_min + 1.0, int(bins) + 1)
    bin_width = float(edges[1] - edges[0])
    xs = np.linspace(edges[0], edges[-1], 3000)

    if assign_mode == "manual":
        assign_mode_title = "Manual assignment"
    elif assign_mode == "predict_proba":
        assign_mode_title = "Probabilistic assignment"
    elif assign_mode == "hybrid":
        assign_mode_title = "Hybrid assignment"
    else:
        assign_mode_title = None

    rows = sorted(comps.iterrows(), key=lambda kv: float(kv[1]["mean_bp"]))
    is_wt = {r["peak_alias"]: bool(r.get("is_wt", False)) or str(r["peak_alias"]).upper() == "WT"
             for _, r in rows}
    itd_rows = [r for _, r in rows if not is_wt[r["peak_alias"]]]

    # One colour per original peak; split children get shades of it.
    roots = sorted({_root_alias(r["peak_alias"]) for r in itd_rows},
                   key=lambda a: min(float(r["mean_bp"]) for r in itd_rows
                                     if _root_alias(r["peak_alias"]) == a))
    base = dict(zip(roots, sns.color_palette("husl", max(len(roots), 1))))
    colour = {}
    for root in roots:
        members = sorted((r for r in itd_rows if _root_alias(r["peak_alias"]) == root),
                         key=lambda r: float(r["mean_bp"]))
        for r, c in zip(members, _shades(base[root], len(members))):
            colour[r["peak_alias"]] = c

    def curve(row):
        sd = max(float(row["sd_bp"]), 1e-6)
        pdf = np.exp(-0.5 * ((xs - float(row["mean_bp"])) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))
        return pdf * n * bin_width * float(row["fraction"])

    def reported(alias):
        return reported_aliases is None or is_wt[alias] or str(alias) in reported_aliases

    def style(alias):
        if is_wt[alias]:
            return dict(color="black", lw=2.5, ls="-", alpha=0.9)
        if not reported(alias):
            return dict(color="#9aa3ad", lw=1.4, ls=":", alpha=0.75)
        return dict(color=colour[alias], lw=2.0, ls="--", alpha=0.95)

    def legend_label(row):
        alias = row["peak_alias"]
        text = f"{alias}: μ={float(row['mean_bp']):.1f}, σ={float(row['sd_bp']):.1f}"
        eff = row.get("effective_allele_freq", np.nan)
        if bool(row.get("is_refined_child", False)):
            text += f", peak reads {eff * 100:.1f}% (split)"
        else:
            text += f", GMM weight {float(row['fraction']) * 100:.1f}%, peak reads {eff * 100:.1f}%"
        if not reported(alias):
            text += "  [not reported]"
        return text

    # A few very long reads would otherwise stretch the axis far past every peak.
    right = max([float(np.percentile(x_all, 99.5))]
                + [float(r["mean_bp"]) + 4 * float(r["sd_bp"]) for _, r in rows]) + 5
    beyond = int((x_all > right).sum())
    off_scale = f"; {beyond:,} longer than {right:.0f} bp off-scale" if beyond else ""

    def histograms(ax):
        if x_all is not x:
            ax.hist(x_all, bins=edges, color="#d9d9d9", edgecolor="none",
                    label=f"All trimmed reads (n={len(x_all):,}{off_scale})")
            ax.hist(x, bins=edges, color="#a6a6a6", edgecolor="none",
                    label=f"Assigned to a peak (n={n:,})")
        else:
            ax.hist(x, bins=edges, color="#bdbdbd", edgecolor="none",
                    label=f"Reads (n={n:,})")

    zoom = None
    if itd_rows:
        lo = min(float(r["mean_bp"]) - 4 * float(r["sd_bp"]) for r in itd_rows)
        hi = max(float(r["mean_bp"]) + 4 * float(r["sd_bp"]) for r in itd_rows)
        mid, half = (lo + hi) / 2, max((hi - lo) / 2, 10.0)
        zoom = (max(mid - half, edges[0]), min(mid + half, edges[-1]))

    if zoom is None:
        fig, axes = plt.subplots(1, 1, figsize=(10, 5), dpi=dpi, squeeze=False)
    else:
        fig, axes = plt.subplots(2, 1, figsize=(10, 9), dpi=dpi, squeeze=False,
                                 gridspec_kw=dict(height_ratios=[1, 1.25]))
    top = axes[0, 0]

    # --- Full view ---
    histograms(top)
    for _, row in rows:
        top.plot(xs, curve(row), label=legend_label(row), **style(row["peak_alias"]))
    if zoom is not None:
        top.axvspan(*zoom, color="#f6dcb4", alpha=0.35, zorder=0)
    top.set_xlim(edges[0], min(right, edges[-1]))
    top.set_xlabel("Read length (bp)")
    top.set_ylabel("Read count")
    top.set_title(f"{title}\nMode: {assign_mode_title or 'GMM model only'}")
    top.legend(loc="upper right", fontsize=7.5, frameon=True)
    top.grid(alpha=0.3)

    # --- ITD region ---
    if zoom is not None:
        ax = axes[1, 0]
        histograms(ax)
        in_zoom = (xs >= zoom[0]) & (xs <= zoom[1])
        peak_height = {}
        for row in itd_rows:
            y = curve(row)
            ax.plot(xs, y, **style(row["peak_alias"]))
            peak_height[row["peak_alias"]] = float(y[in_zoom].max()) if in_zoom.any() else 0.0
        counts_all, _ = np.histogram(x_all, bins=edges)
        centres = (edges[:-1] + edges[1:]) / 2
        visible = (centres >= zoom[0]) & (centres <= zoom[1])
        y_top = 1.6 * max([1.0] + list(counts_all[visible]) + list(peak_height.values()))
        ax.set_xlim(*zoom)
        ax.set_ylim(0, y_top)

        # Labels over each ITD, stacked when two sit within a few bp.
        level, last_x = 0, None
        for row in sorted(itd_rows, key=lambda r: float(r["mean_bp"])):
            alias, mu = row["peak_alias"], float(row["mean_bp"])
            level = level + 1 if last_x is not None and mu - last_x < 5 else 0
            last_x = mu
            info = (validated or {}).get(str(alias))
            if info:
                text = f"{alias} · {info['itd_len']} bp\nAF {info['af'] * 100:.2f}% ({info['reads']:,} reads)"
            elif not reported(alias):
                text = f"{alias}\nnot reported"
            else:
                text = f"{alias}\n{float(row.get('effective_allele_freq', np.nan)) * 100:.1f}% of peak reads"
            # above the taller of the curve and the bars around it
            near = visible & (np.abs(centres - mu) <= 4)
            height = max([peak_height[alias]] + list(counts_all[near]))
            y = min(height / y_top + 0.03 + 0.16 * level, 0.9)
            ink = tuple(0.75 * np.array(mcolors.to_rgb(style(alias)["color"])))
            ax.text(mu, y, text, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=8, color=ink)

        notes = []
        for root in roots:
            members = sorted((r for r in itd_rows if _root_alias(r["peak_alias"]) == root),
                             key=lambda r: float(r["mean_bp"]))
            if len(members) > 1:
                how = sorted({str(r.get("split_by", "") or "") for r in members} - {""})
                notes.append(f"{root} split into {', '.join(r['peak_alias'] for r in members)}"
                             + (f" by {' and '.join(how)}" if how else ""))
        if notes:
            ax.text(0.01, 0.98, "\n".join(notes), transform=ax.transAxes, ha="left", va="top",
                    fontsize=8, color="#333333")
        ax.set_xlabel("Read length (bp)")
        ax.set_ylabel("Read count")
        ax.set_title("ITD region (shaded above)", fontsize=10)
        ax.grid(alpha=0.3)

    fig.tight_layout()

    # --- Save ---
    fname = f"{out_prefix}_plot.png"
    fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"[plot_gmm_itds] Saved: {os.path.abspath(fname)}")

    mpl_logger.setLevel(prev_level)
