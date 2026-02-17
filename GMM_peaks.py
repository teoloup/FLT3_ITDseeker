import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import logging
import seaborn as sns
from collections import Counter
from venv import logger
from sklearn.mixture import GaussianMixture
from collections import Counter

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
    data = [(rid, seq, len(seq)) for rid, seq in reads_dict.items()]
    df = pd.DataFrame(data, columns=["read_id", "read_seq", "read_len"])
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
        for k in range(1, max_itds_detected + 1):
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
            logger.info(f"Merging close GMM peaks: {group} (means: {comps.loc[group, 'mean_bp'].tolist()})")
            total_frac = comps.loc[group, "fraction"].sum()
            weights_norm = comps.loc[group, "fraction"] / total_frac
            merged_mu = (comps.loc[group, "mean_bp"] * weights_norm).sum()
            merged_sd = (comps.loc[group, "sd_bp"] * weights_norm).sum()
            merged_row = {
                "mean_bp": merged_mu,
                "sd_bp": merged_sd,
                "fraction": total_frac,
                "read_count": int(round(total_frac * n))
            }
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

    #give alias to peaks and sort by fraction, identify wt peak as the one closest to 336bp, calculate putative ITD size and store in the dataframe
    # Identify WT peak (closest to 336 bp)
    wt_peak_id = (comps["mean_bp"] - 336).abs().idxmin()

    #Compute putative ITD size (bp difference relative to WT)
    comps["putative_itd_size"] = comps["mean_bp"] - comps.loc[wt_peak_id, "mean_bp"]
    comps["is_wt"] = comps.index == wt_peak_id

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
    return gmm, comps, df, subsets

def plot_gmm_itds(
    reads_df,
    comps,
    *,
    assign_mode=None,
    bins=100,
    out_prefix="gmm_plot",
    title,
    dpi=150
):
    """
    Plot histogram of read lengths (counts) + GMM component curves (model vs effective AF).

    Parameters
    ----------
    reads_df : pd.DataFrame
        Full dataframe from fit_gmm_itds(), must include 'read_len'.
    comps : pd.DataFrame
        GMM component summary with ['mean_bp','sd_bp','fraction','effective_allele_freq','peak_alias'].
    assign_mode : str
        'manual', 'predict_proba', 'hybrid' (for title).
    bins : int
        Histogram bins.
    out_prefix : str
        File prefix for saved figure.
    title : str
        Plot title.
    dpi : int
        Image resolution.
    """
    mpl_logger = logging.getLogger('matplotlib')
    prev_level = mpl_logger.getEffectiveLevel()
    mpl_logger.setLevel(logging.INFO)


    # --- Prepare data ---
    x = reads_df["read_len"].to_numpy(dtype=float)
    n = len(x)
    bin_width = (x.max() - x.min()) / bins
    xs = np.linspace(x.min(), x.max(), 2000)

    fig, ax = plt.subplots(figsize=(9, 5), dpi=dpi)

    if assign_mode == "manual":
        assign_mode_title = "Manual assignment"
    elif assign_mode == "predict_proba":
        assign_mode_title = "Probabilistic assignment"
    elif assign_mode == "hybrid":
        assign_mode_title = "Hybrid assignment"

    # --- Histogram of read counts ---
    ax.hist(
        x,
        bins=bins,
        color="lightgray",
        alpha=0.6,
        edgecolor="none",
        label=f"Reads (n={n})"
    )

    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Read count")
    ax.set_title(f"{title}\nMode: {assign_mode_title or 'GMM model only'}")

    # --- Plot Gaussian components (from GMM means/sds) ---
    colors = sns.color_palette("husl", len(comps))
    mix_curve = np.zeros_like(xs)

    for i, (idx, row) in enumerate(comps.iterrows()):
        mu, sd, frac = row["mean_bp"], row["sd_bp"], row["fraction"]
        eff_frac = row.get("effective_allele_freq", np.nan)
        alias = row.get("peak_alias", f"peak_{i+1}")
        is_wt = row.get("is_wt", False)

        # Gaussian PDF scaled to read counts
        pdf = (1/(sd*np.sqrt(2*np.pi))) * np.exp(-0.5 * ((xs - mu)/sd)**2)
        y = pdf * (n * bin_width) * frac
        mix_curve += y

        # Distinguish WT visually
        color = "black" if is_wt else colors[i]
        lw = 2.5 if is_wt else 2.0
        ls = "-" if is_wt else "--"

        ax.plot(xs, y, color=color, lw=lw, ls=ls, alpha=0.9,
                label=f"{alias}: μ={mu:.1f}, σ={sd:.1f}, "
                      f"model={frac*100:.1f}%, eff={eff_frac*100:.1f}%")

        # Annotate with both AFs
        y_max = y.max()
        ax.text(
            mu, y_max * 1.05,
            f"Model AF {frac*100:.1f}%\nFitted AF {eff_frac*100:.1f}%",
            ha="center", va="bottom",
            fontsize=8,
            color=color,
            fontweight="bold" if is_wt else "normal",
            bbox=dict(boxstyle="round,pad=0.3",
                      facecolor="white", edgecolor=color, alpha=0.8)
        )

    # --- Legend + layout ---
    ax.legend(loc="upper right", fontsize=8, frameon=True)
    ax.grid(alpha=0.3)
    fig.tight_layout()

    # --- Save ---
    fname = f"{out_prefix}_plot.png"
    fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"[plot_gmm_itds] Saved: {os.path.abspath(fname)}")

    mpl_logger.setLevel(prev_level)