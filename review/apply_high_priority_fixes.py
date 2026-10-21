from pathlib import Path
p=Path('GMM_peaks.py')
s=p.read_text()
s=s.replace('    wt_amplicon_length=336,\n    force_k=None,','    wt_amplicon_length=336,\n    wt_peak_tolerance=5.0,\n    force_k=None,',1)
a=s.index('    # Give each peak an alias and flag the WT peak')
b=s.index('    #Assign aliases',a)
s=s[:a]+'''    # A sample need not contain WT. Do not relabel an arbitrary ITD as WT.
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

'''+s[b:]
s=s.replace('''        wt_mean = float(
            comps.loc[(comps["mean_bp"] - wt_amplicon_length).abs().idxmin(), "mean_bp"]
        )''','''        wt_mean = float(wt_amplicon_length)''')
p.write_text(s)
p=Path('Helper_functions.py'); s=p.read_text()
s=s.replace('''        if abs(val) >= z_thresh:
''','''        # Z measures separation from the panel mean, not from the runner-up.
        # Require a probability margin too, regardless of reference count.
        margin = row.get("prob_delta", row.get("delta", 0.0))
        if val >= z_thresh and margin >= delta_thresh:
''',1)
# Count both directions from gapless block boundaries, including deletions.
s=s.replace('''            gaps_near = sum(
                ins_len
                for t_end, _q_end, ins_len in insertions_from_aligned_blocks(aln_blocks)
                if lo <= t_end <= hi
            )''','''            t_blocks, q_blocks = aln_blocks
            gaps_near = 0
            for i in range(len(t_blocks) - 1):
                t_end, t_next = int(t_blocks[i][1]), int(t_blocks[i + 1][0])
                q_end, q_next = int(q_blocks[i][1]), int(q_blocks[i + 1][0])
                # A deletion spans a target interval; count it if that interval
                # overlaps the ITD window, even when it starts before the window.
                if t_end <= hi and t_next >= lo:
                    gaps_near += (t_next - t_end) + (q_next - q_end)''')
s=s.replace('''            # Residual insertions in the read relative to the ITD reference.''','''            # Residual insertions AND deletions relative to the ITD reference.''')
p.write_text(s)
p=Path('Pairwise_aligment_toolkit.py'); s=p.read_text()
s=s.replace('df_second[["read_id", "ref_alias", "metric_value"]]','df_second[["read_id", "ref_alias", "metric_value", "prob_score"]]')
s=s.replace('''    df_best["delta"] = df_best["metric_value"] - df_best["metric_value_second"]''','''    df_best["delta"] = df_best["metric_value"] - df_best["metric_value_second"]
    df_best["prob_delta"] = df_best["prob_score"] - df_best["prob_score_second"]''')
p.write_text(s)
p=Path('Nano_ITDseeker.py'); s=p.read_text()
needle='''    parser.add_argument(
        "--per-peak-read-assignment-mode",'''
s=s.replace(needle,'''    parser.add_argument(
        "--wt-peak-tolerance", type=float, default=5.0,
        help="Maximum distance from expected WT length to label a peak WT (bp; default: 5).",
    )
'''+needle,1)
s=s.replace('''    args = parser.parse_args()
''','''    args = parser.parse_args()
    if args.wt_peak_tolerance < 0 or args.wt_peak_tolerance >= args.min_itd_size:
        parser.error("--wt-peak-tolerance must be nonnegative and smaller than --min-itd-size")
''',1)
s=s.replace('''            force_k=force_k,
            wt_amplicon_length=wt_amplicon_length,''','''            force_k=force_k,
            wt_amplicon_length=wt_amplicon_length,
            wt_peak_tolerance=args.wt_peak_tolerance,''',1)
s=s.replace('''            cluster_wt_peak=cluster_wt_peak,
            min_child_fraction''','''            cluster_wt_peak=cluster_wt_peak,
            wt_peak_tolerance=args.wt_peak_tolerance,
            min_child_fraction''',1)
p.write_text(s)
p=Path('.gitignore'); s=p.read_text(); s+='\n# Local review environments, supplied data and generated validation artifacts\n.review-venv/\n.validation-venv/\nbam_data/\nreview/baseline/\nreview/*.tar\nreview/runs/\nreview/synthetic/\n'; p.write_text(s)
