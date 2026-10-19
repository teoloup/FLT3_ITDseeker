from pathlib import Path
p=Path('haplotype_split.py'); s=p.read_text()
s=s.replace('''    min_child_reads: int,
) -> PeakRefineResult:''','''    min_child_reads: int,
    wt_peak_tolerance: float = 5.0,
) -> PeakRefineResult:''',1)
s=s.replace('''        wt_mean = float(
            comps.loc[(comps["mean_bp"] - wt_amplicon_length).abs().idxmin(), "mean_bp"]
        )''','''        wt_mean = float(wt_amplicon_length)''',1)
s=s.replace('''    split_children: Dict[str, List[Dict]] = {}
''','''    split_children: Dict[str, List[Dict]] = {}
    used_aliases = set(subsets_out)
''',1)
s=s.replace('''            child_alias = f"{parent_alias}_H{i}"
''','''            prefix = "ITD_WT" if parent_alias.upper() == "WT" else parent_alias
            base_alias = f"{prefix}_H{i}"
            child_alias = base_alias
            suffix = 2
            while child_alias in used_aliases:
                child_alias = f"{base_alias}_{suffix}"
                suffix += 1
            used_aliases.add(child_alias)
''',1)
s=s.replace('''        split_children[parent_alias] = children
''','''        if parent_alias.upper() == "WT":
            # WT-like length is a provisional label. All other children must use
            # variant aliases understood by extraction and competitive validation.
            closest = min(children, key=lambda ch: abs(ch["mean_bp"] - wt_amplicon_length))
            if abs(closest["mean_bp"] - wt_amplicon_length) <= wt_peak_tolerance:
                closest["alias"] = "WT"
                wt_mean = closest["mean_bp"]
            else:
                wt_mean = float(wt_amplicon_length)

        split_children[parent_alias] = children
''',1)
s=s.replace('''                    "is_wt": False,
                    "peak_alias": ch["alias"],''','''                    "is_wt": ch["alias"] == "WT",
                    "peak_alias": ch["alias"],''',1)
# Pass through each sequence backend conversion.
s=s.replace('''        min_child_reads=kw["min_child_reads"],
    )''','''        min_child_reads=kw["min_child_reads"],
        wt_peak_tolerance=kw.get("wt_peak_tolerance", 5.0),
    )''')
s=s.replace('''    cluster_wt_peak: bool = False,
    min_child_fraction''','''    cluster_wt_peak: bool = False,
    wt_peak_tolerance: float = 5.0,
    min_child_fraction''',1)
s=s.replace('''                cluster_wt_peak=cluster_wt_peak,
                min_child_fraction''','''                cluster_wt_peak=cluster_wt_peak,
                wt_peak_tolerance=wt_peak_tolerance,
                min_child_fraction''',1)
s=s.replace('''            cluster_wt_peak=cluster_wt_peak,
            min_child_fraction''','''            cluster_wt_peak=cluster_wt_peak,
            wt_peak_tolerance=wt_peak_tolerance,
            min_child_fraction''',1)
p.write_text(s)
