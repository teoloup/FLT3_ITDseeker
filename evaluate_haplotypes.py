#!/usr/bin/env python3
"""
Score a Nano_ITDseeker run against the truth tables from simulate_itd_data.py.

Answers the questions the haplotype work exists to answer:

  - was every true ITD recovered, or were some merged into one call?
  - were spurious haplotypes invented (over-splitting)?
  - how close is the reported allele frequency to the truth?
  - did the consensus stop filling with `N`?

Two matching subtleties, both real rather than cosmetic:

  * A tandem duplication of length L inserted at position P is the same event as
    one inserted at P-L; the aligner picks either. Positions are therefore
    compared modulo the ITD length.
  * For the same reason the extracted insertion sequence can come out cyclically
    rotated, so sequences are compared allowing rotation.

Usage:
    python evaluate_haplotypes.py --run-dir out/ --sample simA \\
        --truth-dir sim_data --scenario A --method isonclust
"""

import argparse
import csv
import glob
import json
import os
import re
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

AMPLICON_START_HG38 = 28033881


# ---------------------------------------------------------------------------
# Sequence comparison
# ---------------------------------------------------------------------------

def is_rotation(a: str, b: str) -> bool:
    return len(a) == len(b) and len(a) > 0 and b in (a + a)


def identity(a: str, b: str) -> float:
    """Ungapped identity over the overlapping prefix, ignoring N."""
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    same = sum(1 for i in range(n) if a[i] == b[i] and a[i] != "N")
    return same / max(len(a), len(b))


def best_rotational_identity(a: str, b: str) -> float:
    """Highest identity over all cyclic rotations of `a`."""
    if not a or not b:
        return 0.0
    best = 0.0
    for i in range(len(a)):
        best = max(best, identity(a[i:] + a[:i], b))
    return best


def same_itd(called_seq: str, called_len: int, called_pos: int,
             true_seq: str, true_len: int, true_pos: int,
             len_tol: int = 2, ident_min: float = 0.85) -> bool:
    if abs(called_len - true_len) > len_tol:
        return False
    # position, allowing the tandem-duplication placement ambiguity
    pos_ok = min(abs(called_pos - true_pos),
                 abs(called_pos - (true_pos - true_len)),
                 abs((called_pos - called_len) - true_pos)) <= len_tol + 2
    if is_rotation(called_seq, true_seq):
        return True
    return pos_ok and best_rotational_identity(called_seq, true_seq) >= ident_min


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_truth(truth_dir: str, scenario: str):
    haps = []
    with open(os.path.join(truth_dir, f"truth_haplotypes_{scenario}.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["itd_len"] or int(r["itd_len"]) == 0:
                continue
            haps.append({
                "name": r["haplotype"],
                "af": float(r["realised_af"]),
                "itd_len": int(r["itd_len"]),
                "genomic_pos": int(r["genomic_ins_pos"]),
                "seq": r["itd_seq"],
            })
    reads = {}
    p = os.path.join(truth_dir, f"truth_reads_{scenario}.tsv")
    if os.path.exists(p):
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                reads[r["read_id"]] = r["haplotype"]
    return haps, reads


def load_calls(run_dir: str, sample: str) -> List[Dict]:
    vcf = os.path.join(run_dir, f"{sample}_FLT3_ITD_calls.vcf")
    calls = []
    if not os.path.exists(vcf):
        return calls
    for line in open(vcf):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        calls.append({
            "alias": f[2],
            "pos": int(f[1]),
            "itd_len": int(info.get("ITD_LEN", 0)),
            "af": float(info.get("AF", 0.0)),
            "af_gmm": float(info.get("AF_GMM", "nan") or "nan"),
            "dp": int(info.get("DP", 0)),
            "seq": f[4][1:],  # ALT minus the anchor base
        })
    return calls


def load_consensus(run_dir: str, sample: str) -> Dict[str, Dict]:
    p = os.path.join(run_dir, "flt3_data", f"{sample}_itd_consensus_seq.tsv")
    out = {}
    if not os.path.exists(p):
        return out
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["peak_alias"]] = {
                "seq": r["consensus_seq"],
                "len": int(r["consensus_len"]),
                "n_reads": int(r["n_total_reads"]),
                "n_unique": int(r["n_unique"]),
                "n_count": r["consensus_seq"].count("N"),
            }
    return out


def peak_read_assignments(run_dir: str, sample: str) -> Dict[str, str]:
    """read_id -> peak alias, from the per-read insertion table."""
    p = os.path.join(run_dir, "flt3_data", f"{sample}_itd_insertions.tsv")
    out = {}
    if not os.path.exists(p):
        return out
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["read_id"]] = r["peak_alias"]
    return out


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def evaluate(run_dir: str, sample: str, truth_dir: str, scenario: str,
             method: str) -> Dict:
    haps, truth_reads = load_truth(truth_dir, scenario)
    calls = load_calls(run_dir, sample)
    cons = load_consensus(run_dir, sample)

    matched: Dict[str, Optional[Dict]] = {h["name"]: None for h in haps}
    used = set()
    for h in haps:
        for i, c in enumerate(calls):
            if i in used:
                continue
            if same_itd(c["seq"], c["itd_len"], c["pos"],
                        h["seq"], h["itd_len"], h["genomic_pos"]):
                matched[h["name"]] = c
                used.add(i)
                break

    extra = [c for i, c in enumerate(calls) if i not in used]

    rows = []
    for h in haps:
        c = matched[h["name"]]
        rows.append({
            "haplotype": h["name"],
            "true_af": round(h["af"], 4),
            "true_len": h["itd_len"],
            "recovered": int(c is not None),
            "called_alias": c["alias"] if c else "",
            "called_len": c["itd_len"] if c else "",
            "called_af": round(c["af"], 4) if c else "",
            "called_af_gmm": round(c["af_gmm"], 4) if c else "",
            "af_error": round(c["af"] - h["af"], 4) if c else "",
            "af_gmm_error": (round(c["af_gmm"] - h["af"], 4)
                             if c and c["af_gmm"] == c["af_gmm"] else ""),
            "len_error": (c["itd_len"] - h["itd_len"]) if c else "",
        })

    n_rec = sum(r["recovered"] for r in rows)
    af_errs = [abs(r["af_error"]) for r in rows if r["af_error"] != ""]
    gmm_errs = [abs(r["af_gmm_error"]) for r in rows if r["af_gmm_error"] != ""]
    total_n = sum(v["n_count"] for v in cons.values())
    peaks_with_n = sum(1 for v in cons.values() if v["n_count"] > 0)

    # read-level agreement, when both truth reads and a peak assignment exist
    cluster_metrics = {}
    assign = peak_read_assignments(run_dir, sample)
    shared = [r for r in assign if r in truth_reads]
    if shared:
        try:
            from sklearn.metrics import (adjusted_rand_score, completeness_score,
                                         homogeneity_score)
            y_true = [truth_reads[r] for r in shared]
            y_pred = [assign[r] for r in shared]
            cluster_metrics = {
                "n_reads_scored": len(shared),
                "ari": round(adjusted_rand_score(y_true, y_pred), 4),
                "homogeneity": round(homogeneity_score(y_true, y_pred), 4),
                "completeness": round(completeness_score(y_true, y_pred), 4),
            }
        except ImportError:
            pass

    return {
        "method": method,
        "sample": sample,
        "scenario": scenario,
        "n_true_itds": len(haps),
        "n_recovered": n_rec,
        "n_missed": len(haps) - n_rec,
        "n_extra_calls": len(extra),
        "extra_aliases": [c["alias"] for c in extra],
        "max_abs_af_error": round(max(af_errs), 4) if af_errs else None,
        "mean_abs_af_error": round(sum(af_errs) / len(af_errs), 4) if af_errs else None,
        "mean_abs_af_gmm_error": (round(sum(gmm_errs) / len(gmm_errs), 4)
                                  if gmm_errs else None),
        "total_consensus_Ns": total_n,
        "peaks_with_Ns": peaks_with_n,
        "n_consensus_peaks": len(cons),
        "cluster_metrics": cluster_metrics,
        "per_haplotype": rows,
        "consensus": cons,
    }


def print_report(res: Dict) -> None:
    print("=" * 96)
    print(f"method={res['method']}  sample={res['sample']}  scenario={res['scenario']}")
    print("=" * 96)
    hdr = (f"{'haplotype':<9} {'true_af':>8} {'len':>4} | {'found':>5} "
           f"{'called_alias':<14} {'called_af':>9} {'af_err':>8} {'af_gmm_err':>10}")
    print(hdr)
    print("-" * 96)
    for r in res["per_haplotype"]:
        mark = "yes" if r["recovered"] else "MISS"
        print(f"{r['haplotype']:<9} {r['true_af']:>8} {r['true_len']:>4} | {mark:>5} "
              f"{str(r['called_alias']):<14} {str(r['called_af']):>9} "
              f"{str(r['af_error']):>8} {str(r['af_gmm_error']):>10}")
    print("-" * 96)
    print(f"  recovered {res['n_recovered']}/{res['n_true_itds']}   "
          f"missed {res['n_missed']}   extra {res['n_extra_calls']} "
          f"{res['extra_aliases'] if res['extra_aliases'] else ''}")
    print(f"  AF error: mean {res['mean_abs_af_error']}  max {res['max_abs_af_error']}"
          f"   (cluster-count AF error: {res['mean_abs_af_gmm_error']})")
    print(f"  consensus Ns: {res['total_consensus_Ns']} across "
          f"{res['peaks_with_Ns']}/{res['n_consensus_peaks']} peaks")
    if res["cluster_metrics"]:
        cm = res["cluster_metrics"]
        print(f"  read clustering: ARI {cm['ari']}  homogeneity {cm['homogeneity']}  "
              f"completeness {cm['completeness']}  (n={cm['n_reads_scored']})")
    print("  per-peak consensus:")
    for alias, c in sorted(res["consensus"].items(), key=lambda kv: -kv[1]["n_reads"]):
        print(f"    {alias:<14} len={c['len']:<4} n={c['n_reads']:<6} "
              f"unique={c['n_unique']:<5} Ns={c['n_count']:<4} {c['seq'][:46]}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--scenario", default="A")
    ap.add_argument("--method", default="unknown")
    ap.add_argument("--json-out", default=None)
    ap.add_argument("--tsv-out", default=None)
    args = ap.parse_args()

    res = evaluate(args.run_dir, args.sample, args.truth_dir, args.scenario, args.method)
    print_report(res)

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(res, fh, indent=2)
    if args.tsv_out:
        with open(args.tsv_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t",
                               fieldnames=list(res["per_haplotype"][0].keys()))
            w.writeheader()
            w.writerows(res["per_haplotype"])


if __name__ == "__main__":
    main()
