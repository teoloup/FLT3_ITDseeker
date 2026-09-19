#!/usr/bin/env python3
"""
In-silico FLT3-ITD read generator.

Produces a coordinate-sorted, indexed BAM plus per-read and per-haplotype truth
tables, so haplotype-separation methods can be scored against a known answer.
The real validation set only has fragment-analysis *sizes*, which cannot tell
you which read came from which ITD.

The BAM is written directly with pysam: no aligner and no hg38 index are needed.
Reads are emitted the way the real BAMs carry them, so the whole pipeline runs
unchanged on the output:

  - header declares chr13 at its hg38 length, so amplicon coordinates resolve
  - roughly half the reads carry the 0x10 flag with SEQ reverse-complemented,
    which `samtools fastq` restores and `cutadapt --rc` then re-tags
  - primers are present on both ends, because the pipeline trims them

The error model is matched to the measured profile of sample 14417, the
ITD-negative control in test_data/: mismatch 0.57%, insertion 0.30%, deletion
0.44%, giving ~99.4% median identity. Indels are biased into homopolymers, and
quality strings are generated consistently with the errors actually injected so
the quality-aware clusterers receive meaningful input.

Usage:
    python simulate_itd_data.py --scenario A --out-dir sim_data
"""

import argparse
import csv
import logging
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pysam

logger = logging.getLogger(__name__)

# --- Reference amplicon, identical to the one Nano_ITDseeker.py uses ---------
DEFAULT_REF_WT = (
    "CTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCG"
    "TTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTT"
    "TAAAGATTTTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAAT"
    "TTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCT"
    "GAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTAAATTGCT"
)

CHROM = "chr13"
CHR13_LENGTH_HG38 = 114364328
AMPLICON_START_HG38 = 28033881  # 1-based, matches Nano_ITDseeker.amplicon_coords

# Measured from test_data/14417_2runs_hg38_RG.bam (ITD-negative control), as
# fractions of aligned bases. Deliberately not rounded up: the point is to
# reproduce the real discrimination problem, not an easier or harder one.
ERR_MISMATCH = 0.00569
ERR_INSERTION = 0.00298
ERR_DELETION = 0.00435

# Indels concentrate in homopolymers on ONT. A base inside a run of >= this many
# identical bases gets its indel probability multiplied by HOMOPOLYMER_FACTOR,
# with the non-homopolymer rate scaled down so the genome-wide rate is preserved.
HOMOPOLYMER_MIN_RUN = 3
HOMOPOLYMER_FACTOR = 4.0

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


@dataclass
class Haplotype:
    """One simulated allele.

    `ins_pos` is a 0-based offset into DEFAULT_REF_WT: the insertion is placed
    between ins_pos-1 and ins_pos, matching how build_itd_reference_per_peak
    constructs its references (ref[:pos] + itd + ref[pos:]).

    `dup_len` bases immediately upstream of ins_pos are duplicated, which is what
    an ITD actually is. `divergent` optionally mutates positions within the
    duplicated copy, to create two ITDs identical in both length and position.
    """

    name: str
    af: float
    ins_pos: Optional[int] = None
    dup_len: int = 0
    divergent: Tuple[int, ...] = field(default_factory=tuple)

    def insert_seq(self, ref: str) -> str:
        if self.ins_pos is None or self.dup_len == 0:
            return ""
        start = max(0, self.ins_pos - self.dup_len)
        dup = list(ref[start:self.ins_pos])
        for off in self.divergent:
            if 0 <= off < len(dup):
                # deterministic substitution: next base in the cycle
                dup[off] = BASES[(BASES.index(dup[off]) + 1) % 4]
        return "".join(dup)

    def allele_seq(self, ref: str) -> str:
        if self.ins_pos is None:
            return ref
        return ref[:self.ins_pos] + self.insert_seq(ref) + ref[self.ins_pos:]

    @property
    def itd_len(self) -> int:
        return self.dup_len if self.ins_pos is not None else 0


# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------
# Positions are chosen inside the region the real ITDs occupy (the amplicon runs
# 0..335; real calls land around local offset 60-270).

SCENARIOS: Dict[str, List[Haplotype]] = {
    # A -- discrimination test. ITD_A and ITD_B are the same length at different
    # positions; ITD_D matches ITD_A in both length and position and differs only
    # in sequence. A length-based splitter cannot separate A from B or D.
    "A": [
        Haplotype("WT", af=0.55),
        Haplotype("ITD_A", af=0.20, ins_pos=180, dup_len=45),
        Haplotype("ITD_B", af=0.12, ins_pos=120, dup_len=45),
        Haplotype("ITD_C", af=0.08, ins_pos=180, dup_len=30),
        Haplotype("ITD_D", af=0.05, ins_pos=180, dup_len=45,
                  divergent=(5, 20, 35)),
    ],
    # B -- sensitivity. A same-length pair down at the AF where the real
    # validation set still has fragment-analysis-confirmed ITDs (0.012).
    "B": [
        Haplotype("WT", af=0.965),
        Haplotype("ITD_A", af=0.020, ins_pos=180, dup_len=45),
        Haplotype("ITD_B", af=0.015, ins_pos=120, dup_len=45),
    ],
}


def homopolymer_weights(seq: str) -> np.ndarray:
    """Per-base indel weighting, raised inside homopolymer runs."""
    w = np.ones(len(seq), dtype=float)
    i = 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        run = j - i + 1
        if run >= HOMOPOLYMER_MIN_RUN:
            w[i:j + 1] = HOMOPOLYMER_FACTOR
        i = j + 1
    # preserve the overall rate
    return w / w.mean()


def apply_errors(seq: str, rng: np.random.Generator) -> Tuple[str, str]:
    """Inject substitutions, insertions and deletions; return (seq, qual).

    Quality is assigned from what actually happened to each emitted base rather
    than drawn independently, so a quality-aware denoiser sees a real signal:
    correct bases get high Q, substituted and inserted bases get low Q.
    """
    weights = homopolymer_weights(seq)
    out_bases: List[str] = []
    out_q: List[int] = []

    for i, base in enumerate(seq):
        w = weights[i]
        r = rng.random()

        # deletion: emit nothing
        if r < ERR_DELETION * w:
            continue

        # substitution
        if r < ERR_DELETION * w + ERR_MISMATCH:
            alt = BASES[(BASES.index(base) + 1 + rng.integers(3)) % 4] if base in BASES else base
            out_bases.append(alt)
            out_q.append(int(rng.integers(4, 12)))
        else:
            out_bases.append(base)
            out_q.append(int(rng.integers(22, 41)))

        # insertion after this base
        if rng.random() < ERR_INSERTION * w:
            out_bases.append(BASES[rng.integers(4)])
            out_q.append(int(rng.integers(3, 10)))

    qual = "".join(chr(min(q, 93) + 33) for q in out_q)
    return "".join(out_bases), qual


def simulate(scenario: str, n_reads: int, seed: int, out_dir: str) -> Dict[str, str]:
    haps = SCENARIOS[scenario]
    total_af = sum(h.af for h in haps)
    if abs(total_af - 1.0) > 1e-6:
        raise ValueError(f"Scenario {scenario} AFs sum to {total_af}, not 1.0")

    rng = np.random.default_rng(seed)
    os.makedirs(out_dir, exist_ok=True)

    ref = DEFAULT_REF_WT
    alleles = {h.name: h.allele_seq(ref) for h in haps}

    # draw haplotype per read
    names = [h.name for h in haps]
    probs = np.array([h.af for h in haps], dtype=float)
    assignments = rng.choice(len(haps), size=n_reads, p=probs)

    bam_path = os.path.join(out_dir, f"sim_{scenario}.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": CHR13_LENGTH_HG38}],
        "RG": [{"ID": "1", "SM": f"SIM_{scenario}", "PL": "Nanopore", "LB": "sim"}],
        "PG": [{"ID": "simulate_itd_data", "PN": "simulate_itd_data.py",
                "CL": f"--scenario {scenario} --reads {n_reads} --seed {seed}"}],
    }

    truth_rows = []
    realised = {h.name: 0 for h in haps}

    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for i, hap_idx in enumerate(assignments):
            hap = haps[hap_idx]
            realised[hap.name] += 1

            seq, qual = apply_errors(alleles[hap.name], rng)
            # reads are emitted with primers attached; the pipeline trims them
            is_reverse = bool(rng.random() < 0.5)
            stored = revcomp(seq) if is_reverse else seq
            stored_qual = qual[::-1] if is_reverse else qual

            read_id = f"sim{scenario}_{i:06d}"
            a = pysam.AlignedSegment()
            a.query_name = read_id
            a.query_sequence = stored
            a.flag = 16 if is_reverse else 0
            a.reference_id = 0
            a.reference_start = AMPLICON_START_HG38 - 1  # pysam is 0-based
            a.mapping_quality = 60
            # A plain match CIGAR is enough: the pipeline only uses the BAM to
            # select reads by region and MAPQ, then re-derives everything from
            # the sequence itself.
            a.cigar = [(0, len(stored))]
            a.query_qualities = pysam.qualitystring_to_array(stored_qual)
            a.set_tag("RG", "1")
            bam.write(a)

            truth_rows.append({
                "read_id": read_id,
                "haplotype": hap.name,
                "is_itd": int(hap.ins_pos is not None),
                "itd_len": hap.itd_len,
                "ins_pos_local": hap.ins_pos if hap.ins_pos is not None else "",
                "strand": "-" if is_reverse else "+",
                "observed_len": len(seq),
            })

    pysam.index(bam_path)

    reads_tsv = os.path.join(out_dir, f"truth_reads_{scenario}.tsv")
    with open(reads_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(truth_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(truth_rows)

    haps_tsv = os.path.join(out_dir, f"truth_haplotypes_{scenario}.tsv")
    with open(haps_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["haplotype", "target_af", "realised_af", "n_reads", "itd_len",
                    "ins_pos_local", "genomic_ins_pos", "itd_seq"])
        for h in haps:
            gpos = "" if h.ins_pos is None else AMPLICON_START_HG38 + h.ins_pos - 1
            w.writerow([
                h.name, f"{h.af:.4f}", f"{realised[h.name] / n_reads:.4f}",
                realised[h.name], h.itd_len,
                "" if h.ins_pos is None else h.ins_pos,
                gpos, h.insert_seq(ref),
            ])

    logger.info("Wrote %s (%d reads)", bam_path, n_reads)
    for h in haps:
        logger.info(
            "  %-7s target_af=%.4f realised=%.4f n=%d itd_len=%d",
            h.name, h.af, realised[h.name] / n_reads, realised[h.name], h.itd_len,
        )

    return {"bam": bam_path, "truth_reads": reads_tsv, "truth_haplotypes": haps_tsv}


def self_check(scenario: str, paths: Dict[str, str], n_reads: int) -> List[str]:
    """Confirm the emitted data says what it claims to. Returns failures."""
    problems = []
    haps = SCENARIOS[scenario]
    ref = DEFAULT_REF_WT

    # realised AF should track target within sampling error
    with open(paths["truth_haplotypes"]) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        tgt, got = float(r["target_af"]), float(r["realised_af"])
        tol = max(0.01, 4.0 * (tgt * (1 - tgt) / n_reads) ** 0.5)
        if abs(tgt - got) > tol:
            problems.append(
                f"{r['haplotype']}: realised AF {got:.4f} vs target {tgt:.4f} (tol {tol:.4f})"
            )

    # the same-length pairs that make the scenario interesting must really be equal
    lengths: Dict[int, List[str]] = {}
    for h in haps:
        if h.itd_len:
            lengths.setdefault(h.itd_len, []).append(h.name)
    shared = {k: v for k, v in lengths.items() if len(v) > 1}
    if scenario == "A" and not shared:
        problems.append("scenario A has no two ITDs sharing a length; it cannot test the failure mode")

    # same-length haplotypes must differ in sequence, or there is nothing to split
    for size, names in shared.items():
        seqs = {n: next(h for h in haps if h.name == n).insert_seq(ref) for n in names}
        if len(set(seqs.values())) != len(seqs):
            problems.append(f"ITDs of length {size} have identical sequences: {names}")

    # read counts must match the BAM
    with open(paths["truth_reads"]) as fh:
        n_truth = sum(1 for _ in csv.DictReader(fh, delimiter="\t"))
    n_bam = sum(1 for _ in pysam.AlignmentFile(paths["bam"], "rb").fetch(until_eof=True))
    if n_truth != n_bam or n_bam != n_reads:
        problems.append(f"read counts disagree: truth={n_truth} bam={n_bam} requested={n_reads}")

    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenario", choices=sorted(SCENARIOS), default="A")
    ap.add_argument("--reads", type=int, default=8000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", default="sim_data")
    ap.add_argument("--log-level", default="INFO", choices=["INFO", "DEBUG"])
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s [%(levelname)s] %(message)s")

    paths = simulate(args.scenario, args.reads, args.seed, args.out_dir)

    problems = self_check(args.scenario, paths, args.reads)
    if problems:
        logger.error("Self-check FAILED:")
        for p in problems:
            logger.error("  - %s", p)
        raise SystemExit(1)
    logger.info("Self-check passed.")
    for k, v in paths.items():
        logger.info("  %-18s %s", k, v)


if __name__ == "__main__":
    main()
