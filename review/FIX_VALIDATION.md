# First correctness fixes and validation

Branch: `fix/asv-validation-correctness`, based on ASV commit `26cdfc1`.
No merge to main and no remote push. Inputs were limited to the five BAMs in
`bam_data/test_bam` and generated synthetic BAMs. `bam_data/all_bam` was not used.

## Changes

- Require a best-versus-runner-up probability margin for multi-reference support;
  a large z-score alone no longer accepts a tie.
- Allow the WT peak to be absent. Added `--wt-peak-tolerance` (5 bp by default),
  and use configured WT length as the missing-WT baseline in both split paths.
- Preserve an eligible WT child after WT clustering and give the remaining
  children usable, collision-resistant `ITD_WT_H*` aliases.
- Count breakpoint-window deletions as well as insertions in validation.
- Pass `--threads` to DADA2 error learning and denoising; fixed runs logged 4 threads.
- Require Cutadapt >=5.2: the installed 4.8 rejected rightmost matching on the
  linked 3-prime adapter. This requirement follows the
  [Cutadapt release notes](https://cutadapt.readthedocs.io/en/v5.2/changes.html).
- Avoid a singular-covariance KDE crash when an insertion-length group is constant.

## Supplied test BAMs

Both `gmm2pass` and `dada2` were run before and after the fixes. All 20 real-test
runs completed successfully. The two methods produced the same reported calls
in this test set. All runs used 4 requested threads and an explicit **1% minimum
reporting AF**, to include the known variants below the default 5% cutoff.
This does not establish a 1% detection limit.

| Sample | Length (bp) | Baseline AF (%) | Fixed AF (%) |
|---|---:|---:|---:|
| 10808 | 45 | 41.235 | 41.209 |
| 11531 | 189 | 8.028 | 7.960 |
| 13697 | 24 | 8.687 | 8.589 |
| 13697 | 30 | 3.937 | 3.902 |
| 13697 | 72 | 8.048 | 8.007 |
| 14219 | 30 | 1.489 | 1.470 |
| 14219 | 81 | 37.012 | 36.768 |
| 14417 | No calls | — | — |

No baseline calls were lost. The largest AF change was 0.244 percentage points.
The negative control remained negative. The 81 bp call already existed in the
baseline; metadata lists 80 bp by fragment analysis. That discrepancy is unresolved.
These are before/after comparisons, not a claim that the pipeline AF equals
fragment-analysis AF. Sample 11531 has no fragment-analysis size/AF in the metadata.

## Synthetic results: important limitations

Two noisy scenarios used 2,000 reads each, seed 42, and the repository simulator's
existing error model. Each ran through both methods before and after the fixes
(8 successful executions). Successful execution does not mean correct recovery.

- **No WT, one true 60 bp ITD:** baseline reported no ITD. Fixed code proceeds
  to discovery, but generates a **56 bp consensus**, with only 136 validated reads
  and AF=100% conditional on retained support. This is **not a successful exact
  recovery**. The missing-WT logic is corrected; noisy insertion alignment and
  consensus remain a separate high-priority problem.
- **Scenario A, four true ITDs (three 45 bp haplotypes plus one 30 bp):** both
  versions and both methods reported only the 30 bp call. DADA2 inferred only one
  ASV in the mixed 45 bp peak (704 reads, 702 unique sequences), and its downstream
  consensus was 24 bp with substantial ambiguity. The reported 30 bp AF was
  14.160% after fixes versus a realised truth of 8.000%, reflecting the depleted
  validation denominator. Same-length haplotype separation is **not validated**.
- A separate error-free 60 bp ITD-only control (400 reads, seed 42, one forced GMM
  component) isolates missing-WT discovery from noisy consensus. Baseline reported
  no ITD; fixed code recovered the 60 bp ITD with 400/400 supporting reads. The
  reconstructed VCF allele exactly matches the generated truth despite equivalent
  breakpoint placement. Results are in `runs/clean_control_summary.json`.

These synthetic failures preclude recommending a merge of the ASV branch to main.
Next work should address insertion normalization/consensus integrity and DADA2
separation on known mixed haplotypes, with read-retention and AF checks. Do not tune
parameters merely to make these particular samples pass.

## Verification and reproducibility

The 11-test `unittest` suite covers true WT and missing WT, both size-baseline
fallbacks, WT child aliases, probability-margin handling, a five-reference alignment
with two tied winners, and actual 60-versus-48 bp read alignment with a positive
control. Run from the repository root:

```text
.review-venv/Scripts/python.exe -m unittest discover -s tests -v
```

Full pipeline runs used WSL Ubuntu, the local `.validation-venv` Python environment,
Cutadapt 5.2, and DADA2 1.38.0 in `/tmp/flt3-itd-review-26cdfc1/dada2`.
The R environment is temporary and may need recreating after a system restart.
The original source snapshot is `review/baseline`, from `git archive 26cdfc1`.
`review/run_test_bams.py` generates noisy scenarios with `prepare`, then runs a
selected method (`gmm2pass` or `dada2`). `review/run_clean_control.py` runs the clean
control. The shell launchers contain this workstation's runtime paths.

Machine-readable comparisons and detailed logs/VCFs are under `review/runs/`:
`gmm2pass_summary.json`, `dada2_summary.json`, and `clean_control_summary.json`.
Raw BAMs, metadata, environments, and generated run artifacts are excluded from Git.

The earlier review remains a historical record. Its medium-priority merged-posterior,
incomplete-cluster-membership, final-haplotype-cap, and low-depth-GMM findings remain
open. The new WT tolerance is a configurable heuristic; neither it nor the softmax
score is a calibrated biological confidence measure.
