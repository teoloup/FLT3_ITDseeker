# Context-aware insertion consensus

Base commit: `7133826`; branch: `fix/asv-validation-correctness`.
This follow-up targets single-haplotype consensus reconstruction, not DADA2 tuning.

## Cause and change

The noisy 60 bp ITD-only challenge yielded insertion boundaries at many positions
(e.g. 595 observations at local boundary 119, 79 at 180). Equivalent duplication
representations carry differently rotated insertion payloads. Previously the code
aligned those payloads without their reference context, discarded low-coverage
MSA columns, and independently attached the resulting 56 bp sequence at the median
read boundary (122). Both length and placement could therefore be wrong.

The new path reconstructs each insertion observation as
`WT[:boundary] + payload + WT[boundary:]`, deduplicates these contextual alleles,
and builds their weighted consensus. Alignment of that consensus to WT recovers
an insertion payload and its boundary together. It does not force the payload to
match WT, pad it to the expected GMM size, or change the MSA voting thresholds.
Consensuses without exactly one insertion of reportable size are skipped with a
warning. Complex multi-insertion haplotypes remain outside this single-insertion
consensus representation.

The consensus table now includes `consensus_ins_pos_ref`,
`raw_median_ins_pos_ref`, and `allele_consensus_len`. The legacy
`median_ins_pos_ref` column aliases the consensus boundary for compatibility.
Validation and reference construction prefer the explicit consensus field.
MSA plots now show the full contextual allele; output insertion sequences remain
insertion payloads. The Python consensus function now requires `ref_seq`.

The requested thread budget is divided across active peak workers, preventing
each member of the per-peak process pool from independently consuming every CPU.
A single peak can use the entire requested budget. Longer contextual alignments
cost more work than isolated insertion alignments.

## Tests

- Eighteen regression tests pass in the Linux pipeline environment, including
  equivalent rotations, coherent sequence/anchor output, a novel non-tandem
  insertion, consensus size bounds, and all earlier validation regressions.
- Synthetic consensus sweep: 400 reads per case, a 60 bp duplication at local
  positions 120, 180, and 240, each at zero, normal, and double the repository
  simulator's error rates. All nine cases recovered the exact full allele,
  not just a 60 bp length. This sweep exercises noisy read alignment, insertion
  extraction, and consensus, without the BAM trimming/GMM/DADA2 stages.
- Full pipeline comparison uses only the five supplied `test_bam` inputs and the
  two existing synthetic challenges, DADA2, four requested threads, and the same
  1% reporting threshold as the earlier report. `all_bam` remains untouched.

Detailed results are recorded in `runs/context_sweep.json`,
`runs/context_unit_tests.log`, and `runs/context_dada2_summary.json`.
The corresponding runners are `sweep_context_consensus.py` and
`run_context_bams.py`; use the same Linux runtime documented in FIX_VALIDATION.md.

This remains an insertion caller: WT flanks are reference context, not independent
consensus evidence about linked SNVs or deletions outside the extracted insertion.
The synthetic error model and limited test set do not establish clinical
sensitivity or an AF detection limit. Same-length haplotype separation and AF
bias from missing candidates remain separate work.

## Full-pipeline results

The complete context-consensus batch kept all calls, sequences, positions and AFs
unchanged relative to `7133826` in the five supplied test BAMs:

| Test BAM | ITD sizes (bp) | AF (%) |
|---|---|---|
| 10808 | 45 | 41.209 |
| 11531 | 189 | 7.960 |
| 13697 | 24, 30, 72 | 8.589, 3.902, 8.007 |
| 14219 | 30, 81 | 1.470, 36.768 |
| 14417 | None | Negative control |

The noisy pure-ITD challenge now recovers the exact 60 bp allele, AF 100%,
validated depth 1536. Previously it produced a wrong 56 bp call with depth 136.

The mixed `sim_A` challenge exposed an unresolved 43 bp N-rich consensus that
initially reached the VCF. A final guard now rejects any contextual consensus
containing non-ACGT bases before creating a candidate reference. This is a
conservative rejection, not successful separation of the three true 45 bp
haplotypes. The final mixed run reports only the 30 bp insertion, AF 14.160%,
depth 1024; its simulated AF is 8%. Missing candidates still bias quantification.
The rejected peak remains visible in the MSA plot and console warning.

`runs/context_dada2_summary.json` records the batch BEFORE this final ambiguity
guard; its mixed-sample 43 bp call is a diagnostic failure, not a final result.
Final reruns and exact command audits are under `runs/final_dada2/`.
`runs/final_unit_tests.log` records 18 passing tests, including the ambiguity
rejection, backend removal, and exact command quoting/routing/append checks.

## DADA2-only cleanup and command audit

Removed isONclust and AmpliCI implementations, CLI options and Docker installs.
Removed the gmm2pass backend/chaining and its unused CLI parameters; the initial
GMM and its options remain. `--haplotype-method` accepts `dada2` and `none`.
The standalone historical GMM refinement helper remains for existing regression
coverage, but is not called by the pipeline. Historical benchmark evidence is
retained. The Dockerfile includes the curl/bzip2/certificate dependencies needed
for its DADA2 environment setup; an image build was not tested (no Docker daemon).

Every external samtools, Cutadapt and DADA2 Rscript launch is recorded in the
output root's `<sample>_commands.jsonl`, and printed at INFO level. JSON preserves
argv, resolved executable, cwd, UTC time, shell=False and stream routing; shell
quotes protect adapter semicolons and paths in the display string. DADA2's actual
Rscript wrapper arguments include every exposed setting and the thread count.
The audit records attempted launches and appends across reruns. It survives
intermediate cleanup; deleted FASTQs must be regenerated before replaying a call.

Final smoke verification (`check_final_validation.py`) passed for both synthetic
BAMs and test BAM 10808. The positive BAM retained the exact 45 bp allele,
position, AF 41.209%, and depth 13531. The audit contained respectively 5, 4,
and 4 external invocations, including the expected DADA2 arguments and four
threads. All three audit files remained after the pipeline deleted its temp
folder. The checked summary is `runs/final_dada2_summary.json`.
