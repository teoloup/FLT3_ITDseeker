# Benchmarks

Measurements behind the non-obvious defaults in this pipeline. Recorded here
rather than under `test_output/`, which is gitignored, so the evidence for a
default travels with the code that uses it.

Regenerate any of these with `simulate_itd_data.py` + `evaluate_haplotypes.py`;
the validation BAMs referenced by the real-data tables are not in the repo.

## Haplotype splitting

| file | what it shows |
|---|---|
| `method_comparison.tsv` | gmm2pass / isonclust / dada2 / amplici on simulated scenario A, where three of four ITDs share a length |
| `real_data_method_comparison.tsv` | the same four methods on the real validation samples |
| `real_11531_four_methods.txt`, `real_13697_14219_four_methods.txt` | the raw per-peak output behind that table |
| `isonclust_parameter_sweep.tsv` | why `--isonclust-aligned-threshold` defaults to 0.80 and not isONclust's own 0.4 |
| `min_child_fraction_experiment.tsv` | whether the size guardrail was what blocked the hardest haplotype (it was not, for two of three tools) |
| `truth_haplotypes_A.tsv`, `simA_*_per_haplotype.tsv` | the simulated truth set and per-haplotype scoring |

**Headline.** On simulation the clusterers beat the length-based baseline 3/4
against 1/4. On the real samples the baseline is already correct on all three,
`dada2` matches it, and `isonclust` and `amplici` each lose an ITD. The two
disagree because the simulation was built to defeat length-based separation,
while real co-occurring ITDs in this set differ in size — the case the GMM
already handles. `gmm2pass` therefore remains the default.

## Competitive validation scoring

| file | what it shows |
|---|---|
| `validation_scoring_diagnosis.tsv` | correct-reference assignment rate for every scorer tried; the shipped one managed 0.717, the adopted one 0.945 |
| `af_accuracy_before_after.tsv` | allele-frequency error before and after that change, on simulated data |
| `real_data_af_shift.tsv` | the same change on real data: every position and length unchanged, maximum AF shift 0.0025 |

## Performance

| file | what it shows |
|---|---|
| `thread_scaling.tsv` | wall time and per-stage time from 1 to 20 workers; the pipeline is 66-86% parallelisable, so 8 threads is the practical limit |
| `msa_panel_size.tsv` | why `--msa-max-unique` defaults to 150 |
| `msa_parallel_timing.tsv`, `msa_parallel_verification.tsv` | the per-peak MSA parallelisation: byte-identical output, and a benefit small enough to be worth stating |

## Validation set

`validation_summary.tsv` — the seven fragment-analysis-confirmed ITDs across
five samples, previous-version calls against current ones. All lengths and
positions match exactly.
