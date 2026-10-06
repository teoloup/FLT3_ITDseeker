# Benchmarks

Measurements behind the non-obvious defaults in this pipeline. Recorded here
rather than under `test_output/`, which is gitignored, so the evidence for a
default travels with the code that uses it.

Regenerate any of these with `simulate_itd_data.py` + `evaluate_haplotypes.py`;
the validation BAMs referenced by the real-data tables are not in the repo.

## Haplotype splitting

The pipeline now splits each GMM length peak with DADA2 only. The other methods
in these tables have been removed from the code, and their results are kept as
the evidence for that choice. `gmm2pass` is the earlier second pass, which fitted
a local `k=1` vs `k=2` length GMM inside each peak. Only the `dada2` rows can be
regenerated with this version.

| file | what it shows |
|---|---|
| `method_comparison.tsv` | gmm2pass / isonclust on simulated scenario A, where three of four ITDs share a length (dada2 and amplici on the same scenario are in `min_child_fraction_experiment.tsv`, at the default 0.15) |
| `real_data_method_comparison.tsv` | the same four methods on the real validation samples |
| `real_11531_four_methods.txt`, `real_13697_14219_four_methods.txt` | the raw per-peak output behind that table |
| `isonclust_parameter_sweep.tsv` | why the removed isONclust backend used an aligned threshold of 0.80 rather than isONclust's own 0.4 |
| `min_child_fraction_experiment.tsv` | whether the size guardrail was what blocked the hardest haplotype (it was not, for two of three tools) |
| `haplotype_method_chaining.tsv` | gmm2pass alone against gmm2pass followed by dada2, on a simulated 80/20 mixture and on 13697 (chaining has since been removed) |
| `unbalanced_mixture_detection.tsv` | `max_minor_fraction` on clean real peaks (0.024 to 0.101) against an 80/20 mixture (0.217) that leaves no `N` in the consensus. This is why DADA2 runs on every peak, and the basis for the 0.15 mixed-consensus warning |
| `truth_haplotypes_A.tsv`, `simA_*_per_haplotype.tsv` | the simulated truth set and per-haplotype scoring |

**Headline.** On simulation the sequence clusterers recovered 3 of 4 ITDs,
against 1 of 4 for the length-based baseline. On the real samples the baseline
was already correct on all three, `dada2` matched it, and `isonclust` and
`amplici` each lost an ITD. The two results differ because the simulation was
built to defeat length-based separation, while the real co-occurring ITDs in
this set differ in size, which the GMM already handles. `dada2` was the only
method that matched the baseline on every real sample and also beat it on
simulation. It is now the default and only backend.

## Competitive validation scoring

| file | what it shows |
|---|---|
| `validation_scoring_diagnosis.tsv` | correct-reference assignment rate for every scorer tried; the shipped one managed 0.717, the adopted one 0.945 |
| `af_accuracy_before_after.tsv` | allele-frequency error before and after that change, on simulated data |
| `real_data_af_shift.tsv` | the same change on real data: every position and length unchanged, maximum AF shift 0.0025 |

## Strand bias

| file | what it shows |
|---|---|
| `strand_bias_depth_sensitivity.tsv` | the same 3-point strand difference in 10808 (odds ratio 0.885) gives Fisher p = 0.0005 at full depth and p = 0.71 at 2% depth, which is why the report also requires an odds-ratio skew before flagging bias |

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
