# Medaka consensus comparison: matched v5.0.0 model tested

All five prepared comparisons completed on 2026-09-23. Medaka changed zero
bases in every draft. It preserved the real 72 bp ITD and the exact synthetic
60 bp allele, but did not resolve the real drafts' terminal differences.
There is no demonstrated accuracy benefit on these cases. Keep polishing
experimental rather than adding it to the default production pipeline.

Medaka can polish a separated haplotype from any clustering method. It requires
reads and a draft full-amplicon sequence; it is not specific to amplicon_sorter
and does not itself separate mixed ITDs into haplotypes.

## What the paper actually tested

[The amplicon_sorter paper, section 3.2.1](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.8603)
used Medaka 1.4.3 as an optional additional polishing step. It improved two of
ten consensuses that initially differed from Sanger. This is evidence to test
polishing, not evidence that it always improves consensus or replaces clustering.

[Current Medaka documentation](https://github.com/nanoporetech/medaka#usage)
requires reads plus a draft and stresses selecting a model appropriate for the
basecaller. Current models are trained for Flye drafts; performance with our
short amplicon drafts must be measured rather than assumed.

## Supported experimental arrangements

1. Existing per-haplotype MSA draft + that haplotype's full reads -> Medaka.
2. Amplicon_sorter draft + the SAME full reads -> Medaka.
3. To remove MSA entirely: obtain a draft another way (e.g. amplicon_sorter's
   consensus or a representative-read/POA draft), then polish it with Medaka.
   A mixed length peak is not a sufficiently resolved haplotype by itself.

Polish the full insertion-containing amplicon, then align it to WT to recover
the insertion sequence and its boundary together. Do not polish isolated ITD
payloads without their flanks, and do not interpret polishing as AF estimation.
Keep competitive validation afterward. A draft containing only WT would also
make an unfair test of recovery of an insertion already detected by the pipeline.

## Installed and prepared

An isolated CPU environment is installed under `/tmp/flt3-medaka-review` in WSL:
Medaka 2.2.2, CPU PyTorch 2.9.1, minimap2 2.30, samtools/htslib 1.22.1,
and bcftools 1.22. The first attempt failed Medaka's dependency check before
inference because bcftools was missing. It was installed and the setup script
was corrected; that failed attempt is preserved as `real72_sorter_failed_preflight`.
The production environment and default backend were not modified.
The installer commands and output are in `runs/medaka/setup.log`.

`benchmark_medaka.py --prepare-only` successfully prepared five cases:

| Case | Reads | Draft |
|---|---:|---|
| real72_sorter | 589 | Amplicon_sorter per-peak consensus |
| real72_msa | Same 589 | Current MSA insertion reconstructed in WT context |
| real72_short_draft | Same 589 | Original 405 bp draft missing three primer-end bases |
| sim60_sorter | 2,000 | Amplicon_sorter 60 bp ITD consensus |
| sim60_msa | Same 2,000 | Current MSA 60 bp ITD draft |

The two synthetic drafts are identical. They are a preservation/regression
control, not an independent comparison of draft algorithms. The simulator does
not implement a particular ONT basecaller's learned error distribution, so
Medaka results on it cannot establish accuracy on real reads.

Prepared FASTA files, read counts and SHA256 hashes are recorded in
`runs/medaka/prepared_cases.json`. The original 72 bp drafts differ at the primer
end, not the insertion. Evaluate insertion sequence/length, anchor consistency,
terminal coverage, changes away from the draft, ambiguous bases and runtime.
Real-read sequence changes alone do not establish which draft is more accurate;
there is no independent full-amplicon truth for those reads.

## Confirmed model and provenance

The user supplied FLO-MIN114, SQK-NBD114-96, SUP at 400 bps and
Dorado software version 7.9.8, followed by the exact model
`dna_r10.4.1_e8.2_400bps_sup@v5.0.0`. The
[Medaka 2.2.2 mapping](https://github.com/nanoporetech/medaka/blob/v2.2.2/medaka/options.py)
selects `r1041_e82_400bps_sup_v5.0.0` for consensus. This explicit model was
used in every run. The supplied model is taken as the intended model for these
test reads; per-run provenance of the merged BAM cannot be independently checked.

The supplied test BAM headers contain generic RG metadata and minimap2/samtools
program records, but no original basecaller model. Medaka's actual model-resolution
command failed with `Failed to parse basecaller models from input file` because
RG `DS` is absent. See `runs/medaka/model_resolution.log`.

No default model was assumed. The run path below has now been exercised end to end.

```bash
# Run using the existing Linux validation Python; the script invokes isolated Medaka.
python review/benchmark_medaka.py --prepare-only
python review/benchmark_medaka.py --model r1041_e82_400bps_sup_v5.0.0
/tmp/flt3-medaka-review/venv/bin/python review/audit_medaka.py review/runs/medaka/runs/c1de56e24879
```

The runner uses native Linux scratch storage, an explicit model, two requested
feature threads, a batch size of eight, and a five-minute per-case timeout.
It records the top-level argv in `commands.jsonl` and enables inherited Bash
tracing so the nested minimap2, samtools and Medaka commands appear in
`console.log`. Distinct model strings use distinct output directories. It
refuses to overwrite an existing run and checks for a single nonempty output
consensus; biological correctness remains a separate assessment.

## Results

Results are in `runs/medaka/runs/c1de56e24879`, with per-case `summary.json`,
`coverage_audit.json`, exact command audit, Bash trace, BAM and probability HDF.

| Case | Reads | Output length | Changed bases | ITD | Seconds |
|---|---:|---:|---:|---:|---:|
| real72_sorter | 589 | 408 | 0 | 72 bp | 20.05 |
| real72_msa | 589 | 408 | 0 | 72 bp | 17.99 |
| real72_short_draft | 589 | 405 | 0 | 72 bp | 17.79 |
| sim60_sorter | 2,000 | 396 | 0 | 60 bp, exact truth | 18.10 |
| sim60_msa | 2,000 | 396 | 0 | 60 bp, exact truth | 17.64 |

The three real outputs have the same insertion sequence and local WT anchor 223.
The synthetic alignment reports anchor 119, an equivalent repeat representation
of the simulated insertion at 180: the reconstructed full allele equals truth
exactly. Compare full alleles, not just raw insertion coordinates.

These are actual inference results, not draft-only fallback: all input reads
mapped as primary alignments, probability HDF positions covered every draft base,
inference completed, and every gap BED was empty. No output contains ambiguous
bases. Median aligned-base depth was 588 on real drafts and 1,996 on synthetic
drafts. `audit_medaka.py` records coverage and HDF evidence independently.

Terminal behavior matters: the sorter draft starts `TTGT...`, while the MSA ITD
reconstructed in WT context starts `CTGT...`. Medaka retained both. Alignment to
the sorter draft gives 342 reads covering its first base, versus only six at the
MSA draft's first base, consistent with alignment/soft-clipping sensitivity at
the boundary. Those reference-dependent counts do not establish independent
truth. The shortened draft stayed three bases shorter despite using the same
589 reads. Polishing did not standardize ends or repair that truncation here.

## Recommendation

Medaka can run after either draft source, but these tests do not justify replacing
MSA or enabling polishing by default. Keep resolving haplotypes first and retain
the contextual ITD extraction and competitive read validation. The 72 bp cluster
split is addressed by clustering within its length peak (see
`AMPLICON_SORTER_REVIEW.md`), not by polishing two separate clusters afterward.
Before broader adoption, test drafts with known internal errors and additional
real haplotypes with independent truth. This experiment covers one real ITD and
one synthetic allele; it does not establish general accuracy or sensitivity.
