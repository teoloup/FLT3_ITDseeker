# Medaka consensus comparison: prepared, awaiting basecaller model

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
Medaka 2.2.2, CPU PyTorch 2.9.1, minimap2 2.30, samtools/htslib 1.22.
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

## Required information before running inference

The supplied test BAM headers contain generic RG metadata and minimap2/samtools
program records, but no original basecaller model. Medaka's actual model-resolution
command failed with `Failed to parse basecaller models from input file` because
RG `DS` is absent. See `runs/medaka/model_resolution.log`.

A matching Dorado/Guppy model name (and flow-cell chemistry) is required from the
user, including whether the merged BAMs combine different basecalling models.
No default model is silently assumed. No Medaka consensus accuracy result has
been produced yet. The run path below has been prepared, not exercised end to end.

```bash
# Run using the existing Linux validation Python; the script invokes isolated Medaka.
python review/benchmark_medaka.py --prepare-only
python review/benchmark_medaka.py --model CONFIRMED_MEDAKA_MODEL --case real72_sorter --case real72_msa --case real72_short_draft
```

The runner uses native Linux scratch storage, an explicit model, two requested
feature threads, a batch size of eight, and a five-minute per-case timeout.
It records the top-level argv in `commands.jsonl` and enables inherited Bash
tracing so the nested minimap2, samtools and Medaka commands appear in
`console.log`. Distinct model strings use distinct output directories. It
refuses to overwrite an existing run and checks for a single nonempty output
consensus; biological correctness remains a separate assessment.
