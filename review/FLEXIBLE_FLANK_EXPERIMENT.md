# Flexible-flank discovery experiment

Branch: `experiment/flexible-flank-rescue`, based on stable `fb9cae3`.
This is separate development/test work. The stable branch and main are unchanged.

## Changes

This branch defaults to 20 bp locator anchors shifted 10 bp toward both WT ends,
with at most three substitutions **per anchor** (15% of a 20 bp anchor). Indels
are not part of the anchor search. Configuration flags are:

```text
--observed-anchor-length 20
--observed-anchor-shift 10
--observed-anchor-max-mismatches 3
```

`--observed-anchor-shift 15` tests the second requested position. All coordinates
below are 0-based, half-open on the 336 bp WT reference:

| Configuration | Left anchor | Right anchor | Exact observed interior |
| --- | --- | --- | --- |
| Stable: 15 bp, shift 0, no mismatches | [25,40) | [300,315) | [40,300) |
| Experiment: 20 bp, shift 10 | [15,35) | [305,325) | [35,305) |
| Experiment: 20 bp, shift 15 | [10,30) | [310,330) | [30,310) |

The outer window start/stop are `25 - shift` and `len(WT) - 21 + shift`.
Settings must keep both anchors inside WT and leave a nonempty interior.
Mismatch allowance must be nonnegative and smaller than anchor length.

The search partitions each anchor into `k+1` exact pieces for `k` permitted
substitutions. Any qualifying hit must contain at least one exact piece; candidate
locations then undergo a bounded Hamming-distance check. Duplicate read sequences
share cached locator results. Discovery does not align every read and introduces
no new dependency.

All matching left/right locations are considered in the required order, without
overlap. A read must have **one unique plausible ordered pair**. Multiple plausible
pairs are rejected even when one is exact and another has substitutions; the code
does not choose the first/last match or a preferred low-mismatch copy. Ambiguity,
missing anchors, and wrong-order outcomes are counted explicitly.

Only the matched locator anchors are replaced with their WT sequences before
reconstructing a seed. The complete observed interior stays exact. This allows
flank-error-bearing reads to support an insertion without silently repairing an
interior substitution or a complex indel. It also deliberately treats any genuine
anchor-region SNP as a locator difference, rather than calling that flank variant.
The assay remains an insertion caller, not a general flank-variant caller.

The audit records exact interior observations, exact-flank observations, normalized
flank observations, total/max substitutions per anchor, and anchor configuration.
The compatibility field `exact_reads` now counts identical **interiors**, not
necessarily exact complete observed flanks. `exact_interior_reads`,
`exact_flank_reads`, and `normalized_flank_reads` make the distinction explicit.

Adjacent duplication identity is retained as an annotation and no longer has an
80% acceptance cutoff. A recurrent low-identity pure insertion can therefore be
proposed. This alone does not establish a tandem-duplication mechanism. Identity
and flank-normalization support remain in consensus candidate metadata and audit
TSV; VCF fields are `DUP_IDENTITY` and `FLANK_NORMALIZED_READS`. These describe
seed evidence, not all reads recruited by competitive validation.

Recurrence >=20, at least two observations per strand, median whole-read Q20,
unambiguous exact interior, complete net insertion size, reconstruction as one pure
insertion, candidate deduplication/cap, WT competition, gap checks, and reporting
AF remain in place. Default minimum insertion size remains 12 bp and reporting AF
remains 5%. The same-size junction-uncertainty annotations remain in place.

## Retained-data comparison

`flexible_anchor_comparison.json` and `.tsv` compare the original stable module,
experimental 15 bp exact anchors, and both 20 bp/three-substitution shifts.
All configurations use the same already-trimmed, >=330 bp FASTQ records and
recurrence threshold. The 11551 negative challenge deliberately uses minimum
insertion size 3 bp; other samples use 12 bp. No GMM fitting, denoising, or final
read alignment is repeated for this discovery comparison.
The stable locator used first/last matches without an ambiguity audit; its TSV
ambiguity field is therefore blank rather than a claim of zero ambiguous matches.

| Sample | Stable coverage | Shift 10 coverage | Shift 15 coverage | Shift 10 / 15 discovery time | Ambiguous pairs, shift 10 / 15 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 11528 | 85.40% | 92.40% | 91.52% | 0.282 / 0.297 s | 0 / 1 |
| 13462 | 87.19% | 93.58% | 92.82% | 0.401 / 0.353 s | 1 / 1 |
| 11551 | 87.39% | 92.70% | 90.71% | 0.588 / 0.500 s | 0 / 0 |
| 10808 | 86.99% | 92.86% | 91.37% | 0.445 / 0.429 s | 1 / 1 |
| 12156 | 85.28% | 92.18% | 91.15% | 0.412 / 0.365 s | 0 / 0 |

Both experimental shifts retain exactly the same qualifying full-allele seed IDs
as stable discovery in these five samples. There are **no new or lost qualifying
seeds**. Removing the identity cutoff also adds no seed in these data; a synthetic
low-identity insertion regression verifies that the cutoff is actually removed.

Better locator coverage does not automatically increase coherent seed recurrence:
the exact interior now spans 10 or 20 additional bases, where sequencing errors
remain unrepaired. The 11528 21 bp seed has 117 stable observations versus 114
(shift 10) and 109 (shift 15). The two 13462 45 bp seeds have 57/31 stable
observations, 57/31 with shift 10, and 54/30 with shift 15. Both remain discoverable.
11551 proposes only its coherent 15 bp seed, with no spurious 4 bp seed even at
minimum size 3. The weak 12156 27 bp signal still does not reach recurrence 20.

Shift 10 is the provisional branch default because it has higher measured coverage
and stronger retained seed support than shift 15 on these data. These measurements
do not establish a clinical sensitivity or specificity gain. Timing is a single
discovery-only run per configuration, not a full-pipeline performance claim.

## Full integration check

One complete 13462 run with the experimental shift-10 defaults, three threads and
review reporting AF 0.005 completed successfully. It produces exactly the same
candidate reference sequences and all 22115 per-read classifications as the saved
stable run. The 45 bp alternatives retain 137 and 74 validated reads, and the
48 bp call retains 238. The weaker 45 bp sequence remains below reporting AF 0.005.
The 57/31 exact-interior seed supports comprise 54/29 exact-flank observations
and 3/2 anchor-normalized observations. The 48 bp seed has 92 observations,
89 with exact flanks and three normalized.

The exported 48 bp and stronger 45 bp VCF calls carry `DUP_IDENTITY=1.0000` and
`0.9778`, respectively, and `FLANK_NORMALIZED_READS=3`. The main HTML keeps only essential sequence uncertainty; normalization
and identity details remain in audit TSV/VCF. `flexible_anchor_integration.json` records
the aggregate verification, commands and single-run metrics; raw read IDs and
sequences are excluded. The integration took 112.29 s, with largest-process RSS
753.2 MiB; this single run is not a matched full-pipeline performance benchmark.

All 62 regression tests pass, including exhaustive-match comparison, coordinate
bounds, three substitutions per anchor, ambiguity rejection, unchanged interior
bases, and acceptance of a synthetic low-identity insertion. Export assertions
also verify the new VCF annotation fields.

Reproduce inside the existing project container mounted at `/workspace`:

```sh
git archive --format=tar --output=review/flexible_anchor_baseline.tar fb9cae3 itdseeker/observed_haplotypes.py
python review/compare_flexible_anchors.py
python -m unittest discover -s tests
```

The saved-data integration verifier is `python review/check_flexible_anchor_integration.py`.

The archive contains trusted tracked stable source and is ignored by Git. Raw
FASTQs and per-read audit data remain local. The experiment comparison outputs
contain aggregate counts and seed hashes, not read IDs or read sequences.
