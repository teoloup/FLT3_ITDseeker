# Recurrent observed haplotype rescue

Historical stable design and measurements at `fb9cae3`. This experimental branch
overrides the anchor defaults and identity filter as described in
[FLEXIBLE_FLANK_EXPERIMENT.md](FLEXIBLE_FLANK_EXPERIMENT.md).

Baseline source: `080b9a7`. Implementation branch: `fix/observed-haplotype-rescue`.

The two discovery failures share one rescue path. Before fitting the GMM, count
complete observed windows between exact outer anchors on every trimmed read.
This includes excluded and WT-assigned reads. After the ordinary GMM/DADA2/MSA
path finishes, add qualifying observed sequence seeds alongside resolved
consensuses. Thus an unresolved consensus can be rescued without emitting its N
or lowering consensus thresholds. Exact full alleles already represented by a
consensus are deduplicated. Every new allele still competes against WT and all
other accepted alleles in the existing final validation.

Discovery requires 20 identical observations, at least two per strand, median
whole-read Q20, unambiguous A/C/G/T, a complete net insertion within the requested
size bounds, exact reconstruction as one pure insertion, and at least 80% identity
to an adjacent WT copy. Whole-allele reconstruction rejects substitutions,
deletions and separate insertion gaps; it does not promote an isolated short gap
inside a longer event. Discovery recurrence and reporting AF are separate.
The existing defaults of 12 bp minimum size and 5% reporting AF remain unchanged.

One-base junction alternatives remain observed sequence hypotheses, rather than
an inferred claim of separate clones. HTML and the audit explain that uncertainty;
the VCF has a `SEQUENCE_UNCERTAINTY` flag. Competition can leave reads ambiguous or
split support between these sequences. Their individual AFs should not be read as
proven independent clone fractions. The final audit includes validated read counts
for candidates below the reporting threshold.

## Cost and reproducibility

Discovery uses string anchors and exact insertion decomposition, without an extra
full read alignment. Counting duplicate observations happens once per full window;
only recurrent seeds undergo decomposition. Extra accepted references increase the
number of comparisons in the existing final competition. A bounded cache (2048
unique sequences per worker chunk) reuses duplicate-sequence alignment metrics
while preserving each read's ID and strand.

`observed_rescue_discovery.json` records a replay on five available FASTQs, with
0.11–0.24 seconds per sample in that run. It found the 21 bp seed in 11528 (117
identical observations), both 45 bp alternatives in 13462 (57 and 31), and the 27 bp
seed in 10808 (26). The weaker 27 bp signal in 12156 did not qualify. With minimum
size deliberately lowered to 3 bp, 11551 produced only the coherent 15 bp seed
(2965 observations), with no 4 bp seed.

For matched full-pipeline runs, `benchmark_observed_rescue.py` uses the original
tracked source, the current source with rescue disabled, and the current source
with rescue enabled. Each runs alone with three threads, minimum size 12 bp,
reporting AF 0.005 (0.5%), and HTML enabled. This is a review threshold, not the
unchanged 5% default. Wall time includes extraction, trimming, GMM/DADA2/MSA,
plots, validation and reports. Linux `ru_maxrss` is the largest waited child
process's peak RSS, not the simultaneous aggregate memory of all workers or the
container. Single runs are measurements, not a stable performance guarantee.

Create the original-source archive from the repository root before benchmarking:

```sh
git archive --format=tar --output=review/observed_rescue_baseline.tar 080b9a7 Nano_ITDseeker.py itdseeker
```

In the project container with the repository mounted at `/workspace`:

```sh
cd /workspace
python review/benchmark_observed_rescue.py
PYTHONPATH=/workspace python review/check_observed_rescue.py
python review/summarize_observed_rescue.py
python -m unittest discover -s tests
```

`observed_rescue_results.json` holds compact commands, timings, memory, call
counts and candidate audit evidence. Full logs, reports, validation tables and
reads remain local under `review/observed_rescue_benchmark/` and are ignored by
Git. Prior fragment recheck summaries and reusable review scripts are committed;
sequencing inputs and per-read audit tables remain local.

To inspect lower-support candidates in an independently chosen review run, use
`--min-allele-frequency 0.001`; this changes reporting only, not discovery or
competition. Candidates below 0.005 already have validated counts in the audit,
so inspecting those counts does not require another alignment.

## Measured rescue cases

| Sample | Original source | Cache, rescue disabled | Cache and rescue | Change vs original | Largest process RSS, original / enabled |
| --- | ---: | ---: | ---: | ---: | ---: |
| 11528 | 154.86 s | 143.38 s | 156.40 s | +1.0% | 796.1 / 796.1 MiB |
| 13462 | 91.60 s | 76.72 s | 108.13 s | +18.0% | 754.4 / 753.9 MiB |
| 11551 | 136.41 s | 115.93 s | 89.05 s | -34.7% | 856.4 / 855.7 MiB |

The 11551 negative challenge retains exactly the same 15 bp call, 6657 supporting
reads, AF 33.508%, and 13210 WT reads. Every read's final classification is
unchanged, and there is no new candidate reference. Final alignment takes 23.99 s
with rescue disabled and 23.93 s enabled. The larger total-time difference comes
from run variability in other phases: GMM fitting alone is approximately 37 s
versus 19 s despite unchanged GMM settings and output. These single runs cannot
establish a guaranteed speedup. Its separate minimum-size-3 discovery replay
also excludes the misleading 4 bp isolated gap.

For 11528, the 21 bp seed has 117 exact observations: 115 were unassigned by the
GMM and two were assigned WT. Existing final competition recruits 278 validated
reads (137 plus, 141 minus), AF 2.592% of 10724 classified/validated reads. All
three original 159/36/33 bp calls remain. The major 159 bp count stays 1278;
36 bp changes from 130 to 128 and 33 bp from 70 to 58 as reads compete against
the newly available 21 bp reference. The AF denominator changes from 10460 to
10724; retaining a call does not imply its old AF should stay fixed.
An explicit read-ID comparison shows the 278 comprise 265 previously unclassified
or filtered reads, 12 previously supporting 33 bp, and one previously supporting
36 bp. One other former 36 bp read becomes unclassified. The 8982 WT reads retain
their WT classification.

For 13462, the unresolved H2 cluster supplies 55 of the 57 exact observations
of one 45 bp sequence and 30 of the 31 observations of the other; the remaining
two and one were GMM-unassigned. Final competition validates 137 and 74 reads,
respectively. Their combined size-level support is 211 / 20457 = 1.031%, without
establishing separate biological clones. The stronger sequence's AF is 0.670%
and passes the review reporting threshold; the weaker is 0.362%, below 0.5%, so
it remains in the validation audit rather than the VCF. Lowering reporting AF
to 0.001 would expose both without changing their competition.

The 48 bp call remains, but its support changes from 416 / 20424 (2.037%) to
238 / 20457 (1.163%) once the previously unavailable 45 bp references compete.
This is the intended reassignment of sequence evidence, not preservation of an
inflated original AF. Neither the recovered 21 bp nor the recovered 45 bp signal
would pass the unchanged default 5% reporting threshold.
Read-ID comparison shows that 176 former 48 bp reads now support the 45 bp
alternatives (114 and 62); 35 previously unclassified/filtered reads join them
(23 and 12). Two other former 48 bp reads become unclassified. All 20008 WT reads
retain their WT classification. The compact results include the full aggregate
transition matrix without storing read IDs.

The all-GMM-components-filtered integration check also passes: forcing
`--max-peak-sd 0` on 11528 leaves zero GMM-assigned reads, yet observed seeds still
reach WT competition and produce the 159/21/36 bp references, VCF and HTML.
The 21 bp candidate gains 279 reads in that artificial candidate set. The normal
run has 278 because it also includes the GMM-derived 33 bp candidate. This check
verifies discovery independence; it is not a replacement for the normal GMM/MSA
candidate set or a performance comparison. All 51 regression tests pass after the
final integration changes.

## Limitations

Exact anchors miss some reads and exact recurrence misses noisy or rare alleles.
This is a conservative supplement, not unrestricted recovery of every fragment
peak. Complex events that cannot be represented as one pure insertion do not
qualify. Read histograms, GMM curves and insertion pileups describe GMM-derived
observations. The GMM and insertion-size plots now add final-call panels showing
every reported sequence, including independent rescues and overlapping same-size
alternatives, without fabricating fitted peaks or read histograms. Independent
seed candidates are also recorded in the consensus table, audit, reference FASTA,
validation support table, VCF and HTML. Default validation uses all trimmed reads;
the optional `--validate-peak-reads-only` continues to restrict that read set.
