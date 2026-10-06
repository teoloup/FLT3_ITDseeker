# amplicon_sorter evaluation

Recommendation: worth developing as an experimental per-peak backend, but not
ready to replace DADA2 by default. It improves separation of different-breakpoint
ITDs in the mixed challenge, while missing a three-substitution haplotype and
requiring stronger handling of duplicate alleles, partial output and runtime.

Status: experimental comparison; production remains DADA2.
Scope: existing synthetic BAMs and supplied `bam_data/test_bam` only.

## Version and setup

- Upstream: https://github.com/avierstr/amplicon_sorter
- Pinned checkout: `bf50ddcceb5d162d95a4fa64ad09c0b3e574001f` (script version 2025-10-09).
- Local clone: `review/tools/amplicon_sorter` (ignored third-party checkout).
- Added `edlib==1.3.9.post1` to the existing Linux validation venv, not production requirements.
- Four processes, no `--all`, no `--random`. `--allreads --maxreads 100000`
  preserves all eligible reads in these small inputs while retaining default
  batch comparisons. `--allreads` and `--all` are distinct switches.
- Same samtools/Cutadapt extraction as the current pipeline: MAPQ 20, retained
  primers, reverse-complement handling, minimum length 330, maximum length 666.
- `--save_fastq` preserves read IDs. Cluster FASTA headers are numeric indexes
  and include an extra consensus record, so they cannot be directly counted or
  mapped as original BAM read IDs.

Exact argument vectors and routing are in each experiment's `commands.jsonl`
and each input directory's `commands.jsonl` under `review/runs/amplicon_sorter`.
`manifest.json` records source hash and environment versions. Commands are logged
before execution; exit status alone is not accepted as evidence of completion.

## Synthetic clustering results

The mixed challenge contains WT plus four true ITDs: A (45 bp), B (45 bp at a
separate location), C (30 bp), and D (45 bp, three substitutions different from A).
There are 1,999 primer-trimmed reads. The results below use the UNMODIFIED upstream
script with its output on native Linux temporary storage, copied back afterward.

| Input/settings | Time | Clustering result |
|---|---:|---|
| Whole sample, defaults (96% consensus merge) | 118 s | WT; A+C+D merged; B |
| Whole sample, `--similar_consensus 99` | 110 s | WT; A+D merged; B; C |
| Current GMM 45 bp peak, `--similar_consensus 99.5` | 50 s | A+D merged; B |

All 1,999 reads were assigned once in both completed whole-sample runs. The
99% run produced exact WT, A, B and C full-amplicon consensuses. Counts were
1,093 WT, 487 A+D (382 A plus 105 D), 259 B, and 160 C. These are cluster counts,
not competitively validated allele frequencies. Defaults absorbed every C read
into the A+D cluster despite recovering an apparently clean A consensus.

For the direct DADA2 comparison, the exact same 704-read GMM peak was used:
amplicon_sorter produced an exact A consensus from 465 reads (364 A + 101 D),
and an exact B consensus from 239 B reads. DADA2 previously returned one ASV
for this peak in about eight seconds. Thus separation improved, at a runtime
cost, but D remained hidden even at the stricter merge threshold. These are
single-run observations; upstream stochastic consensus sampling has no CLI seed.

## Reliability findings

Initial four-process runs writing intermediates under WSL `/mnt/c` exited zero
with incomplete output: mixed default produced only WT; other runs produced no
clusters. These MUST NOT be interpreted as negative biological results.

A diagnostic copy changes only the top-level `except Exception: continue` to
print and re-raise. It exposed `IndexError` at `update_groups`, line 1755,
while reading a malformed temporary record `:0.984` instead of three fields.
The source creates a separate `Lock()` in each writer worker, which cannot
serialize writers against one another. Native Linux temporary storage avoided
this observed failure in the completed repeat runs. Filesystem interaction is
therefore implicated; the exact concurrency cause has not been proven.

Any integration must validate output completion and read IDs, preserve unassigned
reads, and fail loudly on malformed/partial output. Merely checking return code
or finding a consensus FASTA is insufficient. The review harness checks that no
`.group` files remain, a per-input consensus file exists, and IDs are valid and
not assigned twice. These are useful checks, not a comprehensive upstream health
check. No algorithmic patch was applied to the evaluated upstream script.

## Algorithm and integration implications

- Keep the initial GMM. It protects length-distinct minor ITDs from merging and
  limits runtime. Use amplicon_sorter on each candidate peak, rather than blindly
  replacing the entire pipeline with its whole-sample consensus output.
- Start an experimental per-peak backend with a stricter consensus merge threshold
  (99-99.5%), but do not treat that range as validated optimal settings.
- `--similar_consensus` controls a late merge. Initial group formation and the
  hard-coded 95% read/consensus threshold in `finetune()` also affect separation.
  Raising just the merge threshold cannot guarantee separation of near-identical
  haplotypes, as the missing three-substitution D haplotype demonstrates.
- The existing pipeline also has an independent sensitivity limit:
  `--min-subpeak-fraction 0.15`. D represents 101/704 = 14.35% of this parent
  peak, so even perfect D separation would be folded back by that default.
  An integration study must evaluate this guardrail as well as tool thresholds.
- It clusters sequences using edit distance; FASTQ qualities are retained for
  output, rather than entering a DADA2-style learned error model.
- Use its read assignments with our context-aware consensus and competitive
  validation. A high percentage assigned, or a clean consensus, does not establish
  correct haplotype separation or unbiased AF.
- Avoid `--all`. Even without it, grouping highly similar reads can be expensive.
  Never silently cap positive samples and report the resulting cluster fractions
  as full-sample AF. The negative-control capped run is a smoke test only.

## Reproduction

With the existing Linux validation environment activated and `MPLBACKEND=Agg`:

```bash
python review/benchmark_amplicon_sorter.py --sample sim_A --profile default --native
python review/benchmark_amplicon_sorter.py --sample sim_A --profile sc99 --native
python review/benchmark_amplicon_sorter.py --sample sim_A --profile sc995 --peak ITD_2 --native
python review/benchmark_amplicon_sorter.py --sample 13697_2runs_hg38_RG --profile sc99 --native
python review/benchmark_amplicon_sorter.py --sample 14417_2runs_hg38_RG --profile sc99 --native --maxreads 2000
python review/benchmark_amplicon_sorter.py --sample sim_no_wt --profile default --native
python review/replay_amplicon_clusters.py
python review/summarize_amplicon_sorter.py
```

The harness resumes existing summaries; use a fresh experiment folder to repeat
stochastic runs. It enforces a 600-second process-group timeout. The replay adapter
uses the measured 704-read assignments after checking exact input membership;
it exercises downstream production code without registering a permanent backend.
It is a cached-assignment proof of concept, not a timed fresh-tool pipeline run.

## Replay through the current pipeline

The review adapter replayed the two measured 45 bp clusters into the existing
pipeline and left the separate 30 bp peak intact. It verifies that all 704 IDs
match the current GMM peak before proceeding. Contextual consensus and competitive
validation then reported exact full alleles A, B and C:

| Allele | Simulated realised AF | DADA2 pipeline AF | Cluster-replay pipeline AF |
|---|---:|---:|---:|
| A, 45 bp | 19.10% | Missing | 23.707% |
| B, 45 bp | 12.95% | Missing | 12.358% |
| C, 30 bp | 8.00% | 14.160% | 8.512% |
| D, 45 bp / three substitutions | 5.25% | Missing | Missing |

Validated depth increased from 1,024 to 1,586. This is an improvement in this
challenge, not unbiased quantification: A remains inflated while D is absent.
Outputs and exact extraction commands: `runs/amplicon_sorter/pipeline_replay/`.

## Real positive sample

On all 7,690 eligible reads from test BAM 13697, `--similar_consensus 99` with
native Linux output completed in 347 seconds and assigned 7,677 reads once.
It recovered the expected 24, 30 and 72 bp insertion payloads. Cluster sizes:
6,020 WT-like, 672 with 24 bp, 317 with 30 bp, and two clusters of 238 and 430
with the SAME 72 bp insertion at the same boundary. The two 72 bp full-amplicon
consensuses had lengths 405 and 408 and edit distance three; their insertion
payloads were identical. This is evidence to collapse equivalent reconstructed
ITD alleles before competitive validation, rather than interpreting every
cluster as a separate ITD. We did not establish whether their flanking difference
is biological or technical. Raw cluster fractions are not validated AFs.

The replay's validated A-reference reads were 292 true A and 84 true D;
B-reference reads were 196 true B, C-reference reads 135 true C, and WT-reference
reads 879 true WT. This directly demonstrates the unresolved A/D cross-assignment.

## Negative control and runtime boundary

The first 2,000 eligible reads from negative test BAM 14417, using native Linux
output and 99% consensus merging, completed in 187 seconds. It assigned 1,997
reads to one WT-like cluster, with three selected reads unassigned. Its consensus
was 331 bp, five bases shorter than reference, with no insertion of reportable
size. This is a limited negative-control smoke test, not full-sample validation.
The original full 12,707-read WSL-mounted-output run exceeded the 600-second
limit and was terminated as a process group; it is recorded as a timeout, not a
negative result. This supports retaining the GMM and avoiding whole-sample
clustering of deep WT-dominated amplicons.

## Pure-ITD control

The noisy 60 bp ITD-only synthetic BAM completed on native Linux storage in
74 seconds: all 2,000 reads in one cluster and an exact 396 bp full allele
(reference 336 bp plus the true 60 bp insertion). On WSL-mounted output, the
same input instead completed with two 1,000-read clusters having identical
exact consensuses. This is another reason to use native scratch storage and
collapse identical reconstructed alleles before reporting.

All six initial native-storage experiments completed; none contained duplicate read
assignments or unknown read IDs. The exact replay alleles were independently
verified by reconstructing full alleles from VCF anchors and comparing with the
simulation truth, rather than accepting approximate lengths or sequence rotations.

No production code, default backend, or production dependency was changed.

## Follow-up: why the 72 bp insertion split

Read-level diagnosis and parameter tests are now complete. The 405 bp consensus
is exactly the 408 bp consensus with its first three bases (`TTG`) removed:

```text
408 bp: TTGTACCTTTCAGCATTTTGACG...
405 bp:    TACCTTTCAGCATTTTGACG...
```

The ITD sequence and all remaining bases are identical. Using the upstream
comparison function, the two drafts have 99.3% global identity and 100%
semi-global identity (`HW`, the mode used for consensus merging). They therefore
already pass 99% AND 99.5% consensus merging thresholds. The split is not justified
as two different ITDs by these results.

Both groups contain reads with both primer-end patterns. Among reads with the
exact short `TACCTTTCAGCA` prefix, all 128 in the first group and all 108 in the
second already start there in the sequences extracted from the supplied BAM.
The current Cutadapt step did not create this truncation. Because the BAM was
built from previously trimmed FASTQs, this does NOT establish whether upstream
trimming or sequencing caused those endpoints. The groups also differ strongly
in read orientation and input batches: the first contains batches 0, 6 and 7;
the second contains batches 1-5. This is inconsistent with two cleanly separated
insertion alleles and supports an early grouping/read-end effect.

The tool first creates broad gene groups, then processes each separately.
These clusters came from DIFFERENT gene groups (`trimmed_0` and `trimmed_1`),
so the later `--similar_consensus` step never directly compares them. The early
`--length_diff_consensus` gate defaults to 8%. It can retain separate broad groups
when their representative consensuses have different amplicon lengths. Native
whole-sample tests widening that gate support this explanation; the precise
intermediate random draft choices were not retained, so that part of the causal
chain is inferred from code, logs and the perturbation test.

| Experiment | Reads selected | Result | Seconds |
|---|---:|---|---:|
| Whole sample, 99%, default length gate | 7,690 | 24, 30, and two clusters with identical 72 bp insertions | 347 |
| Whole sample, 99%, `--length_diff_consensus 25` | 7,690 | One cluster per 24, 30 and 72 bp insertion | 247 |
| GMM 72 bp peak, 99%, run 1 | 589 | One 72 bp cluster, all 589 reads | 84 |
| GMM 72 bp peak, 99%, run 2 | 589 | One 72 bp cluster, all 589 reads | 84 |
| GMM 72 bp peak, 99.5% | 589 | One 72 bp cluster, all 589 reads | 84 |
| Mixed synthetic sample, 99%, length gate 25 | 1,999 | WT, A+D, B, C, unchanged from the earlier 99% result | 71 |

The wider-gate real run assigned 7,673 reads: 6,027 WT-like, 664 with 24 bp,
317 with 30 bp, and 665 with 72 bp; 17 were unassigned. This is not the same
population as the 589-read GMM-selected peak, so their counts are not directly
interchangeable. These remain clustering tests, not new validated VCF runs.

Recommendation: keep per-peak clustering, which removed this duplication in
three runs without relaxing global merging. `--length_diff_consensus 25` is a
useful tested whole-sample alternative here, not a universal setting: it affects
early broad grouping and also late comparisons, and the code uses a permissive
80% early identity check when this setting exceeds 8%. Preserve distinct ITDs
and deduplicate only equivalent reconstructed insertion alleles.

Evidence: `runs/amplicon_sorter/72bp_split_diagnosis.json`, individual run logs,
and `diagnose_72bp_split.py`. To reproduce the parameter perturbation:

```bash
python review/benchmark_amplicon_sorter.py --sample 13697_2runs_hg38_RG --profile sc99 --native --length-diff-consensus 25
python review/benchmark_amplicon_sorter.py --sample 13697_2runs_hg38_RG --profile sc99 --peak ITD_2 --native
python review/benchmark_amplicon_sorter.py --sample 13697_2runs_hg38_RG --profile sc99 --peak ITD_2 --native --replicate 2
python review/benchmark_amplicon_sorter.py --sample 13697_2runs_hg38_RG --profile sc995 --peak ITD_2 --native
python review/benchmark_amplicon_sorter.py --sample sim_A --profile sc99 --native --length-diff-consensus 25
```
