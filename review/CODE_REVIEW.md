# Code review: feat/asv-haplotype-clustering

Reviewed commit: `26cdfc1` (2026-09-21). Algorithm source files were not changed. This review combines code inspection with executable synthetic checks; it is not a real-BAM validation or an estimate of clinical sensitivity/specificity.

## Findings, in priority order

1. **High: multi-reference z-score accepts tied candidates.** `Pairwise_aligment_toolkit.py:202-207,234-241`; `Helper_functions.py:217-230`.

   The standardization itself is mathematically correct: `(p - mean(p)) / population_sd(p)`. The problem is interpreting the winning value as sufficient evidence. Multi-reference classification ignores the computed best-versus-second delta. Scores `[0.9,0.9,0.1]` produce an ambiguous tie, but adding two irrelevant references gives `[0.9,0.9,0.1,0.1,0.1]`, z=1.22474 for BOTH winners, and an ITD-supporting call despite delta=0. A tie between ITDs is broken by alias ordering. Thus inferred haplotype counts and allele frequencies can depend on how many other candidate references exist. Softmax with beta=500 also makes small score differences look decisive; its values are not calibrated probabilities.

   Recommendation: require explicit winner separation in every case, reject exact ties, and calibrate the score/margin criteria against truth-labelled data. Do not interpret this z-score as a significance test or probability of correctness.

2. **High: an ITD-only sample can be labelled WT and reported negative.** `GMM_peaks.py:246-259`; `Nano_ITDseeker.py:508-513`.

   The closest peak is always labelled WT, with no maximum permissible distance from the configured 336 bp length. A synthetic population of 200 reads of length 396, fitted as one component, becomes WT with putative ITD size 0. The default downstream path excludes WT and has no ITD candidate to validate. The same assumption shifts all other predicted sizes when the true WT component is absent or filtered out.

   Recommendation: allow WT to be absent; use the known reference length as the fallback size baseline and verify candidate WT identity rather than assuming the nearest component must be WT.

3. **High: `--cluster-wt-peak` creates aliases downstream validation rejects.** `haplotype_split.py:161,203-205`; `Helper_functions.py:226-236,299-340`.

   WT children become `WT_H1`, `WT_H2`, etc. They are marked non-WT and sent through candidate construction, but support classification recognizes variants by an ITD name and gap validation accepts only exact WT or ITD-prefixed aliases. A direct probe of `WT_H2` produces `WT-supporting` followed by `unknown_alias:WT_H2`. Hidden ITDs recovered by this option cannot be reported through this path.

   Recommendation: separate biological reference type from display aliases; retain a genuine WT reference and explicitly identify WT-derived variant children.

4. **High: breakpoint validation ignores deletions within the proposed ITD.** `Helper_functions.py:61-84,305-327`.

   Only positive query-minus-target jumps are counted. Missing bases in the read therefore do not contribute to `max_gap_bp`. With the actual FLT3 sequence, a proposed 60 bp duplication and a read containing only 48 bp of that duplication, the real competitive alignment and validation functions produce `ITD_1`, identity=1.0, adjusted score=0.940312, `ITD-supporting`, and `pass`. Gapless-block identity hides the deletion. This may attribute a shorter haplotype to a longer candidate when the correct reference is absent.

   Recommendation: examine both insertion and deletion lengths, plus breakpoint-spanning coverage, using alignment coordinates. Decide and test an explicit tolerance for ONT deletions rather than silently ignoring all deletions.

5. **Medium: probabilistic assignment uses pre-merge component confidence.** `GMM_peaks.py:196-222`.

   After merging close components, predict_proba/hybrid still threshold the largest ORIGINAL component posterior. Probabilities for children belonging to the same merged peak are never summed. In a real fitted two-component synthetic example that merges to one peak, only 512/600 reads survive a 0.85 threshold, although every read has total posterior 1.0 for the sole merged peak. This can distort support and AF. The default manual assignment is unaffected by this particular bug.

   Recommendation: aggregate component posteriors by merged ID before selecting and thresholding the winning merged peak.

6. **Medium: incomplete cluster assignments leave inconsistent read membership.** `haplotype_split.py:116-125,175-182`.

   Unmapped reads are omitted from cluster membership. Once a parent splits, its subset is removed, but unmapped rows retain the now-nonexistent parent alias in reads_df. A 100-read parent with 80 assigned reads becomes 80 subset reads while reads_df still holds 100, including 20 stale aliases. Child fractions also use mapped-read count rather than the documented full parent count. AmpliCI can trigger this through short/ambiguous-read removal or unassigned reads; DADA2 excludes NA assignments too. These reads still reach competitive validation, so this is not simply a claim that all final AF denominators lose them: it specifically breaks subset consistency, discovery support, reported cluster fractions, and chained splitting.

   Recommendation: preserve unmapped reads in an explicit residual group or an explicit unassigned state, and define a consistent parent denominator.

7. **Medium: DADA2 ignores `--threads`.** `haplotype_split.py:392-393`; `dada2_cluster.R:57,66`.

   The Python-to-R command does not pass a thread count, and both R calls use `multithread=TRUE`. The default backend can therefore use available cores even for `--threads 1`. Pass the requested integer to R and both DADA2 calls. This is confirmed by wiring inspection; R was not executed here.

8. **Medium: `--max-itds-to-detect` does not cap final haplotypes.** `Nano_ITDseeker.py:86,404-415,472-474`; `haplotype_split.py:153-174`.

   The option caps first-pass GMM components only. Subsequent splitting/chaining can create more variants, with no final cap. Either document/rename it as a first-pass component limit or implement an explicit final selection policy. This finding is from code inspection.

9. **Lower priority: low-depth input can crash GMM model selection.** `GMM_peaks.py:88-98`.

   The component search is not capped by read count. Three reads under defaults fail at k=4 with `Expected n_samples >= n_components`. One read also cannot be fitted by the estimator. Use explicit minimum-depth handling and cap feasible component counts; insufficient evidence should not be reported as a biological negative.

## Checks that did not reveal a bug

- The hard-coded 336 bp sequence exactly matches the genomic PLUS-strand hg38 interval chr13:28033881-28034216 retrieved from [UCSC's sequence API](https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr13;start=28033880;end=28034216). FLT3 being a minus-strand gene does not invert this stored sequence. For internal insertion boundary p, `start+p-1` and `wt[p-1]` are consistent genomic anchor conventions. hg19 was not independently fetched.
- The VCF writer now prefixes the insertion payload with its REF anchor.
- Merged-peak variance includes both within-component variance and between-component mean spread.
- Minimum/maximum ITD size and MSA threshold CLI arguments are passed into the relevant extraction/consensus calls. Cutadapt's 10% upper-length allowance precedes the exact insertion-size filter; that allowance is not itself a disconnected parameter.
- The alignment-block review was cross-checked against [Biopython's alignment documentation](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html).

## Evidence and limits

Run from the clone directory:

```powershell
.\.review-venv\Scripts\python.exe review\review_checks.py
.\.review-venv\Scripts\python.exe review\alignment_checks.py
```

Saved outputs: `review_checks_output.txt`, `alignment_checks_output.txt`; reference response: `reference_orientation.json`. The scripts exercise production functions; they are reproduction probes, not a complete regression suite. Environment: Python 3.9.5, NumPy 1.26.4, pandas 2.3.2, SciPy 1.13.1, scikit-learn 1.6.1, Biopython 1.85. The local `.review-venv` uses system site packages with local NumPy/scikit-learn installations.

No BAM files were present. Full pipeline execution remains untested: samtools, cutadapt and Rscript were not found on PATH, pysam/pymuscle5 were absent from the checked Python environment, and the Docker daemon was unavailable. No source fixes were applied.

For the BAM phase, prioritize a negative control, a clear single ITD, little/absent WT, and two same-length or nearby-length ITDs. Compare discovery and validation read retention, candidate ties, breakpoint insertions/deletions, consensus/reference identity, and AF against orthogonal truth where available. Quantification currently conditions on reads retained through GMM and validation; BAM testing should measure whether retention differs between haplotypes. Also distinguish insufficient/no coverage from a tested negative, since several current early-exit paths write no-ITD outputs.
