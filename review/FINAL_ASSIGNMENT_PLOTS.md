# Main plots from final competitive assignments

Report cleanup: only the full read-length view plus its ITD zoom is embedded.
The redundant net-length image remains a saved diagnostic artifact. Generic
discovery/anchor/identity annotations remain in audit TSV/VCF, while the HTML
shows a short sequence-uncertainty note only when relevant. Sequence and technical
images are collapsed, and plot legends use readable size/AF/support labels.

The previous plots kept GMM-derived candidate distributions and added final call
markers. That could leave a rejected/mixed H2 cluster labelled "not reported"
while rescued 45 bp sequences were reported. The global insertion histogram also
kept pre-validation H1/H2 support rather than the reads assigned to final alleles.

Main plots are now regenerated at the end from the final validation table:

1. Require `valid_support=True`, `filter_reason=pass`, and an explicit supporting
   call. Keep WT support and ITD support for the AF-filtered reported aliases only.
2. Join these read IDs to actual trimmed read lengths. Missing lengths, duplicate
   passed assignments, or disagreement with reported support counts raise an error.
3. Draw actual 1 bp histogram bins separately for each final allele. The full
   read-length panel includes WT; its zoom contains reported ITD reads only.
   Same-size alternatives have separate colors and line styles, with size A/B,
   AF and supporting read counts in the legends. Record IDs remain in collapsed sequence details. No Gaussian is invented or fitted.
4. The second plot uses **actual trimmed read length minus WT reference length**.
   This saved diagnostic is a measured net read-length change, not an extracted insertion event.
   Dotted markers show exact consensus-called insertion sizes. Sequencing indels
   and trimming can move observed histogram peaks away from those sizes.

The legacy filenames `*_itd_gmm_fit_plot.png` and `*_itd_size_distribution.png`
remain for report compatibility; their captions and axis labels describe the new
data accurately. The genuine earlier-stage GMM plot is stored separately as
`*_candidate_gmm_fit_plot.png` and appears only in a collapsed diagnostic section.
Candidate-stage insertion observations/pileups remain explicitly labelled as
earlier evidence, not final support. Neither calling nor AF calculation changes.

## Verified reports

The same report the user opened was regenerated in place:

`review/observed_rescue_benchmark/13462/merge_review_af001/results/13462_itd_report.html`

Its final distributions use 20008 WT reads, 137 reads for 45 bp A, 74 for 45 bp B,
and 238 for 48 bp, at reporting AF 0.001. The two 45 bp forms remain uncertain
sequence alternatives rather than established distinct biological clones.

The 11528 example at AF 0.005 was refreshed at:

`review/observed_rescue_benchmark/11528/enabled/results/11528_itd_report.html`

It uses WT 8982 and final ITD support 1278 (159 bp), 278 (21 bp), 128 (36 bp),
and 58 (33 bp). All histogram values come from retained measured FASTQ lengths,
matched by read ID to the saved validation table. No read alignment was repeated
and the VCFs were not rewritten. Aggregate source/count/bin evidence is in
`final_assignment_plot_check_13462.json` and `final_assignment_plot_check_11528.json`.

Both 13462 main figures were visually checked. Automated tests additionally
verify the actual histogram input arrays: a 45 bp call supported by measured
377/380 bp reads produces 377/380 bp read-length values and 41/44 bp net values,
not a fabricated 381/381 bp or 45/45 bp spike. Invalid, filtered and unreported
candidate rows do not appear. All 66 regression tests pass.

Refresh the exact retained report inside the project container:

```sh
python review/refresh_final_assignment_report.py
python review/refresh_final_assignment_report.py --sample 11528 --mode enabled --min-af 0.005
python review/render_called_peak_plots.py
```

The standalone linked preview page is also regenerated from the same final
assignment sources. Raw read IDs and sequences are not included in the committed
aggregate evidence; the retained per-read files and reports stay local.
