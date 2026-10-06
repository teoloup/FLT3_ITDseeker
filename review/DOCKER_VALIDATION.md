# Docker test image, 2026-09-23

Built, validated and published: `teoloup/hematology-aml-flt3-itd:asv-fixes-20260923`
(linux/amd64). After the user logged in, the push retry succeeded. A separate
`docker buildx imagetools inspect` request confirmed the published tag and digest.

Published manifest digest:
`sha256:d9ec387724d0f44045919dadcdd888107becc6b8f8107807f383d96f073bac5e`.
The published digest matches the locally validated image.
Runtime Python/R sources are from `d8540b6`; Docker packaging and targeted GMM
regressions were committed in `e6f4381`. No main merge is needed to test this image.

## Reported failure

The old image merged nearby length components, assigned zero reads, and crashed
when plotting an empty array. Current code already includes mixture total
variance (within-component variance plus between-component mean variance) and
an empty-input plot guard. Averaging child SDs can create an unrealistically
narrow merged peak around a fractional length, rejecting every integer-length
read. This explains a plausible path to the supplied log; sample 11446 itself
was not available among the supplied test BAMs, so its exact outcome is untested.

New regression tests use discrete peaks at 330, 331, 332, 333 and 336 bp and verify
that the merged SD reflects their spread and all 200 reads remain assigned.
Another test verifies that plotting an empty read table returns without crashing.

## Image and validation

`.dockerignore` allows only build inputs and runtime sources, excluding BAMs,
metadata, review artifacts, virtual environments and Git history. The image
contains a `Nano_ITDseeker` launcher as well as its existing Docker entrypoint.
Amplicon_sorter and Medaka remain review experiments, outside the image.

The dependency checks passed: Python 3.10, samtools 1.16.1, Cutadapt 5.2,
DADA2 1.38.0, numpy 1.26.4, Biopython 1.88, pysam 0.24.1 and sklearn 1.7.2.
MUSCLE source resolved to `08c9a4294e7dd4d1b93673ea5ddc48ce977e956a`.
The Dockerfile uses dependency constraints rather than a fully frozen lockfile;
use the published digest to reproduce this particular image after publication.

All 20 regression tests passed both locally and inside the image. The launcher
help check passed. Three complete container runs with `--min-allele-frequency
0.02 --html-report` matched previously validated variant calls exactly:

| Input | Expected and observed calls | AF | DP |
|---|---|---:|---:|
| sim_no_wt | Exact synthetic 60 bp allele | 100% | 1536 |
| 10808_hg38_RG | Exact known 45 bp allele | 41.209% | 13531 |
| 14417_2runs_hg38_RG | No ITD calls | - | - |

HTML reports and persisted samtools/Cutadapt/DADA2 command audits were checked.
Only synthetic and `bam_data/test_bam` inputs were used. The old published image
and `latest` tag were not changed. With user approval, unused Docker build cache
older than 24 hours was pruned, reclaiming 1.095 GB; no images or volumes were removed.

## Commands and evidence

Run from the repository root:

```powershell
docker build --platform linux/amd64 --build-arg VCS_REF=d8540b6 -t teoloup/hematology-aml-flt3-itd:asv-fixes-20260923 .
./review/test_docker_image.ps1
./.review-venv/Scripts/python.exe review/check_docker_results.py
docker push teoloup/hematology-aml-flt3-itd:asv-fixes-20260923
```

The test runner records Docker argv in `runs/docker_20260923/docker_commands.jsonl`.
Build/push logs, dependency listing, unit-test output, per-sample outputs and
`checked_results.json` are in that same directory. Pipeline command audits are
in each sample output folder. The failed push log and successful `push_retry.log`
are retained there.

The user's test command can now be run as:

```bash
apptainer exec \
  docker://teoloup/hematology-aml-flt3-itd:asv-fixes-20260923 \
  Nano_ITDseeker \
  -b 11446_2runs_hg38_sorted.bam -o flt3_asv_fixes -s 11446_2runs \
  -g hg38 -t 16 --min-allele-frequency 0.02 --html-report
```

Use `singularity` instead of `apptainer` if that is the runtime installed on the
cluster. The launcher was verified with Docker; Apptainer conversion was not
available to test here. A separate output directory preserves the old run.
