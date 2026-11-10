# amphid on nemits2 v0.1.0: evaluation scripts

The checks behind stephenrdoyle/nf-amphid-pipeline#10 (amphid support for this release) and this repo's issues #1 (the 40-candidate cap in `its2_classify.py`) and #2 (weak records that stop *H. placei* calls). All scripts are against `work/build_v0.1.0`, with test reads in `work/testdata/` and DADA2 results in `work/test/`; `work/` is not committed, so these must be rebuilt or fetched first. amphid is expected at `../../../nf-amphid-pipeline` (branch `nemits2-support`).

Scripts find their inputs and write their outputs relative to this folder, so they run from here. What they regenerate in place (`reads150/`, `reftop*.tsv`, `amph.*.tsv`) is git-ignored; the reported versions are in `results/`.

## amphid on the test sets

| script | what it does |
|---|---|
| `run_project.py <project> <150\|250> [args]` | amphid's `assign_db_pairs.py` on every sample of `work/testdata/<project>`, into `acc/<project>.<len>`. Cuts 2x150 from the raw 2x251 reads into `reads150/`, which isn't kept here (567 MB). Uses `nc_primers.tsv` and `--primer-offset 3`. Set `TAG=` to write into a separate folder. |
| `score.py acc/<set>` | read-weighted species / genus share, read fates, top taxa |
| `acc/` | amphid's per-sample outputs for the four sets (`clusters.tsv`, `assign_db.json`), as reported in PR #10 |

## amphid against the reference classifier

| script | what it does |
|---|---|
| `cmp_ref.py <clusters.tsv>...` | `classify_asvs.py` (as shipped, TOP=40) on amphid's own sequences, compared call by call |
| `ref_top.py <TOP> <work/test/prefix>...` | `classify_asvs.py` on the DADA2 ASVs with the candidate cap set to TOP; prints the gap-rule species share (issue #1) |
| `cmp_top.py acc/<set> <TOP>` | the same, on amphid's sequences, compared call by call with amphid |
| `exhaustive.py acc/<set> <n>` | amphid's span scoring of the n most abundant sequences against every record, compared with amphid's calls |
| `disagree.py acc/<set> work/test/<set>` | amphid's species calls against assignTaxonomy, IdTaxa and the reference on exactly matching DADA2 ASVs |
| `disagree_genus.py` | the same at genus level |
| `amphid_loo.py <build> <prefix> [--cap N] [-N N]` | the leave-one-out of `bin/its2_leave_one_out.py` (same queries and removal), scored by both the reference and amphid; `--cap 0 -N 5000` is uncapped |
| `v4_regress.py` | amphid only: vespa fixture outputs of `main` vs the branch, byte for byte. Needs `orig_bin/` (`git show main:bin/<file>` for every file in amphid's `bin/`) |

## Weak records (issue #2)

"Weak" means a partial or unmeasured record with one accession; "robust" means any other record.

| script | what it does |
|---|---|
| `redundancy.py` | records per species, and within-species nearest-neighbour distances for the top 10 species |
| `placei.py` | *H. contortus* records by distance to *H. placei*; the 7 within 3 edits |
| `weak_census.py` | weak records by whether their species has a robust record, and their distance to it |
| `narrow_sets.py` | writes `narrow_sets.json`: suspect (more than 5 edits from every robust record of the species) and redundant (within 1 edit) weak records |
| `seven.py` | which class each of the 7 *H. contortus* records falls in |
| `weak.py`, `weak_narrow.py` `acc/<set> <n>` | real-data species share with weak, suspect, or suspect + redundant records not counted as another species |
| `loo_weak.py`, `loo_narrow.py` | the leave-one-out (uncapped) under the same variants |

## results/

| file | contents |
|---|---|
| `loo_amphid.txt` | first leave-one-out comparison (TOP 60; an older tally that counted sub-family calls as family or higher) |
| `loo_capped.txt`, `loo_uncapped.txt` | the leave-one-out with and without the 60-candidate cap |
| `loo_weak.txt`, `loo_narrow.txt` | the weak-record variants (issue #2) |
| `reftop40.*.tsv`, `reftop1000.*.tsv` | `classify_asvs.py` on the DADA2 ASVs at TOP 40 and 1000 |
| `reftop1000.amph.*.tsv` | the same at TOP 1000 on amphid's sequences |

The scripts that run a multiprocessing pool guard it with `if __name__ == "__main__"`. On macOS a pool started from a script piped on stdin hangs, so run them as files.
