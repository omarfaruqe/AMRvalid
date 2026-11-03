## Submission

This is a new submission.

AMRvalid evaluates antimicrobial resistance phenotype predictions made by
external models and tools, with methods that account for the population
structure of bacterial isolates. It does not train models.

## Test environments

* local: macOS (aarch64), R 4.6.2
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)
* win-builder: R-devel (2026-11-02 r90827 ucrt)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

  On win-builder the NOTE also listed "et" and "al" as possibly misspelled
  words in DESCRIPTION. They come from the reference "Collins et al. (2024)"
  and are spelled correctly.

## Notes for reviewers

* The example data in `inst/extdata` (48 KB) are derived from public NCBI
  Pathogen Detection data; the script that builds them is in `data-raw/`
  in the source repository and excluded from the package build.
* Functions with a `seed` argument set the random seed only when the user
  supplies one, and restore the previous random state on exit.
* `amr_report()` writes a file only when the user gives a `file` path;
  examples and tests write to `tempdir()`.
