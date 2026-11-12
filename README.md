# AMRvalid

<!-- badges: start -->
[![R-CMD-check](https://github.com/omarfaruqe/AMRvalid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/omarfaruqe/AMRvalid/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24151080.svg)](https://doi.org/10.5281/zenodo.24151080)
<!-- badges: end -->

**Population-aware evaluation of antimicrobial resistance (AMR) phenotype predictions.**

AMRvalid evaluates predictions made by *other* tools and models (machine-learning
classifiers, rule-based genotype interpreters such as AMRFinderPlus, ResFinder,
Kleborate or AMRgen) against phenotypic reference results. It reports performance
per species-antimicrobial stratum, in both machine-learning terms (sensitivity,
specificity, AUC, MCC, Brier score) and clinical microbiology terms (categorical
agreement, very major, major and minor errors, essential agreement), with interval
estimates that remain honest when resistant counts are small.

AMRvalid does **not** train models, detect resistance genes, compute genomic
distances or interpret raw MICs into S/I/R. It uses the classes of the
[AMR](https://amr-for-r.org) package (`ab`, `mo`, `sir`, `mic`) for antimicrobial,
organism, category and MIC identifiers.

> Status: 1.0.0. The API is stable; any future changes will go through deprecation first.

## Installation

```r
# development version from GitHub
# install.packages("remotes")
remotes::install_github("omarfaruqe/AMRvalid")
```

## Example

```r
library(AMRvalid)

preds <- data.frame(
  isolate_id  = sprintf("iso%02d", 1:8),
  species     = "Klebsiella pneumoniae",
  drug        = "meropenem",
  truth       = c("R", "R", "R", "S", "S", "S", "S", "R"),
  probability = c(0.91, 0.75, 0.30, 0.10, 0.20, 0.65, 0.05, 0.88),
  model_id    = "my_model_v1"
)

x <- amr_import_predictions(preds)
amr_validate(x)
amr_evaluate(x, threshold = 0.5)
```

## Real example data

```r
# Public NCBI Pathogen Detection data: K. pneumoniae meropenem and
# ciprofloxacin phenotypes with rule-based predictions from AMRFinderPlus calls
x <- amr_import_predictions(amr_example_data())
amr_evaluate(x, weight_by = "genomic_cluster")
```

See `vignette("ncbi-kpneumoniae", package = "AMRvalid")` for a full worked example.

## Population-aware evaluation

```r
# Did lineages, dates or countries leak from training into test?
amr_validate_splits(x, train_manifest = train_ids, claimed_design = "group")

# How does performance change as test isolates get further from training?
# `dist` is a matrix or long table of genomic distances (Mash, SKA, SNPs...).
amr_distance_eval(x, distances = dist, train_manifest = train_ids, threshold = 0.5)

# Intervals that resample whole lineages rather than single isolates
amr_bootstrap(x, resample_by = "lineage", B = 2000, threshold = 0.5, seed = 1)

# Two models on the same isolates, with paired cluster-bootstrap differences
amr_paired_compare(x, models = c("old", "new"), resample_by = "lineage", threshold = 0.5)
```

## Roadmap

| Release | Scope |
|---|---|
| 0.1.0 | Data contract, validation, ML and AST metrics, small-denominator intervals (done) |
| 0.2.0 | Split-integrity checks, distance-to-training evaluation, cluster bootstrap, paired comparison, simulation study (done) |
| 0.3.0 | Calibration, prevalence adjustment, cluster weighting, provenance and TRIPOD+AI-aligned reporting, public example data and real-data vignette (done) |
| 1.0.0 | API review (done), CRAN submission (in progress); software paper next |

The full design is in [`docs/AMRvalid_Revised_Architecture_v3.docx`](https://github.com/omarfaruqe/AMRvalid/blob/main/docs/AMRvalid_Revised_Architecture_v3.docx).

## Licence

MIT © Md Omar Faruqe

## Citation

If you use AMRvalid, please cite it. `citation("AMRvalid")` gives the reference
in text and BibTeX form. The DOI
[10.5281/zenodo.24151080](https://doi.org/10.5281/zenodo.24151080) always points
to the latest release.
