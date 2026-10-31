# AMRvalid

<!-- badges: start -->
[![R-CMD-check](https://github.com/omarfaruqe/AMRvalid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/omarfaruqe/AMRvalid/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
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

> Status: early development (0.0.0.9000). The API will change.

## Installation

```r
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
validate_amr_data(x)
amr_species_eval(x, threshold = 0.5)
```

## Roadmap

| Release | Scope |
|---|---|
| 0.1.0 | Data contract, validation, ML and AST metrics, small-denominator intervals |
| 0.2.0 | Split-integrity checks, distance-to-training evaluation, cluster bootstrap, paired comparison |
| 0.3.0 | Calibration, prevalence adjustment, provenance and TRIPOD+AI-aligned reporting |

The full design is in [`docs/AMRvalid_Revised_Architecture_v3.docx`](docs/AMRvalid_Revised_Architecture_v3.docx).

## Licence

MIT © Md Omar Faruqe
