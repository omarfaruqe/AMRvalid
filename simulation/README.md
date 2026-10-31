# Simulation study

Reproducible simulation behind the AMRvalid methods paper. It is excluded
from the R package build (see `.Rbuildignore`).

## Questions

1. **Optimism.** How much does each validation design overstate performance
   on lineages never seen in training?
2. **Interval coverage.** How often do bootstrap intervals that resample
   isolates, versus whole lineages, contain the true performance on new
   lineages?

## Design

Each replicate calls `amr_simulate()` with 40 development lineages of 25
isolates and 200 external lineages, which serve as the ground truth for
performance on new lineages. Resistance depends on a causal marker (log-odds
effect 3) and a clonally inherited lineage effect whose standard deviation
defines three scenarios: weak (0.5), moderate (1.5) and strong (2.5).

Two fixed predictors are evaluated under two designs:

| Predictor | Behaviour |
|---|---|
| `knn` | Proportion resistant among the 5 nearest training isolates; memorises lineages |
| `marker` | Proportion resistant among training isolates with the same marker status; learns the mechanism |

| Design | Test set |
|---|---|
| `random` | 20% of isolates, chosen at random |
| `lineage` | All isolates from 20% of lineages (8 lineages) |

For each predictor and design, the internal test AUC is compared with the
external AUC, and 95% percentile intervals from `amr_bootstrap()` are computed
with `group = "isolate_id"` and `group = "lineage"`.

## Running

From the package root, after installing the package dependencies:

```bash
Rscript simulation/run_study.R 200 7 200
```

The arguments are the number of replicates per scenario, CPU cores and
bootstrap replicates. Results are written to `simulation/results/`:

| File | Contents |
|---|---|
| `replicates.csv` | One row per replicate, scenario, predictor and design |
| `summary.csv` | Mean optimism, interval coverage and width per cell |
| `fig_optimism.png` | Distribution of optimism by predictor, design and scenario |
| `fig_coverage.png` | Coverage of isolate- and lineage-level intervals under the lineage-held-out design |
| `run_info.txt` | Settings, versions and run time |

## Interpreting coverage

Under the random design, the internal test set contains relatives of
training isolates, so its target is not performance on new lineages; low
coverage there reflects optimism, not only interval width. Coverage is
therefore summarised for the lineage-held-out design, where the internal test
lineages are new lineages too.

## Results (200 replicates per scenario, B = 200)

Run on 31 October 2026 (see `results/run_info.txt`).

**Optimism** (mean internal AUC minus true AUC on new lineages; 5th–95th
percentile in brackets):

| Scenario | knn / random | knn / lineage | marker / random | marker / lineage |
|---|---|---|---|---|
| Weak lineage effect (sd 0.5) | +0.20 [+0.13, +0.27] | +0.01 | +0.00 | −0.01 |
| Moderate (sd 1.5) | +0.26 [+0.18, +0.33] | +0.01 | +0.00 | −0.00 |
| Strong (sd 2.5) | +0.33 [+0.26, +0.40] | +0.00 | +0.00 | −0.01 |

A random split made a predictor with no ability to generalise (true AUC 0.50)
look like a useful one (internal AUC 0.70–0.83). The lineage-held-out design
removed the bias, at the price of noisier estimates from only 8 test lineages.

**Coverage of the true AUC** under the lineage-held-out design (nominal 95%):

| Scenario | Predictor | Isolate bootstrap | Lineage bootstrap |
|---|---|---|---|
| Weak | knn | 61% | 92% |
| Weak | marker | 87% | 90% |
| Moderate | knn | 54% | 92% |
| Moderate | marker | 64% | 91% |
| Strong | knn | 50% | 93% |
| Strong | marker | 52% | 88% |

Isolate-level intervals missed the true value up to half the time. Lineage-level
intervals were close to nominal; their remaining shortfall (88–93% rather than
95%) reflects having only 8 independent lineages, which is why
`amr_bootstrap()` warns when there are fewer than 20 groups.
