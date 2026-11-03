# AMRvalid 0.2.0.9000

* `amr_calibration()`: calibration-in-the-large (O/E ratio and intercept),
  calibration slope, expected and maximum calibration error, binned and
  smoothed reliability curves, optional grouping column, and a `plot()`
  method drawing reliability diagrams.
* `amr_distance_eval()` now reports the calibration intercept and slope in
  each distance band for score predictions.
* `target_prevalence` in `amr_metrics()` and `amr_species_eval()` adds PPV
  and NPV at the prevalence of the setting where the predictor will be used;
  `amr_ppv_npv()` does the conversion directly.
* `weight_by` in `amr_metrics()`, `amr_ast_agreement()` and
  `amr_species_eval()` gives cluster-weighted metrics in which each lineage
  or genomic cluster counts once, so outbreak clones do not dominate.
* `amr_provenance()`: what was evaluated and how, each item labelled as
  checked by AMRvalid, supplied by the analyst or not assessable, and mapped
  to the TRIPOD+AI checklist (version 11 January 2024).
* `amr_report()`: Markdown or HTML report with data summary, validation and
  split checks, performance, calibration, distance-to-training results,
  group-aware intervals, detected limitations and TRIPOD+AI coverage.
* Simulation study: each bootstrap is now seeded, so results no longer depend
  on the number of CPU cores. Re-run results: isolate-level intervals covered
  the true AUC in 49-87% of replicates, lineage-level intervals in 86-93%;
  optimism is unchanged.

# AMRvalid 0.2.0

Population-aware evaluation.

* `amr_validate_splits()`: isolate, lineage, cluster, sequence-type, temporal
  and geographic overlap between training and test isolates, with checks
  against a claimed validation design.
* `amr_distance_eval()` and `amr_nearest_train_distance()`: performance by
  distance to the nearest training isolate, with a trend summary of how
  performance decays as test isolates become less similar to training data.
* `amr_bootstrap()`: percentile intervals from resampling isolates,
  lineages or genomic clusters, with warnings when there are few groups and
  when intervals have zero width.
* `amr_paired_compare()`: paired bootstrap differences between two models on
  matched isolates, plus an exact McNemar test.
* `amr_simulate()`: synthetic lineages, genomic distances and predictions from
  a lineage-memorising and a mechanism-based predictor under random and
  lineage-held-out splits, with an external set of unseen lineages.
* New vignette, "Why population-aware evaluation matters".
* Simulation study (`simulation/`, not part of the installed package): random
  splits overstated a lineage-memorising predictor by +0.20 to +0.33 AUC, and
  isolate-level bootstrap intervals covered the true AUC in only 50-87% of
  replicates, against 88-93% for lineage-level intervals.
* Subsets and reshaped results now print as plain data frames.
* Figures use a colourblind-safe palette.
* Citation information via `citation("AMRvalid")`.

# AMRvalid 0.1.0

Core evaluation.

* Data contract (`amr_contract()`, JSON Schema), `amr_import_predictions()`,
  `validate_amr_data()` and `amr_data_summary()`. Species and antimicrobial
  names must match the AMR package exactly by default, because its fuzzy
  matching can silently map a typo to a different antimicrobial.
* `amr_metrics()`: sensitivity, specificity, PPV, NPV, accuracy, balanced
  accuracy, F1, MCC, ROC-AUC with DeLong interval, Brier score and log loss.
  No default threshold; I/SDD handled by an explicit `intermediate` policy.
* `amr_ast_agreement()`: categorical agreement, very major, major and minor
  errors (class or all-isolate denominators), essential agreement and bias
  for MIC predictions.
* `amr_species_eval()`: per model, dataset, species and antimicrobial
  evaluation; test split only by default; refuses to run on invalid data.
* `amr_prop_ci()`: Wilson and Clopper-Pearson intervals.
