# AMRvalid 0.0.0.9000

* Package scaffold, licence and continuous-integration checks.
* Data contract (`amr_contract()`, JSON Schema), `amr_import_predictions()`,
  `validate_amr_data()` and `amr_data_summary()`.
* `amr_metrics()`: sensitivity, specificity, PPV, NPV, accuracy, balanced
  accuracy, F1, MCC, ROC-AUC with DeLong interval, Brier score and log loss.
  No default threshold; I/SDD handled by an explicit `intermediate` policy.
* `amr_ast_agreement()`: categorical agreement, very major, major and minor
  errors (class or all-isolate denominators), essential agreement and bias
  for MIC predictions.
* `amr_species_eval()`: per model, dataset, species and antimicrobial
  evaluation; test split only by default; refuses to run on invalid data.
* `amr_prop_ci()`: Wilson and Clopper-Pearson intervals.
* `amr_validate_splits()`: isolate, lineage, cluster, sequence-type, temporal
  and geographic overlap between training and test isolates, with checks
  against a claimed validation design.
* `amr_distance_eval()` and `amr_nearest_train_distance()`: performance by
  distance to the nearest training isolate, with a trend summary of how
  performance decays as test isolates become less similar to training data.
* `amr_bootstrap()`: percentile intervals from resampling isolates,
  lineages or genomic clusters, with warnings when there are few groups.
* `amr_paired_compare()`: paired bootstrap differences between two models on
  matched isolates, plus an exact McNemar test.
* `amr_simulate()`: synthetic lineages, genomic distances and predictions from
  a lineage-memorising and a mechanism-based predictor under random and
  lineage-held-out splits, with an external set of unseen lineages.
