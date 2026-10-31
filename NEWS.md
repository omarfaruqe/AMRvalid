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
