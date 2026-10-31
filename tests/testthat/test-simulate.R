test_that("simulation returns importable predictions with the expected structure", {
  sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10, external_lineages = 10, seed = 1)
  expect_named(sim, c("predictions", "train_manifest", "isolates", "distances", "parameters"))
  expect_equal(nrow(sim$isolates), 300)
  expect_setequal(unique(sim$predictions$model_id),
                  c("knn/random", "knn/lineage", "marker/random", "marker/lineage"))
  expect_equal(dim(sim$distances), c(300, 300))
  x <- amr_import_predictions(sim$predictions)
  expect_true(validate_amr_data(x)$ok)
  expect_true(all(x$probability >= 0 & x$probability <= 1))
})

test_that("designs have the advertised train/test separation", {
  sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10, external_lineages = 0, seed = 2)
  x <- amr_import_predictions(sim$predictions)
  s <- amr_validate_splits(x, train_manifest = sim$train_manifest, group = "lineage")
  cons <- s$consistent_designs
  expect_true(all(cons$isolate))
  expect_true(cons$lineage[cons$model_id == "knn/lineage"])
  expect_false(cons$lineage[cons$model_id == "knn/random"])
  expect_equal(sum(sim$predictions$model_id == "knn/random"), 40)
})

test_that("simulation is reproducible and restores the random state", {
  set.seed(5); before <- .Random.seed
  a <- amr_simulate(n_lineages = 10, isolates_per_lineage = 5, external_lineages = 5, seed = 3)
  expect_identical(.Random.seed, before)
  b <- amr_simulate(n_lineages = 10, isolates_per_lineage = 5, external_lineages = 5, seed = 3)
  expect_identical(a$predictions, b$predictions)
})

test_that("random splits overstate the lineage-memorising predictor", {
  sim <- amr_simulate(n_lineages = 40, isolates_per_lineage = 20, external_lineages = 60,
                      return_distances = FALSE, seed = 11)
  x <- amr_import_predictions(sim$predictions)
  e <- suppressMessages(amr_species_eval(x, validate = FALSE))
  auc <- function(model, dataset) e$estimate[e$metric == "roc_auc" & e$model_id == model &
                                               e$dataset_id == dataset]
  optimism_knn <- auc("knn/random", "internal") - auc("knn/random", "external")
  optimism_marker <- auc("marker/random", "internal") - auc("marker/random", "external")
  expect_gt(optimism_knn, 0.05)
  expect_gt(optimism_knn, optimism_marker)
  expect_null(sim$distances)
})
