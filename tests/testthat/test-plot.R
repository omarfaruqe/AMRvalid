plot_ok <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(expr)
}

test_that("plot methods draw without error and return their input", {
  x <- amr_import_predictions(amr_example_data())
  e <- amr_evaluate(x)
  expect_identical(plot_ok(plot(e)), e)
  expect_identical(plot_ok(plot(e, metrics = c("very_major_error", "major_error", "minor_error"))), e)
  eb <- amr_evaluate(x, by = "country")
  expect_identical(plot_ok(plot(eb[eb$country %in% c("USA", "Brazil"), ], metrics = "sensitivity")),
                   eb[eb$country %in% c("USA", "Brazil"), ])
  b <- amr_bootstrap(x, resample_by = "genomic_cluster", B = 20, seed = 1)
  expect_identical(plot_ok(plot(b)), b)
  mero <- x[x$drug == "MEM", ]
  p <- amr_paired_compare(mero, models = c("carbapenemase genes", "carbapenemase or porin + ESBL/AmpC"),
                          B = 20, seed = 1)
  expect_identical(plot_ok(plot(p, metrics = c("sensitivity", "specificity"))), p)
})

test_that("distance and score plots draw, and bad metrics are reported", {
  sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10, external_lineages = 20, seed = 1)
  x <- amr_import_predictions(sim$predictions)
  d <- amr_distance_eval(x[x$model_id %in% c("knn/random", "marker/random"), ], sim$distances,
                         train_manifest = sim$train_manifest)
  expect_identical(plot_ok(plot(d)), d)
  expect_identical(plot_ok(plot(d, metric = "brier")), d)
  expect_error(plot_ok(plot(d, metric = "nope")), "not in the results")
  e <- suppressMessages(amr_evaluate(x))
  expect_identical(plot_ok(plot(e)), e)
  expect_error(plot_ok(plot(e, metrics = "sensitivity")), "not in the results")
})

test_that("colours follow sorted model identifiers and are capped", {
  cols <- AMRvalid:::model_colours(c("b", "a", "b"))
  expect_equal(names(cols), c("a", "b"))
  expect_equal(unname(cols), AMRvalid:::amr_palette[1:2])
  expect_error(AMRvalid:::model_colours(letters[1:9]), "at most 8")
})
