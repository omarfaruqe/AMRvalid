test_that("example data loads and follows the data contract", {
  d <- amr_example_data()
  expect_s3_class(d, "data.frame")
  expect_true(all(c("isolate_id", "species", "drug", "truth", "prediction", "model_id",
                    "genomic_cluster", "collection_year", "country") %in% names(d)))
  expect_setequal(unique(d$drug), c("meropenem", "ciprofloxacin"))
  expect_equal(length(unique(d$model_id)), 4)
  expect_true(all(d$truth %in% c("S", "I", "R")))
  expect_true(all(d$prediction %in% c("S", "R")))
  expect_true(is.character(d$genomic_cluster))
  x <- amr_import_predictions(d)
  expect_equal(unique(x$prediction_type), "binary")
  v <- amr_validate(x)
  expect_true(v$ok)
})

test_that("each predictor covers the same isolates for its antimicrobial", {
  d <- amr_example_data()
  for (drug in unique(d$drug)) {
    ids <- split(d$isolate_id[d$drug == drug], d$model_id[d$drug == drug])
    expect_length(ids, 2)
    expect_setequal(ids[[1]], ids[[2]])
  }
})
