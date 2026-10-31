split_toy <- function() {
  data.frame(
    isolate_id = c("a", "b", "c", "d", "e", "f"),
    species = "K. pneumoniae",
    drug = "meropenem",
    truth = c("R", "S", "R", "S", "R", "S"),
    probability = c(0.9, 0.1, 0.8, 0.3, 0.7, 0.2),
    split = c("train", "train", "train", "test", "test", "test"),
    lineage = c("ST258", "ST11", "ST258", "ST15", "ST258", "ST307"),
    collection_date = c("2019-01-01", "2019-06-01", "2020-01-01", "2021-01-01", "2019-03-01", "2022-01-01"),
    country = c("CA", "CA", "US", "FR", "CA", "DE")
  )
}

test_that("shared lineages are counted and contradict a group design", {
  s <- amr_validate_splits(amr_import_predictions(split_toy()), claimed_design = "group")
  expect_s3_class(s, "amr_split_check")
  row <- s$summary[s$summary$check == "lineage", ]
  expect_equal(row$n_shared, 1)
  expect_equal(row$prop_test_isolates_affected, 1 / 3)
  expect_false(s$ok)
  expect_true("group" %in% s$violated)
  expect_true(s$consistent_designs$isolate)
  expect_output(print(s), "PROBLEMS FOUND")
})

test_that("temporal and geographic overlap are detected", {
  s <- amr_validate_splits(amr_import_predictions(split_toy()))
  expect_equal(s$summary$n_shared[s$summary$check == "collection_date"], 1)
  expect_false(s$consistent_designs$temporal)
  expect_equal(s$summary$n_shared[s$summary$check == "country"], 1)
  expect_false(s$consistent_designs$geographic)
  expect_true(s$ok)
})

test_that("a clean lineage-held-out split passes its claim", {
  d <- split_toy()
  d$lineage[5] <- "ST147"
  s <- amr_validate_splits(amr_import_predictions(d), claimed_design = "group")
  expect_true(s$ok)
  expect_true(s$consistent_designs$lineage)
})

test_that("a training manifest replaces train rows and can be model-specific", {
  d <- split_toy()[4:6, ]
  d$split <- NULL
  d$model_id <- "m1"
  manifest <- data.frame(isolate_id = c("x1", "x2", "d"), model_id = "m1",
                         lineage = c("ST11", "ST258", "ST15"))
  s <- amr_validate_splits(amr_import_predictions(d), train_manifest = manifest)
  expect_false(s$ok)
  expect_equal(s$summary$n_shared[s$summary$check == "isolate"], 1)
  expect_equal(s$summary$n_shared[s$summary$check == "lineage"], 2)
})

test_that("checks that cannot be assessed are reported, and claims on them fail", {
  d <- split_toy()[4:6, ]
  d$split <- NULL
  manifest <- data.frame(isolate_id = c("x1", "x2"))
  s <- amr_validate_splits(amr_import_predictions(d), train_manifest = manifest, claimed_design = "group")
  expect_true(is.na(s$summary$n_shared[s$summary$check == "lineage"]))
  expect_false(s$ok)
  expect_match(s$violated, "not assessable")
})

test_that("missing training information is an error", {
  d <- split_toy()
  d$split <- "test"
  expect_error(amr_validate_splits(amr_import_predictions(d)), "No training isolates")
  expect_error(amr_validate_splits(amr_import_predictions(split_toy()), group = "nope"), "not found")
  expect_error(amr_validate_splits(amr_import_predictions(split_toy()), train_manifest = data.frame(id = 1)),
               "isolate_id")
})
