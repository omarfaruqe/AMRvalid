test_that("clean data passes with only warnings or notes", {
  v <- validate_amr_data(amr_import_predictions(toy_scores()), min_n = 5, min_resistant = 2)
  expect_s3_class(v, "amr_validation")
  expect_true(v$ok)
  expect_false(any(v$issues$level == "error"))
})

test_that("duplicates, bad probabilities and conflicting truth are errors", {
  d <- toy_scores()
  d <- rbind(d, d[1, ])
  d$probability[2] <- 1.5
  d$truth[11] <- "S"
  v <- validate_amr_data(amr_import_predictions(d))
  expect_false(v$ok)
  expect_true(all(c("duplicate_key", "probability_range", "conflicting_truth") %in% v$issues$check))
})

test_that("train/test overlap within a model is detected", {
  d <- toy_scores()
  d$split <- "test"
  d2 <- d[1:2, ]
  d2$split <- "train"
  v <- validate_amr_data(amr_import_predictions(rbind(d, d2)))
  expect_true("split_overlap" %in% v$issues$check)
  expect_false(v$ok)
})

test_that("unknown split labels are errors", {
  d <- toy_scores()
  d$split <- "holdout"
  expect_true("split_label" %in% validate_amr_data(amr_import_predictions(d))$issues$check)
})

test_that("single-class and few-resistant strata are warned", {
  d <- toy_scores()
  d$truth <- "S"
  v <- validate_amr_data(amr_import_predictions(d))
  expect_true("single_class" %in% v$issues$check)
  d <- toy_scores()
  d$truth <- c("R", rep("S", 9))
  expect_true("few_resistant" %in% validate_amr_data(amr_import_predictions(d))$issues$check)
})

test_that("intermediate results and mixed breakpoints are flagged", {
  d <- toy_scores()
  d$truth[1] <- "I"
  d$breakpoint_standard <- c(rep("EUCAST 2024", 5), rep("EUCAST 2025", 5))
  v <- validate_amr_data(amr_import_predictions(d))
  expect_true(all(c("intermediate", "mixed_breakpoints") %in% v$issues$check))
})

test_that("missing required columns stop early", {
  d <- toy_scores()
  d$truth <- NULL
  v <- validate_amr_data(amr_import_predictions(d))
  expect_false(v$ok)
  expect_equal(v$issues$check[1], "required_column")
})

test_that("data summary counts phenotypes per stratum", {
  d <- toy_scores()
  d$truth[1] <- "I"
  s <- amr_data_summary(amr_import_predictions(d))
  expect_equal(nrow(s), 1)
  expect_equal(s$n_R, 3)
  expect_equal(s$n_I_SDD, 1)
  expect_equal(s$prevalence_R, 3 / 9)
  expect_equal(s$drug_name, "Meropenem")
})
