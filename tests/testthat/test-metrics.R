get <- function(m, name) m$estimate[m$metric == name]

test_that("classification metrics match a hand-calculated confusion table", {
  # threshold 0.5: TP = 3, FN = 1, FP = 1, TN = 5
  m <- amr_metrics(amr_import_predictions(toy_scores()), threshold = 0.5)
  expect_equal(get(m, "sensitivity"), 3 / 4)
  expect_equal(get(m, "specificity"), 5 / 6)
  expect_equal(get(m, "ppv"), 3 / 4)
  expect_equal(get(m, "npv"), 5 / 6)
  expect_equal(get(m, "accuracy"), 8 / 10)
  expect_equal(get(m, "balanced_accuracy"), (3 / 4 + 5 / 6) / 2)
  expect_equal(get(m, "f1"), 0.75)
  expect_equal(get(m, "mcc"), (3 * 5 - 1 * 1) / sqrt(4 * 4 * 6 * 6))
  expect_equal(m$denominator[m$metric == "sensitivity"], 4)
})

test_that("probability metrics match hand calculations", {
  d <- toy_scores()
  m <- amr_metrics(amr_import_predictions(d), threshold = 0.5)
  expect_equal(get(m, "roc_auc"), 21 / 24)
  y <- as.integer(d$truth == "R")
  expect_equal(get(m, "brier"), mean((d$probability - y)^2))
  expect_equal(get(m, "log_loss"), -mean(y * log(d$probability) + (1 - y) * log(1 - d$probability)))
  lo <- m$lower[m$metric == "roc_auc"]; hi <- m$upper[m$metric == "roc_auc"]
  expect_true(lo < 21 / 24 && hi > 21 / 24 && lo >= 0 && hi <= 1)
})

test_that("perfect and reversed predictions give AUC 1 and 0", {
  d <- toy_scores()
  d$probability <- ifelse(d$truth == "R", 0.9, 0.1)
  expect_equal(get(amr_metrics(amr_import_predictions(d)), "roc_auc"), 1)
  d$probability <- 1 - d$probability
  expect_equal(get(amr_metrics(amr_import_predictions(d)), "roc_auc"), 0)
})

test_that("scores without a threshold return threshold-free metrics only", {
  m <- amr_metrics(amr_import_predictions(toy_scores()))
  expect_false("sensitivity" %in% m$metric)
  expect_true(all(c("roc_auc", "brier", "log_loss") %in% m$metric))
  expect_error(amr_metrics(amr_import_predictions(toy_scores()), threshold = 2), "threshold")
})

test_that("intermediate policy changes the evaluable set", {
  d <- toy_scores()
  d$probability <- NULL
  d$truth[1] <- "I"
  d$prediction <- c("R", "R", "R", "S", "R", "S", "S", "S", "S", "S")
  x <- amr_import_predictions(d)
  expect_equal(get(amr_metrics(x), "n"), 9)
  expect_equal(get(amr_metrics(x, intermediate = "resistant"), "sensitivity"), 3 / 4)
  expect_equal(get(amr_metrics(x, intermediate = "susceptible"), "specificity"), 5 / 7)
})

test_that("single-class strata give NA rather than errors", {
  d <- toy_scores()
  d$truth <- "S"
  m <- amr_metrics(amr_import_predictions(d), threshold = 0.5)
  expect_true(is.na(get(m, "sensitivity")))
  expect_true(is.na(get(m, "roc_auc")))
})
