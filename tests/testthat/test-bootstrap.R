# 20 lineages of 5 isolates; whole lineages are either predicted correctly
# or incorrectly, so isolates within a lineage are strongly dependent.
clustered_toy <- function() {
  lin <- rep(sprintf("L%02d", 1:20), each = 5)
  truth <- rep(c("R", "S", "R", "S", "S"), 20)
  wrong_lineage <- lin %in% sprintf("L%02d", c(2, 5, 9, 14, 17))
  prob <- ifelse(truth == "R", 0.8, 0.2)
  prob[wrong_lineage] <- 1 - prob[wrong_lineage]
  data.frame(isolate_id = sprintf("iso%03d", 1:100), species = "K. pneumoniae", drug = "meropenem",
             truth = truth, probability = prob, lineage = lin)
}

test_that("fast metrics match amr_metrics", {
  x <- amr_import_predictions(toy_scores())
  m <- amr_metrics(x, threshold = 0.5)
  y <- as.integer(toy_scores()$truth == "R")
  f <- AMRvalid:::fast_metrics(y, as.integer(toy_scores()$probability >= 0.5), toy_scores()$probability)
  for (nm in c("sensitivity", "specificity", "ppv", "npv", "accuracy", "mcc", "f1", "roc_auc", "brier", "log_loss")) {
    expect_equal(unname(f[nm]), m$estimate[m$metric == nm], info = nm)
  }
  expect_equal(unname(f["very_major_error"]), 1 - unname(f["sensitivity"]))
})

test_that("bootstrap is reproducible and leaves the global RNG untouched", {
  x <- amr_import_predictions(clustered_toy())
  set.seed(99); before <- .Random.seed
  a <- amr_bootstrap(x, B = 100, threshold = 0.5, seed = 7)
  expect_identical(.Random.seed, before)
  b <- amr_bootstrap(x, B = 100, threshold = 0.5, seed = 7)
  expect_equal(a$lower, b$lower)
  expect_s3_class(a, "amr_bootstrap")
  expect_equal(a$estimate[a$metric == "accuracy"], 0.75)
})

test_that("cluster bootstrap gives wider intervals than isolate bootstrap", {
  x <- amr_import_predictions(clustered_toy())
  iso <- amr_bootstrap(x, resample_by = "isolate_id", B = 1000, threshold = 0.5, seed = 1)
  lin <- amr_bootstrap(x, resample_by = "lineage", B = 1000, threshold = 0.5, seed = 1)
  w <- function(r) diff(unlist(r[r$metric == "accuracy", c("lower", "upper")]))
  expect_gt(w(lin), 1.5 * w(iso))
  expect_equal(unique(lin$n_groups), 20)
  expect_equal(unique(iso$n_groups), 100)
})

test_that("few groups are flagged and missing groups become singletons", {
  d <- clustered_toy()
  d$lineage <- rep(c("A", "B", "C", "D"), 25)
  d$lineage[1:4] <- NA
  r <- amr_bootstrap(amr_import_predictions(d), resample_by = "lineage", B = 50, threshold = 0.5, seed = 1)
  expect_true(all(r$few_groups))
  expect_equal(unique(r$n_groups), 8)
  expect_output(print(r), "Fewer than 20")
})

test_that("zero-width intervals from perfect predictions are flagged", {
  d <- clustered_toy()
  d$probability <- ifelse(d$truth == "R", 0.9, 0.1)
  r <- amr_bootstrap(amr_import_predictions(d), B = 50, threshold = 0.5, seed = 1)
  expect_true(r$zero_width[r$metric == "very_major_error"])
  expect_true(r$zero_width[r$metric == "brier"])
  expect_output(print(r), "zero width")
  ok <- amr_bootstrap(amr_import_predictions(clustered_toy()), B = 50, threshold = 0.5, seed = 1)
  expect_false(any(ok$zero_width))
})

test_that("bootstrap argument checks", {
  x <- amr_import_predictions(clustered_toy())
  expect_error(amr_bootstrap(x, resample_by = "nope"), "not found")
  expect_warning(old <- amr_bootstrap(x, group = "lineage", B = 20, seed = 1), "deprecated")
  expect_equal(old$lower, amr_bootstrap(x, resample_by = "lineage", B = 20, seed = 1)$lower)
  expect_error(amr_bootstrap(x, B = 3), "at least 10")
  r <- amr_bootstrap(x, B = 20, seed = 1)
  expect_false("sensitivity" %in% r$metric)
  expect_true("roc_auc" %in% r$metric)
})

paired_toy <- function() {
  d <- clustered_toy()
  a <- d; a$model_id <- "old"
  b <- d; b$model_id <- "new"
  b$probability <- ifelse(b$truth == "R", 0.8, 0.2)
  rbind(a, b)
}

test_that("identical models give zero differences", {
  d <- clustered_toy()
  a <- d; a$model_id <- "A"
  b <- d; b$model_id <- "B"
  r <- amr_paired_compare(amr_import_predictions(rbind(a, b)), B = 50, threshold = 0.5, seed = 1)
  expect_true(all(r$difference == 0, na.rm = TRUE))
  expect_true(all(r$lower == 0 & r$upper == 0, na.rm = TRUE))
})

test_that("a better model shows a positive paired difference and McNemar p", {
  r <- amr_paired_compare(amr_import_predictions(paired_toy()), models = c("old", "new"),
                          resample_by = "lineage", B = 500, threshold = 0.5, seed = 1)
  expect_s3_class(r, "amr_paired_compare")
  acc <- r[r$metric == "accuracy", ]
  expect_equal(acc$difference, 0.25)
  expect_gt(acc$lower, 0)
  expect_equal(acc$mcnemar_p, stats::binom.test(25, 25)$p.value)
  expect_equal(acc$n_matched, 100)
  expect_output(print(r), "new - old")
})

test_that("unmatched isolates are excluded and thresholds can differ by model", {
  d <- paired_toy()
  d <- d[-c(1, 2), ]
  r <- amr_paired_compare(amr_import_predictions(d), models = c("old", "new"),
                          threshold = c(old = 0.5, new = 0.9), B = 50, seed = 1)
  expect_equal(unique(r$n_matched), 98)
  expect_equal(unique(r$n_unmatched), 2)
  expect_equal(r$model_b[r$metric == "sensitivity"], 0)
  expect_error(amr_paired_compare(amr_import_predictions(d), threshold = c(0.5, 0.6), B = 20),
               "named")
})

test_that("paired comparison argument checks", {
  x <- amr_import_predictions(clustered_toy())
  expect_error(amr_paired_compare(x, B = 20), "Exactly two")
  expect_error(amr_paired_compare(amr_import_predictions(paired_toy()), models = c("old", "x"), B = 20),
               "not found")
})
