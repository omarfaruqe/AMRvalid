cal_toy <- function(n = 5000, transform = identity, seed = 1) {
  set.seed(seed)
  p_true <- stats::runif(n, 0.05, 0.95)
  y <- stats::rbinom(n, 1, p_true)
  data.frame(isolate_id = sprintf("iso%05d", seq_len(n)), species = "E. coli", drug = "ciprofloxacin",
             truth = ifelse(y == 1, "R", "S"), probability = transform(p_true))
}
s_of <- function(d, ...) amr_calibration(amr_import_predictions(d), ...)$summary

test_that("well-calibrated predictions give slope 1, intercept 0 and O/E 1", {
  s <- s_of(cal_toy())
  expect_true(s$slope_lower < 1 && s$slope_upper > 1)
  expect_true(s$intercept_lower < 0 && s$intercept_upper > 0)
  expect_equal(s$oe_ratio, 1, tolerance = 0.05)
  expect_lt(s$ece, 0.03)
})

test_that("over-confident predictions give a slope below 1", {
  s <- s_of(cal_toy(transform = function(p) stats::plogis(2 * stats::qlogis(p))))
  expect_equal(s$slope, 0.5, tolerance = 0.08)
  expect_lt(s$slope_upper, 1)
})

test_that("predictions that are too low give a positive intercept", {
  s <- s_of(cal_toy(transform = function(p) stats::plogis(stats::qlogis(p) - 1)))
  expect_equal(s$intercept, 1, tolerance = 0.1)
  expect_gt(s$oe_ratio, 1)
})

test_that("ECE and MCE match a hand calculation", {
  d <- data.frame(isolate_id = c("a", "b", "c", "d"), species = "E. coli", drug = "ciprofloxacin",
                  truth = c("S", "R", "R", "R"), probability = c(0.1, 0.1, 0.9, 0.9))
  cal <- amr_calibration(amr_import_predictions(d), bins = 2, min_events = 1)
  expect_equal(cal$summary$ece, (2 * 0.4 + 2 * 0.1) / 4)
  expect_equal(cal$summary$mce, 0.4)
  expect_equal(nrow(cal$bins), 2)
  expect_equal(cal$bins$observed, c(0.5, 1))
})

test_that("equal-count binning and grouping by a column work", {
  d <- cal_toy(n = 600)
  d$lineage <- rep(c("L1", "L2", "L3"), 200)
  cal <- amr_calibration(amr_import_predictions(d), binning = "equal_count", by = "lineage")
  expect_equal(nrow(cal$summary), 3)
  expect_setequal(cal$summary$lineage, c("L1", "L2", "L3"))
  expect_true(all(table(cal$bins$lineage) == 10))
  expect_output(print(cal), "Ideal")
})

test_that("degenerate inputs are handled and reported", {
  d <- cal_toy(n = 50)
  d$truth <- "S"
  s <- s_of(d)
  expect_true(is.na(s$slope))
  expect_true(s$sparse)
  expect_match(s$note, "one phenotype class")
  d2 <- cal_toy(n = 50)
  d2$probability <- 0.3
  expect_match(s_of(d2)$note, "equal")
  b <- cal_toy(n = 20)
  b$probability <- NULL
  b$prediction <- b$truth
  expect_error(amr_calibration(amr_import_predictions(b)), "probabilities")
  expect_error(amr_calibration(amr_import_predictions(cal_toy(50)), bins = 1), "bins")
  expect_error(amr_calibration(amr_import_predictions(cal_toy(50)), by = "nope"), "not found")
})

test_that("plot draws without error and returns the object", {
  cal <- amr_calibration(amr_import_predictions(cal_toy(n = 500)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_identical(plot(cal), cal)
})

test_that("distance evaluation reports calibration per band", {
  sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10, external_lineages = 20, seed = 3)
  x <- amr_import_predictions(sim$predictions)
  r <- amr_distance_eval(x[x$model_id == "marker/random", ], sim$distances,
                         train_manifest = sim$train_manifest, bands = 2)
  cal <- r$by_band[r$by_band$family == "calibration", ]
  expect_setequal(unique(cal$metric), c("calibration_intercept", "calibration_slope"))
  combos <- unique(r$by_band[r$by_band$family == "ml", c("dataset_id", "band")])
  expect_equal(nrow(cal), 2 * nrow(combos))
  expect_output(print(r), "calibration_slope")
})
