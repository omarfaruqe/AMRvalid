test_that("Wilson interval matches prop.test without continuity correction", {
  ref <- suppressWarnings(stats::prop.test(3, 4, correct = FALSE)$conf.int)
  ci <- amr_prop_ci(3, 4)
  expect_equal(ci$estimate, 0.75)
  expect_equal(c(ci$lower, ci$upper), as.numeric(ref), tolerance = 1e-10)
})

test_that("Clopper-Pearson interval matches binom.test", {
  ref <- stats::binom.test(1, 12)$conf.int
  ci <- amr_prop_ci(1, 12, method = "clopper-pearson")
  expect_equal(c(ci$lower, ci$upper), as.numeric(ref), tolerance = 1e-10)
  edge <- amr_prop_ci(c(0, 5), c(5, 5), method = "clopper-pearson")
  expect_equal(edge$lower[1], 0)
  expect_equal(edge$upper[2], 1)
})

test_that("zero denominators give NA and bad input errors", {
  ci <- amr_prop_ci(0, 0)
  expect_true(is.na(ci$estimate) && is.na(ci$lower) && is.na(ci$upper))
  expect_error(amr_prop_ci(5, 4), "between 0")
  expect_error(amr_prop_ci(1, 4, conf = 1.2), "conf")
})
