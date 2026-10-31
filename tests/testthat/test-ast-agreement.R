get <- function(m, name) m$estimate[m$metric == name]

ast_toy <- function() {
  data.frame(
    isolate_id = sprintf("iso%02d", 1:6),
    species = "E. coli",
    drug = "ciprofloxacin",
    truth = c("R", "R", "I", "S", "S", "S"),
    prediction = c("R", "S", "R", "S", "R", "S")
  )
}

test_that("categorical agreement and error rates match hand calculations", {
  a <- amr_ast_agreement(amr_import_predictions(ast_toy()))
  expect_equal(get(a, "categorical_agreement"), 3 / 6)
  expect_equal(get(a, "very_major_error"), 1 / 2)
  expect_equal(get(a, "major_error"), 1 / 3)
  expect_equal(get(a, "minor_error"), 1 / 6)
  expect_equal(a$denominator[a$metric == "very_major_error"], 2)
})

test_that("the 'all' denominator divides errors by all isolates", {
  a <- amr_ast_agreement(amr_import_predictions(ast_toy()), vme_denominator = "all")
  expect_equal(get(a, "very_major_error"), 1 / 6)
  expect_equal(get(a, "major_error"), 1 / 6)
})

test_that("scores need a threshold for categories", {
  x <- amr_import_predictions(toy_scores())
  expect_error(amr_ast_agreement(x), "threshold")
  a <- amr_ast_agreement(x, threshold = 0.5)
  expect_equal(get(a, "very_major_error"), 1 / 4)
  expect_equal(get(a, "major_error"), 1 / 6)
})

test_that("essential agreement and bias handle on- and off-scale MICs", {
  d <- ast_toy()
  d$truth <- NULL
  d$prediction <- NULL
  d$truth_mic <- c(">8", "16", "8", "4", "0.5", "<=0.125")
  d$predicted_mic <- c("16", "8", "8", "1", "0.5", "<=0.125")
  a <- amr_ast_agreement(amr_import_predictions(d))
  expect_equal(get(a, "essential_agreement"), 5 / 6)
  expect_equal(get(a, "bias_log2"), mean(c(-1, 0, -2, 0)))
  d$truth_mic <- NULL
  expect_error(amr_ast_agreement(amr_import_predictions(d)), "truth_mic")
})
