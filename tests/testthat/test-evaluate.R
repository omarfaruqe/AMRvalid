two_drugs <- function() {
  d1 <- toy_scores()
  d2 <- toy_scores()
  d2$drug <- "ciprofloxacin"
  d2$probability <- rev(d2$probability)
  rbind(d1, d2)
}

test_that("evaluation is stratified by species and drug", {
  e <- amr_evaluate(amr_import_predictions(two_drugs()), threshold = 0.5)
  expect_s3_class(e, "amr_eval")
  expect_setequal(unique(e$drug_name), c("Meropenem", "Ciprofloxacin"))
  sens <- e$estimate[e$metric == "sensitivity" & e$drug == "MEM"]
  expect_equal(sens, 3 / 4)
  expect_true(all(c("ml", "ast") %in% e$family))
  expect_output(print(e), "amr_eval")
})

test_that("only the test split is evaluated by default", {
  d <- toy_scores()
  d$split <- c(rep("test", 8), "validation", "validation")
  e <- amr_evaluate(amr_import_predictions(d), threshold = 0.5)
  expect_equal(e$estimate[e$metric == "n"], 8)
  e_all <- amr_evaluate(amr_import_predictions(d), threshold = 0.5, split = NULL)
  expect_equal(e_all$estimate[e_all$metric == "n"], 10)
})

test_that("invalid data stops evaluation", {
  d <- toy_scores()
  d$probability[1] <- 3
  expect_error(amr_evaluate(amr_import_predictions(d), threshold = 0.5),
               "Validation failed:\n\\* probability_range")
  expect_silent(try(amr_evaluate(amr_import_predictions(d), threshold = 0.5), silent = TRUE))
})

test_that("missing threshold gives a message and threshold-free metrics", {
  expect_message(e <- amr_evaluate(amr_import_predictions(toy_scores())), "threshold")
  expect_false("sensitivity" %in% e$metric)
  expect_false("ast" %in% e$family)
})

test_that("MIC predictions are evaluated with AST metrics only", {
  d <- data.frame(isolate_id = c("a", "b", "c"), species = "E. coli", drug = "ciprofloxacin",
                  truth_mic = c("1", "2", "0.5"), predicted_mic = c("1", "8", "0.5"))
  e <- amr_evaluate(amr_import_predictions(d))
  expect_equal(unique(e$family), "ast")
  expect_equal(e$estimate[e$metric == "essential_agreement"], 2 / 3)
})

test_that("subsets and reshaped results print as plain data frames", {
  e <- amr_evaluate(amr_import_predictions(toy_scores()), threshold = 0.5)
  sub <- e[e$metric == "accuracy", c("model_id", "estimate")]
  expect_output(print(sub), "estimate")
  x <- amr_import_predictions(toy_scores())
  expect_output(print(x[, c("isolate_id", "truth_sir")]), "isolate_id")
})

test_that("the old name amr_species_eval() still works with a deprecation warning", {
  x <- amr_import_predictions(toy_scores())
  expect_warning(e <- amr_species_eval(x, threshold = 0.5), "amr_evaluate")
  expect_equal(as.data.frame(e), as.data.frame(amr_evaluate(x, threshold = 0.5)))
})

test_that("by = evaluates subgroups within each stratum", {
  d <- toy_scores()
  d$country <- c("CA", "CA", "US", "US", "CA", "CA", "US", "US", NA, "CA")
  x <- amr_import_predictions(d)
  e <- amr_evaluate(x, threshold = 0.5, by = "country")
  expect_setequal(unique(e$country), c("CA", "US", "(missing)"))
  n <- e$estimate[e$metric == "n"]
  names(n) <- e$country[e$metric == "n"]
  expect_equal(unname(n[c("CA", "US", "(missing)")]), c(5, 4, 1))
  # CA: isolates 1, 2, 5, 6, 10 -> truth R, R, S, S, S; predictions R, R, R, S, S
  ca <- e[e$country == "CA", ]
  expect_equal(ca$estimate[ca$metric == "sensitivity"], 1)
  expect_equal(ca$estimate[ca$metric == "specificity"], 2 / 3)
  expect_output(print(e), "by country")
  expect_error(amr_evaluate(x, by = "nope"), "not found")
  expect_error(amr_evaluate(x, by = "drug"), "already strata")
})

test_that("by = accepts several columns and combines with weighting", {
  d <- toy_scores()
  d$country <- rep(c("CA", "US"), 5)
  d$period <- rep(c("early", "late"), each = 5)
  d$lineage <- rep(c("L1", "L2"), each = 5)
  e <- amr_evaluate(amr_import_predictions(d), threshold = 0.5, by = c("country", "period"),
                    weight_by = "lineage")
  expect_equal(nrow(unique(e[, c("country", "period")])), 4)
  expect_equal(sum(e$estimate[e$metric == "n"]), 10)
})

test_that("factor by-columns keep their level order", {
  d <- toy_scores()
  d$period <- factor(rep(c("late", "early"), each = 5), levels = c("early", "late"))
  e <- amr_evaluate(amr_import_predictions(d), threshold = 0.5, by = "period")
  expect_equal(unique(e$period), c("early", "late"))
})
