two_drugs <- function() {
  d1 <- toy_scores()
  d2 <- toy_scores()
  d2$drug <- "ciprofloxacin"
  d2$probability <- rev(d2$probability)
  rbind(d1, d2)
}

test_that("evaluation is stratified by species and drug", {
  e <- amr_species_eval(amr_import_predictions(two_drugs()), threshold = 0.5)
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
  e <- amr_species_eval(amr_import_predictions(d), threshold = 0.5)
  expect_equal(e$estimate[e$metric == "n"], 8)
  e_all <- amr_species_eval(amr_import_predictions(d), threshold = 0.5, split = NULL)
  expect_equal(e_all$estimate[e_all$metric == "n"], 10)
})

test_that("invalid data stops evaluation", {
  d <- toy_scores()
  d$probability[1] <- 3
  expect_error(suppressMessages(capture.output(amr_species_eval(amr_import_predictions(d), threshold = 0.5))),
               "Validation failed")
})

test_that("missing threshold gives a message and threshold-free metrics", {
  expect_message(e <- amr_species_eval(amr_import_predictions(toy_scores())), "threshold")
  expect_false("sensitivity" %in% e$metric)
  expect_false("ast" %in% e$family)
})

test_that("MIC predictions are evaluated with AST metrics only", {
  d <- data.frame(isolate_id = c("a", "b", "c"), species = "E. coli", drug = "ciprofloxacin",
                  truth_mic = c("1", "2", "0.5"), predicted_mic = c("1", "8", "0.5"))
  e <- amr_species_eval(amr_import_predictions(d))
  expect_equal(unique(e$family), "ast")
  expect_equal(e$estimate[e$metric == "essential_agreement"], 2 / 3)
})
