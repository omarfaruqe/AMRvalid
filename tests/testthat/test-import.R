test_that("import coerces identifiers to AMR classes and detects score type", {
  x <- amr_import_predictions(toy_scores())
  expect_s3_class(x, "amr_data")
  expect_s3_class(x$species, "mo")
  expect_s3_class(x$drug, "ab")
  expect_s3_class(x$truth_sir, "sir")
  expect_equal(as.character(x$drug[1]), "MEM")
  expect_equal(unique(x$prediction_type), "score")
  expect_equal(unique(x$model_id), "model_1")
})

test_that("mapping renames columns", {
  d <- toy_scores()
  names(d)[names(d) == "truth"] <- "phenotype"
  names(d)[names(d) == "isolate_id"] <- "Sample"
  x <- amr_import_predictions(d, mapping = c(isolate_id = "Sample", truth_sir = "phenotype"))
  expect_true(all(c("isolate_id", "truth_sir") %in% names(x)))
  expect_error(amr_import_predictions(d, mapping = c(isolate_id = "nope")), "not found")
})

test_that("0/1 and logical predictions become S/R", {
  d <- toy_scores()
  d$probability <- NULL
  d$prediction <- c(1, 1, 1, 0, 1, 0, 0, 0, 0, 0)
  x <- amr_import_predictions(d)
  expect_equal(as.character(x$prediction[1:4]), c("R", "R", "R", "S"))
  expect_equal(unique(x$prediction_type), "binary")
  d$prediction <- d$prediction == 1
  expect_equal(as.character(amr_import_predictions(d)$prediction[4]), "S")
  d$prediction <- c(2, 1, 1, 0, 1, 0, 0, 0, 0, 0)
  expect_error(amr_import_predictions(d), "only 0, 1")
})

test_that("S/I/R predictions are detected", {
  d <- toy_scores()
  d$probability <- NULL
  d$prediction <- c("R", "I", "R", "S", "R", "S", "S", "S", "S", "S")
  expect_equal(unique(amr_import_predictions(d)$prediction_type), "sir")
})

test_that("MIC predictions are detected and coerced", {
  d <- toy_scores()
  d$probability <- NULL
  d$truth_mic <- c(">8", "16", "8", "4", "0.5", "<=0.125", "0.25", "<=0.125", "0.5", "1")
  d$predicted_mic <- c("16", "8", "8", "1", "0.5", "<=0.125", "1", "0.25", "0.5", "1")
  x <- amr_import_predictions(d)
  expect_equal(unique(x$prediction_type), "mic")
  expect_s3_class(x$predicted_mic, "mic")
})

test_that("fuzzy matches are rejected in strict mode and reported otherwise", {
  d <- toy_scores()
  d$drug[1] <- "notadrugatall"
  d$species[2] <- "zzzzbug"
  x <- amr_import_predictions(d)
  expect_equal(nrow(x), 10)
  expect_true(is.na(x$drug[1]))
  expect_true(is.na(x$species[2]))
  co <- attr(x, "coercion")
  expect_equal(co$status[co$input == "notadrugatall"], "rejected")
  v <- validate_amr_data(x)
  expect_false(v$ok)
  expect_true(all(c("drug", "species") %in% v$issues$check))

  y <- amr_import_predictions(d, strict = FALSE)
  expect_false(is.na(y$drug[1]))
  expect_true("approximate_match" %in% validate_amr_data(y)$issues$check)
})

test_that("exact names, codes, abbreviations and short names are accepted", {
  d <- toy_scores()[1:5, ]
  d$drug <- c("meropenem", "MEM", "Meropenem", "ciprofloxacin", "CIP")
  d$species <- c("Klebsiella pneumoniae", "K. pneumoniae", "K pneumoniae", "E. coli", "B_ESCHR_COLI")
  co <- attr(amr_import_predictions(d), "coercion")
  expect_true(all(co$status == "exact"))
})

test_that("non-data-frame input errors", {
  expect_error(amr_import_predictions(1:3), "data frame")
})
