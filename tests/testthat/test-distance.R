# Training isolates sit at positions 1-20 and test isolates at 21-60 on a
# line, so a test isolate's nearest-training distance is its position - 20.
# Predictions are correct for near isolates and increasingly wrong far away.
distance_toy <- function() {
  pos <- 1:60
  ids <- sprintf("iso%02d", pos)
  truth <- rep(c("R", "S"), 30)
  set.seed(42)
  wrong <- pos > 20 & stats::runif(60) < pmin(0.9, (pos - 20) / 45)
  prob <- ifelse(truth == "R", 0.85, 0.15)
  prob[wrong] <- 1 - prob[wrong]
  preds <- data.frame(isolate_id = ids, species = "E. coli", drug = "ciprofloxacin",
                      truth = truth, probability = prob,
                      split = ifelse(pos <= 20, "train", "test"))
  m <- as.matrix(stats::dist(pos))
  dimnames(m) <- list(ids, ids)
  list(preds = preds, matrix = m, pos = pos, ids = ids)
}

test_that("nearest-training distances agree for matrix and long-table input", {
  toy <- distance_toy()
  te <- toy$ids[21:60]; tr <- toy$ids[1:20]
  a <- amr_nearest_train_distance(te, tr, toy$matrix)
  expect_equal(a$distance, 1:40)
  expect_true(all(a$nearest_train_id == "iso20"))
  long <- as.data.frame(as.table(toy$matrix), stringsAsFactors = FALSE)
  names(long) <- c("isolate_1", "isolate_2", "distance")
  long <- long[long$isolate_1 < long$isolate_2, ]
  b <- amr_nearest_train_distance(te, tr, long)
  expect_equal(b$distance, a$distance)
})

test_that("isolates in training get distance 0 and unknown isolates NA", {
  toy <- distance_toy()
  r <- amr_nearest_train_distance(c("iso05", "unknown"), toy$ids[1:20], toy$matrix)
  expect_equal(r$distance[1], 0)
  expect_true(is.na(r$distance[2]))
  expect_error(amr_nearest_train_distance("a", "b", matrix(1)), "row and column names")
  expect_error(amr_nearest_train_distance("a", "b", data.frame(x = 1)), "isolate_1")
})

test_that("band evaluation and trend detect decay with distance", {
  toy <- distance_toy()
  x <- amr_import_predictions(toy$preds)
  r <- amr_distance_eval(x, toy$matrix, threshold = 0.5)
  expect_s3_class(r, "amr_distance_eval")
  expect_equal(nlevels(r$nearest$band), 3)
  expect_equal(sum(table(r$nearest$band)), 40)
  acc <- r$by_band$estimate[r$by_band$metric == "accuracy"]
  expect_true(acc[1] > acc[3])
  expect_true(r$trend$effect < 1)
  expect_true(r$trend$upper < 1)
  expect_output(print(r), "Trend with distance")
})

test_that("scores without a threshold use the Brier trend", {
  toy <- distance_toy()
  r <- amr_distance_eval(amr_import_predictions(toy$preds), toy$matrix)
  expect_match(r$trend$model, "squared error")
  expect_true(r$trend$effect > 0)
})

test_that("custom band boundaries are used and must cover all distances", {
  toy <- distance_toy()
  x <- amr_import_predictions(toy$preds)
  r <- amr_distance_eval(x, toy$matrix, bands = c(0, 5, 20, Inf), threshold = 0.5)
  expect_equal(as.vector(table(r$nearest$band)), c(4, 15, 21))
  expect_error(amr_distance_eval(x, toy$matrix, bands = c(2, 10, Inf), threshold = 0.5), "cover")
})

test_that("precomputed nearest distances and training manifests work", {
  toy <- distance_toy()
  p <- toy$preds[21:60, ]
  p$split <- NULL
  pc <- data.frame(isolate_id = p$isolate_id, nearest_distance = 1:40)
  r1 <- amr_distance_eval(amr_import_predictions(p), pc, threshold = 0.5)
  expect_equal(r1$nearest$distance, 1:40)
  r2 <- amr_distance_eval(amr_import_predictions(p), toy$matrix,
                          train_manifest = data.frame(isolate_id = toy$ids[1:20]), threshold = 0.5)
  expect_equal(r2$nearest$distance, 1:40)
  expect_error(amr_distance_eval(amr_import_predictions(p), toy$matrix), "No training isolates")
})

test_that("isolates without distances are excluded with a warning", {
  toy <- distance_toy()
  m <- toy$matrix[-60, -60]
  expect_warning(r <- amr_distance_eval(amr_import_predictions(toy$preds), m, threshold = 0.5),
                 "1 test isolate")
  expect_equal(sum(!is.na(r$nearest$band)), 39)
})

test_that("weight_by is passed to band metrics and weights the trend", {
  toy <- distance_toy()
  d <- toy$preds
  # clusters of 4 consecutive isolates
  d$genomic_cluster <- paste0("C", (seq_len(nrow(d)) - 1) %/% 4)
  x <- amr_import_predictions(d)
  r0 <- amr_distance_eval(x, toy$matrix, threshold = 0.5)
  r1 <- amr_distance_eval(x, toy$matrix, threshold = 0.5, weight_by = "genomic_cluster")
  expect_true("n_weight_groups" %in% r1$by_band$metric)
  expect_false("n_weight_groups" %in% r0$by_band$metric)
  # weighting by clusters of 4 makes the trend interval wider (fewer units)
  w0 <- log(r0$trend$upper) - log(r0$trend$lower)
  w1 <- log(r1$trend$upper) - log(r1$trend$lower)
  expect_gt(w1, 1.5 * w0)
  expect_output(print(r1), "weighted by genomic_cluster")
  r2 <- amr_distance_eval(x, toy$matrix, threshold = 0.5, target_prevalence = 0.1)
  expect_true("ppv_at_target_prevalence" %in% r2$by_band$metric)
  expect_error(amr_distance_eval(x, toy$matrix, weight_by = "nope"), "not found")
})

test_that("unweighted distance results are unchanged by the weighting code", {
  toy <- distance_toy()
  r <- amr_distance_eval(amr_import_predictions(toy$preds), toy$matrix, threshold = 0.5)
  y <- as.integer(toy$preds$truth == ifelse(toy$preds$probability >= 0.5, "R", "S"))[21:60]
  dz <- (1:40) / stats::sd(1:40)
  ref <- stats::glm(y ~ dz, family = stats::binomial())
  expect_equal(r$trend$effect, unname(exp(stats::coef(ref)[2])))
})
