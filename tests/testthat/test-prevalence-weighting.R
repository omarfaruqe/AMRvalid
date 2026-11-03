get <- function(m, name) m$estimate[m$metric == name]

test_that("amr_ppv_npv matches Bayes' rule", {
  r <- amr_ppv_npv(0.9, 0.8, c(0.1, 0.5))
  expect_equal(r$ppv, c(0.9 * 0.1 / (0.9 * 0.1 + 0.2 * 0.9), 0.9 * 0.5 / (0.9 * 0.5 + 0.2 * 0.5)))
  expect_equal(r$npv[1], 0.8 * 0.9 / (0.8 * 0.9 + 0.1 * 0.1))
  expect_error(amr_ppv_npv(1.2, 0.8, 0.1), "lie in")
})

test_that("target prevalence adds adjusted PPV and NPV with conservative bounds", {
  x <- amr_import_predictions(toy_scores())
  m <- amr_metrics(x, threshold = 0.5, target_prevalence = 0.1)
  se <- 3 / 4; sp <- 5 / 6
  expect_equal(get(m, "ppv_at_target_prevalence"), se * 0.1 / (se * 0.1 + (1 - sp) * 0.9))
  row <- m[m$metric == "ppv_at_target_prevalence", ]
  sl <- m$lower[m$metric == "sensitivity"]; pl <- m$lower[m$metric == "specificity"]
  expect_equal(row$lower, amr_ppv_npv(sl, pl, 0.1)$ppv)
  expect_true(row$lower <= row$estimate && row$estimate <= row$upper)
  expect_error(amr_metrics(x, threshold = 0.5, target_prevalence = 1), "strictly")
})

test_that("unweighted results are unchanged by the weighting code", {
  x <- amr_import_predictions(toy_scores())
  m <- amr_metrics(x, threshold = 0.5)
  expect_equal(get(m, "sensitivity"), 3 / 4)
  expect_equal(m$lower[m$metric == "sensitivity"], amr_prop_ci(3, 4)$lower)
  expect_false("n_weight_groups" %in% m$metric)
})

# Lineage A: 8 resistant isolates, all predicted correctly (an outbreak clone).
# Lineages B-E: one resistant isolate each, all missed.
outbreak_toy <- function() {
  data.frame(
    isolate_id = sprintf("iso%02d", 1:16),
    species = "K. pneumoniae", drug = "meropenem",
    truth = c(rep("R", 12), rep("S", 4)),
    prediction = c(rep("R", 8), rep("S", 4), rep("S", 4)),
    lineage = c(rep("A", 8), "B", "C", "D", "E", "F", "G", "H", "I")
  )
}

test_that("cluster weighting stops an outbreak clone dominating", {
  x <- amr_import_predictions(outbreak_toy())
  m_iso <- amr_metrics(x)
  m_lin <- amr_metrics(x, weight_by = "lineage")
  expect_equal(get(m_iso, "sensitivity"), 8 / 12)
  expect_equal(get(m_lin, "sensitivity"), 1 / 5)
  expect_equal(get(m_lin, "n_weight_groups"), 9)
  expect_equal(m_lin$denominator[m_lin$metric == "sensitivity"], 5)
  expect_equal(get(m_lin, "prevalence_R"), 5 / 9)
  a <- amr_ast_agreement(x, weight_by = "lineage")
  expect_equal(get(a, "very_major_error"), 4 / 5)
  expect_equal(get(amr_ast_agreement(x), "very_major_error"), 4 / 12)
})

test_that("weighted intervals use the effective sample size", {
  x <- amr_import_predictions(outbreak_toy())
  m <- amr_metrics(x, weight_by = "lineage")
  # resistant rows: 8 with weight 1/8 and 4 with weight 1 -> n_eff = 25 / 4.125
  neff <- 5^2 / (8 * (1 / 8)^2 + 4)
  ref <- amr_prop_ci(0.2 * neff, neff)
  expect_equal(m$lower[m$metric == "sensitivity"], ref$lower)
})

test_that("weighted AUC equals the unweighted AUC when every group has one isolate", {
  d <- toy_scores()
  d$lineage <- paste0("L", seq_len(nrow(d)))
  x <- amr_import_predictions(d)
  expect_equal(get(amr_metrics(x, weight_by = "lineage"), "roc_auc"), get(amr_metrics(x), "roc_auc"))
  expect_true(is.na(amr_metrics(x, weight_by = "lineage")$lower[amr_metrics(x, weight_by = "lineage")$metric == "roc_auc"]))
})

test_that("amr_evaluate passes weighting and prevalence through", {
  x <- amr_import_predictions(outbreak_toy())
  e <- amr_evaluate(x, weight_by = "lineage", target_prevalence = 0.2)
  expect_equal(e$estimate[e$metric == "sensitivity"], 1 / 5)
  expect_true("ppv_at_target_prevalence" %in% e$metric)
  expect_output(print(e), "weighted by lineage")
  expect_error(amr_metrics(x, weight_by = "nope"), "not found")
})
