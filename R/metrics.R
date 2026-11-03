#' Classification and probability metrics for one set of predictions
#'
#' Computes machine-learning performance metrics for the rows of an
#' `amr_data` object, treating resistance (R) as the positive class. Use
#' [amr_evaluate()] to compute them per species-antimicrobial stratum;
#' pooling strata is rarely meaningful.
#'
#' @param x An `amr_data` object.
#' @param threshold Probability at or above which a score is called R.
#'   Required for classification metrics when predictions are scores; there
#'   is deliberately no default of 0.5. Choose it on training or validation
#'   data, never on the test set.
#' @param intermediate How to treat I and SDD in binary metrics, for both
#'   reference and prediction: `"exclude"` (default) drops those rows,
#'   `"susceptible"` counts them as S and `"resistant"` as R. NI results are
#'   always excluded.
#' @param ci Interval method for proportions, passed to [amr_prop_ci()].
#' @param conf Confidence level.
#' @param target_prevalence Optional proportion resistant in the setting
#'   where the predictor will be used. Adds `ppv_at_target_prevalence` and
#'   `npv_at_target_prevalence`, computed from sensitivity and specificity,
#'   because PPV and NPV depend on prevalence and public collections often
#'   over-represent resistant isolates.
#' @param split Which partition to evaluate when a `split` column exists.
#'   Defaults to `"test"`, so training predictions are never counted as
#'   performance by accident. Use `NULL` to evaluate all rows.
#' @param weight_by Optional column (e.g. `"genomic_cluster"` or
#'   `"lineage"`) for cluster-weighted metrics: every group contributes a
#'   total weight of 1, shared equally by its isolates, so a large outbreak
#'   clone counts as much as a single sporadic isolate. Rows with a missing
#'   label are their own group.
#'
#' @return A data frame with one row per metric and columns `metric`,
#'   `estimate`, `lower`, `upper`, `numerator` and `denominator`. Metrics
#'   without a closed-form interval (balanced accuracy, F1, MCC, Brier
#'   score, log loss) have `NA` bounds; use [amr_bootstrap()] for
#'   group-aware intervals.
#'
#' @details
#' ROC-AUC is computed with the Mann-Whitney statistic (ties count one
#' half) and its interval with the DeLong variance. Row-level intervals
#' assume isolates are independent, which is optimistic when isolates are
#' clonally related.
#'
#' Prevalence-adjusted values use
#' PPV = sens x p / (sens x p + (1 - spec) x (1 - p)) and
#' NPV = spec x (1 - p) / (spec x (1 - p) + (1 - sens) x p). Both increase
#' with sensitivity and specificity, so their bounds combine the lower (and
#' upper) bounds of the two; these bounds are conservative.
#'
#' With `weight_by`, numerators and denominators are sums of weights, and
#' proportion intervals use the Kish effective sample size
#' (sum of weights squared / sum of squared weights). The AUC becomes a
#' weighted Mann-Whitney statistic without a closed-form interval.
#'
#' @references
#' DeLong ER, DeLong DM, Clarke-Pearson DL (1988). Comparing the areas
#' under two or more correlated receiver operating characteristic curves: a
#' nonparametric approach. *Biometrics* 44:837-845.
#' @export
#' @examples
#' preds <- data.frame(
#'   isolate_id = sprintf("iso%02d", 1:8),
#'   species = "K. pneumoniae",
#'   drug = "meropenem",
#'   truth = c("R", "R", "R", "S", "S", "S", "S", "R"),
#'   probability = c(0.91, 0.75, 0.30, 0.10, 0.20, 0.65, 0.05, 0.88),
#'   lineage = c("ST258", "ST258", "ST258", "ST11", "ST15", "ST15", "ST307", "ST147")
#' )
#' x <- amr_import_predictions(preds)
#' amr_metrics(x, threshold = 0.5)
#' # PPV and NPV where 10% of isolates are resistant
#' amr_metrics(x, threshold = 0.5, target_prevalence = 0.1)
#' # each lineage counts once
#' amr_metrics(x, threshold = 0.5, weight_by = "lineage")
amr_metrics <- function(x, threshold = NULL, intermediate = c("exclude", "susceptible", "resistant"),
                        ci = c("wilson", "clopper-pearson"), conf = 0.95,
                        target_prevalence = NULL, weight_by = NULL, split = "test") {
  check_amr_data(x)
  x <- select_split(x, split)
  intermediate <- match.arg(intermediate)
  ci <- match.arg(ci)
  ptype <- x$prediction_type[1]
  if (identical(ptype, "mic")) {
    stop("MIC predictions have no binary class without breakpoints; use amr_ast_agreement() ",
         "for essential agreement, or interpret predicted MICs into S/I/R first (e.g. AMR::as.sir()).",
         call. = FALSE)
  }
  check_threshold(threshold, ptype)
  if (!is.null(target_prevalence) &&
      (!is.numeric(target_prevalence) || length(target_prevalence) != 1 ||
       target_prevalence <= 0 || target_prevalence >= 1)) {
    stop("`target_prevalence` must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  truth <- collapse_sir(x$truth_sir, intermediate)
  pred <- predicted_class(x, threshold, intermediate)
  prob <- if (identical(ptype, "score")) x$probability else NULL
  w_all <- cluster_weights(x, weight_by, !is.na(truth))
  weighted <- !is.null(weight_by)
  wprop <- function(num_mask, den_mask, w) {
    num <- sum(w[num_mask & den_mask]); den <- sum(w[den_mask])
    if (!weighted) return(c(num, den, amr_prop_ci(num, den, method = ci, conf = conf)))
    wd <- w[den_mask]
    neff <- if (den > 0) den^2 / sum(wd^2) else 0
    est <- if (den > 0) num / den else NA_real_
    tab <- amr_prop_ci(if (den > 0) est * neff else 0, neff, method = ci, conf = conf)
    tab$estimate <- est
    c(num, den, tab)
  }

  out <- list()
  keep <- !is.na(truth) & (if (is.null(prob)) TRUE else !is.na(prob))
  y <- truth[keep]
  n_r <- sum(y == "R")
  out[[1]] <- metric_row("n", length(y), NA, NA, length(y), NA)
  out[[2]] <- metric_row("n_resistant", n_r, NA, NA, n_r, NA)
  if (weighted) {
    n_groups <- length(unique(group_labels(x, weight_by)[keep]))
    out[[3]] <- metric_row("n_weight_groups", n_groups, NA, NA, n_groups, NA)
  }
  wk <- w_all[keep]
  prev <- if (length(y)) sum(wk[y == "R"]) / sum(wk) else NA
  out[[length(out) + 1]] <- metric_row("prevalence_R", prev, NA, NA, sum(wk[y == "R"]), sum(wk))

  sens_tab <- spec_tab <- NULL
  if (!is.null(pred)) {
    ok <- !is.na(truth) & !is.na(pred)
    t <- truth[ok]; p <- pred[ok]; w <- w_all[ok]
    props <- list(
      sensitivity = list(t == "R" & p == "R", t == "R"),
      specificity = list(t == "S" & p == "S", t == "S"),
      ppv = list(t == "R" & p == "R", p == "R"),
      npv = list(t == "S" & p == "S", p == "S"),
      accuracy = list(t == p, rep(TRUE, length(t)))
    )
    for (nm in names(props)) {
      v <- wprop(props[[nm]][[1]], props[[nm]][[2]], w)
      out[[length(out) + 1]] <- metric_row(nm, v$estimate, v$lower, v$upper, v[[1]], v[[2]])
      if (nm == "sensitivity") sens_tab <- v
      if (nm == "specificity") spec_tab <- v
    }
    tp <- sum(w[t == "R" & p == "R"]); fn <- sum(w[t == "R" & p == "S"])
    tn <- sum(w[t == "S" & p == "S"]); fp <- sum(w[t == "S" & p == "R"])
    sens <- safe_div(tp, tp + fn); spec <- safe_div(tn, tn + fp)
    prec <- safe_div(tp, tp + fp)
    out[[length(out) + 1]] <- metric_row("balanced_accuracy", (sens + spec) / 2, NA, NA, NA, NA)
    out[[length(out) + 1]] <- metric_row("f1", safe_div(2 * prec * sens, prec + sens), NA, NA, NA, NA)
    mcc_den <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    out[[length(out) + 1]] <- metric_row("mcc", safe_div(as.numeric(tp) * tn - as.numeric(fp) * fn, mcc_den),
                                         NA, NA, NA, NA)
    if (!is.null(target_prevalence)) {
      pv <- prevalence_adjusted(c(sens_tab$estimate, sens_tab$lower, sens_tab$upper),
                                c(spec_tab$estimate, spec_tab$lower, spec_tab$upper), target_prevalence)
      out[[length(out) + 1]] <- metric_row("ppv_at_target_prevalence", pv$ppv[1], pv$ppv[2], pv$ppv[3], NA, NA)
      out[[length(out) + 1]] <- metric_row("npv_at_target_prevalence", pv$npv[1], pv$npv[2], pv$npv[3], NA, NA)
    }
  }

  if (!is.null(prob)) {
    ok <- truth %in% c("R", "S") & !is.na(prob)
    yy <- as.integer(truth[ok] == "R"); pp <- prob[ok]; ww <- w_all[ok]
    if (weighted) {
      auc <- c(weighted_auc(pp, yy, ww), NA_real_, NA_real_)
    } else {
      auc <- auc_delong(pp, yy, conf)
    }
    out[[length(out) + 1]] <- metric_row("roc_auc", auc[1], auc[2], auc[3], NA, length(yy))
    eps <- 1e-15
    pc <- pmin(pmax(pp, eps), 1 - eps)
    brier <- if (length(yy)) stats::weighted.mean((pp - yy)^2, ww) else NA_real_
    logloss <- if (length(yy)) -stats::weighted.mean(yy * log(pc) + (1 - yy) * log(1 - pc), ww) else NA_real_
    out[[length(out) + 1]] <- metric_row("brier", brier, NA, NA, NA, length(yy))
    out[[length(out) + 1]] <- metric_row("log_loss", logloss, NA, NA, NA, length(yy))
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "settings") <- list(threshold = threshold, intermediate = intermediate, ci = ci, conf = conf,
                                target_prevalence = target_prevalence, weight_by = weight_by)
  res
}

#' PPV and NPV at a given prevalence
#'
#' Converts sensitivity and specificity into the positive and negative
#' predictive values expected where the proportion of resistant isolates is
#' `prevalence`.
#'
#' @param sensitivity,specificity Numeric vectors of values in \[0, 1\].
#' @param prevalence Proportion resistant in the target setting.
#' @return A data frame with `prevalence`, `ppv` and `npv`.
#' @export
#' @examples
#' # a test with 95% sensitivity and 90% specificity
#' amr_ppv_npv(0.95, 0.90, prevalence = c(0.01, 0.1, 0.5))
amr_ppv_npv <- function(sensitivity, specificity, prevalence) {
  if (any(c(sensitivity, specificity, prevalence) < 0 | c(sensitivity, specificity, prevalence) > 1,
          na.rm = TRUE)) {
    stop("Sensitivity, specificity and prevalence must lie in [0, 1].", call. = FALSE)
  }
  n <- max(length(sensitivity), length(specificity), length(prevalence))
  se <- rep_len(sensitivity, n); sp <- rep_len(specificity, n); p <- rep_len(prevalence, n)
  data.frame(prevalence = p,
             ppv = safe_div(se * p, se * p + (1 - sp) * (1 - p)),
             npv = safe_div(sp * (1 - p), sp * (1 - p) + (1 - se) * p))
}

# Point estimate and conservative bounds; sens and spec are c(est, lower, upper).
prevalence_adjusted <- function(sens, spec, prevalence) {
  r <- amr_ppv_npv(sens, spec, prevalence)
  list(ppv = r$ppv, npv = r$npv)
}

# Weight per row: 1 for every evaluable row, or 1 / (group size among
# evaluable rows) when `weight_by` is given.
cluster_weights <- function(x, weight_by, evaluable) {
  n <- nrow(x)
  if (is.null(weight_by)) return(rep(1, n))
  if (!weight_by %in% names(x)) stop("Weighting column '", weight_by, "' not found.", call. = FALSE)
  g <- group_labels(x, weight_by)
  w <- rep(0, n)
  sizes <- table(g[evaluable])
  w[evaluable] <- 1 / as.vector(sizes[g[evaluable]])
  w
}

# Weighted Mann-Whitney AUC: sum over (R, S) pairs of w_i w_j psi / sums.
weighted_auc <- function(p, y, w) {
  pos <- y == 1; neg <- y == 0
  if (!any(pos) || !any(neg)) return(NA_real_)
  psi <- outer(p[pos], p[neg], function(a, b) (a > b) + 0.5 * (a == b))
  sum(outer(w[pos], w[neg]) * psi) / (sum(w[pos]) * sum(w[neg]))
}

# Helpers -------------------------------------------------------------------

check_amr_data <- function(x) {
  if (!inherits(x, "amr_data")) {
    stop("`x` must be an amr_data object; use amr_import_predictions() first.", call. = FALSE)
  }
  if (!"truth_sir" %in% names(x) && !identical(x$prediction_type[1], "mic")) {
    stop("A 'truth_sir' column is required for these metrics.", call. = FALSE)
  }
  invisible(TRUE)
}

check_threshold <- function(threshold, ptype) {
  if (identical(ptype, "score") && !is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0 || threshold > 1) {
      stop("`threshold` must be a single number in [0, 1].", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Reduce S/SDD/I/R/NI to S/R (or NA) according to the intermediate policy.
collapse_sir <- function(v, intermediate) {
  v <- as.character(v)
  out <- ifelse(v %in% c("S", "R"), v, NA_character_)
  mid <- v %in% c("I", "SDD")
  if (intermediate == "susceptible") out[mid] <- "S"
  if (intermediate == "resistant") out[mid] <- "R"
  out
}

# Predicted S/R class, or NULL when scores are given without a threshold.
predicted_class <- function(x, threshold, intermediate) {
  ptype <- x$prediction_type[1]
  if (identical(ptype, "score")) {
    if (is.null(threshold)) return(NULL)
    return(ifelse(is.na(x$probability), NA_character_,
                  ifelse(x$probability >= threshold, "R", "S")))
  }
  collapse_sir(x$prediction, intermediate)
}

safe_div <- function(a, b) ifelse(b == 0 | is.na(b), NA_real_, a / b)

metric_row <- function(metric, estimate, lower, upper, numerator, denominator) {
  data.frame(metric = metric, estimate = as.numeric(estimate), lower = as.numeric(lower),
             upper = as.numeric(upper), numerator = as.numeric(numerator),
             denominator = as.numeric(denominator), stringsAsFactors = FALSE)
}

# ROC-AUC with DeLong confidence interval. Returns c(auc, lower, upper).
auc_delong <- function(score, y, conf) {
  pos <- score[y == 1]; neg <- score[y == 0]
  m <- length(pos); n <- length(neg)
  if (m == 0 || n == 0) return(c(NA_real_, NA_real_, NA_real_))
  psi <- function(a, b) (a > b) + 0.5 * (a == b)
  v10 <- vapply(pos, function(p) mean(psi(p, neg)), numeric(1))
  v01 <- vapply(neg, function(q) mean(psi(pos, q)), numeric(1))
  auc <- mean(v10)
  if (m < 2 || n < 2) return(c(auc, NA_real_, NA_real_))
  se <- sqrt(stats::var(v10) / m + stats::var(v01) / n)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  c(auc, max(0, auc - z * se), min(1, auc + z * se))
}
