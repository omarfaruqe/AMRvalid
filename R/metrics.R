#' Classification and probability metrics for one set of predictions
#'
#' Computes machine-learning performance metrics for the rows of an
#' `amr_data` object, treating resistance (R) as the positive class. Use
#' [amr_species_eval()] to compute them per species-antimicrobial stratum;
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
#'
#' @return A data frame with one row per metric and columns `metric`,
#'   `estimate`, `lower`, `upper`, `numerator` and `denominator`. Metrics
#'   without a closed-form interval (balanced accuracy, F1, MCC, Brier
#'   score, log loss) have `NA` bounds; group-aware bootstrap intervals are
#'   planned.
#'
#' @details
#' ROC-AUC is computed with the Mann-Whitney statistic (ties count one
#' half) and its interval with the DeLong variance. Row-level intervals
#' assume isolates are independent, which is optimistic when isolates are
#' clonally related.
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
#'   probability = c(0.91, 0.75, 0.30, 0.10, 0.20, 0.65, 0.05, 0.88)
#' )
#' amr_metrics(amr_import_predictions(preds), threshold = 0.5)
amr_metrics <- function(x, threshold = NULL, intermediate = c("exclude", "susceptible", "resistant"),
                        ci = c("wilson", "clopper-pearson"), conf = 0.95) {
  check_amr_data(x)
  intermediate <- match.arg(intermediate)
  ci <- match.arg(ci)
  ptype <- x$prediction_type[1]
  if (identical(ptype, "mic")) {
    stop("MIC predictions have no binary class without breakpoints; use amr_ast_agreement() ",
         "for essential agreement, or interpret predicted MICs into S/I/R first (e.g. AMR::as.sir()).",
         call. = FALSE)
  }
  check_threshold(threshold, ptype)

  truth <- collapse_sir(x$truth_sir, intermediate)
  pred <- predicted_class(x, threshold, intermediate)
  prob <- if (identical(ptype, "score")) x$probability else NULL

  out <- list()
  keep <- !is.na(truth) & (if (is.null(prob)) TRUE else !is.na(prob))
  y <- truth[keep]
  n_r <- sum(y == "R"); n_s <- sum(y == "S")
  out[[1]] <- metric_row("n", length(y), NA, NA, length(y), NA)
  out[[2]] <- metric_row("n_resistant", n_r, NA, NA, n_r, NA)
  out[[3]] <- metric_row("prevalence_R", if (length(y)) n_r / length(y) else NA, NA, NA, n_r, length(y))

  if (!is.null(pred)) {
    ok <- !is.na(truth) & !is.na(pred)
    t <- truth[ok]; p <- pred[ok]
    tp <- sum(t == "R" & p == "R"); fn <- sum(t == "R" & p == "S")
    tn <- sum(t == "S" & p == "S"); fp <- sum(t == "S" & p == "R")
    props <- list(
      sensitivity = c(tp, tp + fn),
      specificity = c(tn, tn + fp),
      ppv = c(tp, tp + fp),
      npv = c(tn, tn + fn),
      accuracy = c(tp + tn, tp + tn + fp + fn)
    )
    for (nm in names(props)) {
      v <- props[[nm]]
      ci_tab <- amr_prop_ci(v[1], v[2], method = ci, conf = conf)
      out[[length(out) + 1]] <- metric_row(nm, ci_tab$estimate, ci_tab$lower, ci_tab$upper, v[1], v[2])
    }
    sens <- safe_div(tp, tp + fn); spec <- safe_div(tn, tn + fp)
    prec <- safe_div(tp, tp + fp)
    out[[length(out) + 1]] <- metric_row("balanced_accuracy", (sens + spec) / 2, NA, NA, NA, NA)
    out[[length(out) + 1]] <- metric_row("f1", safe_div(2 * prec * sens, prec + sens), NA, NA, NA, NA)
    mcc_den <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    out[[length(out) + 1]] <- metric_row("mcc", safe_div(as.numeric(tp) * tn - as.numeric(fp) * fn, mcc_den),
                                         NA, NA, NA, NA)
  }

  if (!is.null(prob)) {
    ok <- truth %in% c("R", "S") & !is.na(prob)
    yy <- as.integer(truth[ok] == "R"); pp <- prob[ok]
    auc <- auc_delong(pp, yy, conf)
    out[[length(out) + 1]] <- metric_row("roc_auc", auc[1], auc[2], auc[3], NA, length(yy))
    brier <- if (length(yy)) mean((pp - yy)^2) else NA_real_
    eps <- 1e-15
    pc <- pmin(pmax(pp, eps), 1 - eps)
    logloss <- if (length(yy)) -mean(yy * log(pc) + (1 - yy) * log(1 - pc)) else NA_real_
    out[[length(out) + 1]] <- metric_row("brier", brier, NA, NA, NA, length(yy))
    out[[length(out) + 1]] <- metric_row("log_loss", logloss, NA, NA, NA, length(yy))
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "settings") <- list(threshold = threshold, intermediate = intermediate, ci = ci, conf = conf)
  res
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
