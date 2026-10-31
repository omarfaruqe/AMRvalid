#' Antimicrobial susceptibility testing agreement metrics
#'
#' Computes the agreement metrics used to compare a susceptibility test (or
#' a genome-based prediction) with a reference method.
#'
#' @inheritParams amr_metrics
#' @param vme_denominator Denominator for very major and major error rates:
#'   `"class"` (default) divides very major errors by reference-resistant
#'   isolates and major errors by reference-susceptible isolates;
#'   `"all"` divides both by all evaluable isolates. Guidance documents
#'   differ, so always report which one you used.
#'
#' @details
#' Category metrics compare predicted and reference categories (S, SDD, I,
#' R); NI and missing values are excluded:
#'
#' * **Categorical agreement (CA)**: predicted category equals reference.
#' * **Very major error (VME)**: predicted S, reference R.
#' * **Major error (ME)**: predicted R, reference S.
#' * **Minor error (mE)**: one of the two is I or SDD and the other is S
#'   or R. The denominator is all evaluable isolates.
#'
#' For scores, the category is R when the probability is at or above
#' `threshold`, otherwise S.
#'
#' MIC metrics need `predicted_mic` and `truth_mic`:
#'
#' * **Essential agreement (EA)**: predicted MIC within one doubling
#'   dilution of the reference MIC. Off-scale reference values are handled
#'   conservatively: a reference of `<=a` agrees with predictions `<= 2a`,
#'   and `>a` agrees with predictions `>= a`. Predicted off-scale values are
#'   compared at their bound.
#' * **Bias**: mean of log2(predicted / reference) over on-scale pairs;
#'   positive values mean predictions are too high.
#'
#' Acceptance thresholds (for example from ISO 20776-2 or regulatory
#' guidance) are not applied here; compare the interval bounds against the
#' criteria that apply to your setting.
#'
#' @return A data frame with columns `metric`, `estimate`, `lower`,
#'   `upper`, `numerator` and `denominator`.
#' @export
#' @examples
#' preds <- data.frame(
#'   isolate_id = sprintf("iso%02d", 1:6),
#'   species = "E. coli",
#'   drug = "ciprofloxacin",
#'   truth = c("R", "R", "I", "S", "S", "S"),
#'   prediction = c("R", "S", "R", "S", "R", "S")
#' )
#' amr_ast_agreement(amr_import_predictions(preds))
amr_ast_agreement <- function(x, threshold = NULL, vme_denominator = c("class", "all"),
                              ci = c("wilson", "clopper-pearson"), conf = 0.95) {
  check_amr_data(x)
  vme_denominator <- match.arg(vme_denominator)
  ci <- match.arg(ci)
  ptype <- x$prediction_type[1]
  check_threshold(threshold, ptype)
  out <- list()
  add_prop <- function(name, num, den) {
    tab <- amr_prop_ci(num, den, method = ci, conf = conf)
    out[[length(out) + 1]] <<- metric_row(name, tab$estimate, tab$lower, tab$upper, num, den)
  }

  if (!identical(ptype, "mic")) {
    if (identical(ptype, "score")) {
      if (is.null(threshold)) {
        stop("Score predictions need a `threshold` to define categories.", call. = FALSE)
      }
      pred <- ifelse(is.na(x$probability), NA_character_, ifelse(x$probability >= threshold, "R", "S"))
    } else {
      pred <- as.character(x$prediction)
    }
    truth <- as.character(x$truth_sir)
    cats <- c("S", "SDD", "I", "R")
    ok <- truth %in% cats & pred %in% cats
    t <- truth[ok]; p <- pred[ok]
    n <- length(t)
    mid <- c("I", "SDD")
    vme <- sum(t == "R" & p == "S")
    me <- sum(t == "S" & p == "R")
    minor <- sum((t %in% mid & p %in% c("S", "R")) | (t %in% c("S", "R") & p %in% mid))
    out[[1]] <- metric_row("n", n, NA, NA, n, NA)
    add_prop("categorical_agreement", sum(t == p), n)
    add_prop("very_major_error", vme, if (vme_denominator == "class") sum(t == "R") else n)
    add_prop("major_error", me, if (vme_denominator == "class") sum(t == "S") else n)
    add_prop("minor_error", minor, n)
  }

  if (all(c("predicted_mic", "truth_mic") %in% names(x))) {
    ea <- essential_agreement(x$predicted_mic, x$truth_mic)
    add_prop("essential_agreement", ea$agree, ea$n)
    out[[length(out) + 1]] <- metric_row("bias_log2", ea$bias, NA, NA, NA, ea$n_on_scale)
  } else if (identical(ptype, "mic")) {
    stop("MIC predictions need a 'truth_mic' column for essential agreement.", call. = FALSE)
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "settings") <- list(threshold = threshold, vme_denominator = vme_denominator,
                                ci = ci, conf = conf)
  res
}

# Parse an AMR 'mic' vector into numeric value and censoring operator.
parse_mic <- function(m) {
  s <- gsub("\\s", "", as.character(m))
  op <- ifelse(grepl("^<", s), "<", ifelse(grepl("^>", s), ">", "="))
  val <- suppressWarnings(as.numeric(gsub("^[<>=]+", "", s)))
  list(value = val, op = op)
}

essential_agreement <- function(pred, truth) {
  p <- parse_mic(pred); t <- parse_mic(truth)
  ok <- !is.na(p$value) & !is.na(t$value)
  pv <- p$value[ok]; tv <- t$value[ok]; top <- t$op[ok]; pop <- p$op[ok]
  agree <- ifelse(top == "<", pv <= 2 * tv,
                  ifelse(top == ">", pv >= tv,
                         abs(log2(pv) - log2(tv)) <= 1 + 1e-9))
  on_scale <- top == "=" & pop == "="
  bias <- if (any(on_scale)) mean(log2(pv[on_scale]) - log2(tv[on_scale])) else NA_real_
  list(agree = sum(agree), n = length(agree), bias = bias, n_on_scale = sum(on_scale))
}
