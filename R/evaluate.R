#' Evaluate predictions per stratum
#'
#' Splits an `amr_data` object by model, dataset, species and antimicrobial
#' (and optionally further columns) and computes [amr_metrics()] and, where
#' possible, [amr_ast_agreement()] in each stratum. There is no pooled
#' headline score.
#'
#' @inheritParams amr_metrics
#' @inheritParams amr_ast_agreement
#' @param split Which partition to evaluate when a `split` column exists.
#'   Defaults to `"test"`, so training predictions are not reported as
#'   performance by accident. Use `NULL` to evaluate all rows.
#' @param validate If `TRUE` (default), run [amr_validate()] first and
#'   stop if it reports errors.
#' @param by Optional column name(s) to evaluate subgroups within each
#'   stratum, e.g. `"country"` or a period column you have added. Missing
#'   values form their own subgroup, labelled `"(missing)"`. Factor columns
#'   keep their level order in the output. Subgroup results
#'   are descriptive: small subgroups give wide intervals, and differences
#'   between subgroups are not evidence of a cause.
#' @param target_prevalence,weight_by Passed to [amr_metrics()] (and
#'   `weight_by` to [amr_ast_agreement()]): PPV and NPV at a target
#'   prevalence, and cluster-weighted metrics.
#'
#' @return An object of class `amr_eval`: a long data frame with the
#'   stratum columns (and any `by` columns), `species_name`, `drug_name`, `family`
#'   (`"ml"` or `"ast"`) and the metric columns.
#' @export
#' @examples
#' preds <- data.frame(
#'   isolate_id = sprintf("iso%02d", 1:8),
#'   species = "K. pneumoniae",
#'   drug = rep(c("meropenem", "ciprofloxacin"), each = 4),
#'   truth = c("R", "S", "S", "R", "R", "R", "S", "S"),
#'   probability = c(0.9, 0.2, 0.4, 0.7, 0.8, 0.3, 0.1, 0.6)
#' )
#' x <- amr_import_predictions(preds)
#' amr_evaluate(x, threshold = 0.5)
#'
#' # subgroups, e.g. by country of isolation
#' x$country <- c("CA", "CA", "US", "US", "CA", "US", "CA", "US")
#' amr_evaluate(x, threshold = 0.5, by = "country")
amr_evaluate <- function(x, threshold = NULL, intermediate = c("exclude", "susceptible", "resistant"),
                         vme_denominator = c("class", "all"),
                         ci = c("wilson", "clopper-pearson"), conf = 0.95,
                         split = "test", validate = TRUE, target_prevalence = NULL,
                         weight_by = NULL, by = NULL) {
  check_amr_data(x)
  if (!is.null(by)) {
    if (!is.character(by)) stop("`by` must be a character vector of column names.", call. = FALSE)
    missing_by <- setdiff(by, names(x))
    if (length(missing_by)) stop("`by` column(s) not found: ", paste(missing_by, collapse = ", "), call. = FALSE)
    clash <- intersect(by, c(stratum_columns(x), "metric", "estimate"))
    if (length(clash)) stop("`by` cannot use ", paste(clash, collapse = ", "), "; these are already strata.",
                            call. = FALSE)
  }
  intermediate <- match.arg(intermediate)
  vme_denominator <- match.arg(vme_denominator)
  ci <- match.arg(ci)
  if (validate) {
    v <- amr_validate(x, min_n = 0, min_resistant = 0)
    if (!v$ok) {
      print(v)
      stop("Validation failed; fix the errors above or set `validate = FALSE`.", call. = FALSE)
    }
  }
  if (!is.null(split) && "split" %in% names(x)) {
    x <- x[x$split %in% split, , drop = FALSE]
    if (nrow(x) == 0) stop("No rows with split = '", paste(split, collapse = "', '"), "'.", call. = FALSE)
  }
  ptype <- x$prediction_type[1]
  if (identical(ptype, "score") && is.null(threshold)) {
    message("No `threshold` given: reporting threshold-free metrics only (AUC, Brier score, log loss).")
  }

  for (b in by) {
    lev <- if (is.factor(x[[b]])) levels(x[[b]]) else sort(unique(stats::na.omit(as.character(x[[b]]))))
    v <- as.character(x[[b]])
    v[is.na(v) | v == ""] <- "(missing)"
    x[[b]] <- factor(v, levels = unique(c(setdiff(lev, ""), "(missing)")))
  }
  cols <- c(stratum_columns(x), by)
  key <- interaction(as.data.frame(x)[, cols, drop = FALSE], drop = TRUE, lex.order = TRUE, sep = "\r")
  pieces <- lapply(split(seq_len(nrow(x)), key), function(idx) {
    sub <- x[idx, , drop = FALSE]
    res <- list()
    if (!identical(ptype, "mic")) {
      m <- amr_metrics(sub, threshold = threshold, intermediate = intermediate, ci = ci, conf = conf,
                       target_prevalence = target_prevalence, weight_by = weight_by, split = NULL)
      m$family <- "ml"
      res[[1]] <- m
    }
    has_categories <- !identical(ptype, "score") || !is.null(threshold)
    if (has_categories && (!identical(ptype, "mic") || "truth_mic" %in% names(sub))) {
      a <- amr_ast_agreement(sub, threshold = threshold, vme_denominator = vme_denominator, weight_by = weight_by,
                             split = NULL,
                             ci = ci, conf = conf)
      a <- a[a$metric != "n", , drop = FALSE]
      a$family <- "ast"
      res[[length(res) + 1]] <- a
    }
    met <- do.call(rbind, res)
    first <- as.data.frame(sub)[1, cols, drop = FALSE]
    first[] <- lapply(first, as.character)
    cbind(first[rep(1, nrow(met)), , drop = FALSE], met, row.names = NULL)
  })
  out <- do.call(rbind, pieces)
  out$species_name <- suppressMessages(AMR::mo_name(out$species, language = "en"))
  out$drug_name <- suppressMessages(AMR::ab_name(out$drug, language = "en"))
  front <- c(cols, "species_name", "drug_name", "family", "metric")
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]
  rownames(out) <- NULL
  structure(out, class = c("amr_eval", "data.frame"),
            settings = list(threshold = threshold, intermediate = intermediate,
                            vme_denominator = vme_denominator, ci = ci, conf = conf,
                            split = split, prediction_type = ptype,
                            target_prevalence = target_prevalence, weight_by = weight_by, by = by))
}

#' @export
print.amr_eval <- function(x, digits = 3, ...) {
  if (!all(c("model_id", "dataset_id", "species", "drug", "species_name", "drug_name", "metric", "estimate", "lower", "upper", "numerator", "denominator") %in% names(x)) || is.null(attr(x, "settings"))) {
    return(print(as.data.frame(unclass_df(x)), ...))
  }
  s <- attr(x, "settings")
  strata_cols <- intersect(c("model_id", "dataset_id", "species", "drug", s$by), names(x))
  cat("<amr_eval> ", nrow(unique(as.data.frame(x)[, strata_cols, drop = FALSE])),
      " stratum/strata | type: ", s$prediction_type,
      if (!is.null(s$by)) paste0(" | by ", paste(s$by, collapse = ", ")),
      " | threshold: ", if (is.null(s$threshold)) "none" else s$threshold,
      " | I/SDD: ", s$intermediate,
      if (!is.null(s$weight_by)) paste0(" | weighted by ", s$weight_by),
      if (!is.null(s$target_prevalence)) paste0(" | target prevalence ", s$target_prevalence),
      " | ", s$conf * 100, "% ", s$ci, " intervals\n", sep = "")
  d <- as.data.frame(x)
  fmt <- function(v) ifelse(is.na(v), "", formatC(v, digits = digits, format = "f"))
  shown <- data.frame(
    model = d$model_id, species = d$species_name, drug = d$drug_name,
    d[, intersect(s$by, names(d)), drop = FALSE],
    metric = d$metric,
    estimate = fmt(d$estimate),
    interval = ifelse(is.na(d$lower), "", paste0("[", fmt(d$lower), ", ", fmt(d$upper), "]")),
    count = ifelse(is.na(d$denominator), "",
                   paste0(ifelse(is.na(d$numerator), "", paste0(d$numerator, "/")), d$denominator)),
    stringsAsFactors = FALSE
  )
  print(shown, row.names = FALSE, right = FALSE)
  invisible(x)
}
