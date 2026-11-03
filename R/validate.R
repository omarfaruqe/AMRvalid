#' Validate an AMRvalid prediction table
#'
#' Runs structural and scientific checks before evaluation. Blocking problems
#' are reported as errors; problems that need attention but do not prevent
#' analysis are warnings; notes that affect interpretation are info.
#'
#' @param x An `amr_data` object from [amr_import_predictions()].
#' @param min_n Strata with fewer evaluable isolates than this trigger a
#'   warning.
#' @param min_resistant Strata with fewer reference-resistant isolates than
#'   this trigger a warning, because very major error rates are then
#'   unstable.
#'
#' @return An object of class `amr_validation`: a list with `ok` (`TRUE` if
#'   there are no errors) and `issues` (a data frame with columns `level`,
#'   `check`, `message` and `n_rows`).
#' @export
#' @examples
#' preds <- data.frame(
#'   isolate_id = c("a", "a", "b"),
#'   species = "E. coli",
#'   drug = "ciprofloxacin",
#'   truth = c("R", "R", "S"),
#'   probability = c(0.9, 0.8, 1.4)
#' )
#' amr_validate(amr_import_predictions(preds))
amr_validate <- function(x, min_n = 10, min_resistant = 5) {
  if (!inherits(x, "amr_data")) {
    stop("`x` must be an amr_data object; use amr_import_predictions() first.", call. = FALSE)
  }
  d <- as.data.frame(x)
  issues <- list()
  add <- function(level, check, message, n_rows = NA_integer_) {
    issues[[length(issues) + 1]] <<- data.frame(
      level = level, check = check, message = message, n_rows = as.integer(n_rows),
      stringsAsFactors = FALSE
    )
  }

  # Required fields ----------------------------------------------------------
  for (col in c("isolate_id", "species", "drug")) {
    if (!col %in% names(d)) add("error", "required_column", paste0("Missing required column '", col, "'."))
  }
  if (!any(c("truth_sir", "truth_mic") %in% names(d))) {
    add("error", "required_column", "Supply a reference phenotype: 'truth_sir' (or 'truth') or 'truth_mic'.")
  }
  if (!any(c("probability", "prediction", "predicted_mic") %in% names(d))) {
    add("error", "required_column", "Supply predictions: 'probability', 'prediction' or 'predicted_mic'.")
  }
  if (any(vapply(issues, function(i) i$level == "error", logical(1)))) {
    return(new_validation(issues))
  }

  # Identifiers ----------------------------------------------------------------
  n <- sum(is.na(d$isolate_id) | d$isolate_id == "")
  if (n > 0) add("error", "isolate_id", "Missing or empty isolate identifiers.", n)
  co <- attr(x, "coercion")
  describe <- function(field) {
    if (is.null(co)) return("")
    b <- co[co$field == field & co$status != "exact", , drop = FALSE]
    if (nrow(b) == 0) return("")
    lab <- paste0("'", b$input, "'", ifelse(is.na(b$matched_name), "", paste0(" (closest: ", b$matched_name, ")")))
    paste0(": ", paste(utils::head(lab, 3), collapse = ", "), if (nrow(b) > 3) ", ...")
  }
  n <- sum(is.na(d$species) | as.character(d$species) == "UNKNOWN")
  if (n > 0) add("error", "species", paste0("Species not recognised exactly by the AMR package", describe("species"), "."), n)
  n <- sum(is.na(d$drug))
  if (n > 0) add("error", "drug", paste0("Antimicrobial not recognised exactly by the AMR package", describe("drug"), "."), n)
  if (!is.null(co) && any(co$status == "approximate")) {
    a <- co[co$status == "approximate", , drop = FALSE]
    add("warning", "approximate_match",
        paste0("Names matched approximately; check them: ",
               paste(utils::head(paste0("'", a$input, "' -> ", a$matched_name), 3), collapse = ", "),
               if (nrow(a) > 3) ", ..."), nrow(a))
  }

  # Reference phenotype ------------------------------------------------------
  truth_missing <- rep(TRUE, nrow(d))
  if ("truth_sir" %in% names(d)) truth_missing <- truth_missing & is.na(d$truth_sir)
  if ("truth_mic" %in% names(d)) truth_missing <- truth_missing & is.na(d$truth_mic)
  n <- sum(truth_missing)
  if (n > 0) add("warning", "truth_missing", "Rows without a reference phenotype will be excluded from evaluation.", n)

  if ("truth_sir" %in% names(d)) {
    n <- sum(as.character(d$truth_sir) %in% c("I", "SDD"))
    if (n > 0) {
      add("info", "intermediate",
          "Reference contains I or SDD results; binary metrics need an explicit policy (see `intermediate` in amr_metrics()).", n)
    }
    n <- sum(as.character(d$truth_sir) %in% "NI")
    if (n > 0) add("warning", "non_interpretable", "Reference contains NI (non-interpretable) results; they are excluded from evaluation.", n)
  }

  # Predictions --------------------------------------------------------------
  ptype <- d$prediction_type[1]
  if (is.na(ptype)) {
    add("error", "prediction_type", "Could not detect the prediction type; set `prediction_type` in amr_import_predictions().")
  }
  if ("probability" %in% names(d)) {
    p <- d$probability
    n <- sum(!is.na(p) & (!is.finite(p) | p < 0 | p > 1))
    if (n > 0) add("error", "probability_range", "Probabilities must lie in [0, 1].", n)
  }
  if (identical(ptype, "score")) {
    n <- sum(is.na(d$probability))
    if (n > 0) add("warning", "prediction_missing", "Rows without a probability will be excluded.", n)
  }
  if (identical(ptype, "mic")) {
    n <- sum(is.na(d$predicted_mic))
    if (n > 0) add("warning", "prediction_missing", "Rows without a predicted MIC will be excluded.", n)
    if (!"truth_mic" %in% names(d)) {
      add("warning", "essential_agreement", "MIC predictions without 'truth_mic': essential agreement cannot be computed.")
    }
  }
  if (ptype %in% c("binary", "sir")) {
    n <- sum(is.na(d$prediction))
    if (n > 0) add("warning", "prediction_missing", "Rows without a predicted class will be excluded.", n)
  }

  # Duplicates and conflicts -------------------------------------------------
  keys <- key_columns(d)
  dup <- duplicated(d[, keys, drop = FALSE]) | duplicated(d[, keys, drop = FALSE], fromLast = TRUE)
  if (any(dup)) {
    add("error", "duplicate_key",
        paste0("Duplicate prediction records for the key (", paste(keys, collapse = ", "), ")."),
        sum(dup))
  }
  if ("truth_sir" %in% names(d)) {
    ident <- intersect(c("isolate_id", "species", "drug", "dataset_id"), names(d))
    tr <- d[!is.na(d$truth_sir), c(ident, "truth_sir"), drop = FALSE]
    tr$truth_sir <- as.character(tr$truth_sir)
    tr <- unique(tr)
    conflict <- duplicated(tr[, ident, drop = FALSE])
    if (any(conflict)) {
      add("error", "conflicting_truth",
          "The same isolate and antimicrobial has different reference phenotypes across rows.",
          sum(conflict))
    }
  }

  # Splits -------------------------------------------------------------------
  if ("split" %in% names(d)) {
    bad <- !is.na(d$split) & !d$split %in% valid_splits
    if (any(bad)) {
      add("error", "split_label",
          paste0("Unknown split labels: ", paste(unique(d$split[bad]), collapse = ", "),
                 ". Use ", paste(valid_splits, collapse = ", "), "."), sum(bad))
    }
    by_model <- split(d[, c("isolate_id", "split")], d$model_id)
    n_overlap <- sum(vapply(by_model, function(m) {
      length(intersect(m$isolate_id[m$split %in% "train"], m$isolate_id[m$split %in% "test"]))
    }, integer(1)))
    if (n_overlap > 0) {
      add("error", "split_overlap", "Isolates appear in both the train and test partitions of the same model.", n_overlap)
    }
  }

  # Breakpoints --------------------------------------------------------------
  if ("breakpoint_standard" %in% names(d)) {
    bp <- unique(d[!is.na(d$breakpoint_standard), c("species", "drug", "breakpoint_standard")])
    mixed <- duplicated(bp[, c("species", "drug")])
    if (any(mixed)) {
      add("warning", "mixed_breakpoints",
          "More than one breakpoint standard or version within a species-drug stratum.", sum(mixed))
    }
  } else {
    add("info", "breakpoint_standard", "No 'breakpoint_standard' column; record the breakpoint authority and version for reproducibility.")
  }

  # Strata -------------------------------------------------------------------
  if ("truth_sir" %in% names(d)) {
    ev <- d[as.character(d$truth_sir) %in% c("S", "R"), , drop = FALSE]
    if (nrow(ev) > 0) {
      grp <- interaction(ev[, stratum_columns(ev)], drop = TRUE, sep = " | ")
      counts <- tapply(as.character(ev$truth_sir), grp, function(v) c(n = length(v), r = sum(v == "R")))
      n_vec <- vapply(counts, `[[`, numeric(1), "n")
      r_vec <- vapply(counts, `[[`, numeric(1), "r")
      single <- names(counts)[r_vec == 0 | r_vec == n_vec]
      if (length(single) > 0) {
        add("warning", "single_class",
            paste0("Strata with only one phenotype class (sensitivity or specificity undefined): ",
                   paste(utils::head(single, 3), collapse = "; "), if (length(single) > 3) "; ..."),
            length(single))
      }
      small <- names(counts)[n_vec < min_n]
      if (length(small) > 0) {
        add("warning", "small_stratum", paste0("Strata with fewer than ", min_n, " S/R isolates."), length(small))
      }
      few_r <- names(counts)[r_vec > 0 & r_vec < min_resistant]
      if (length(few_r) > 0) {
        add("warning", "few_resistant",
            paste0("Strata with fewer than ", min_resistant,
                   " resistant isolates; very major error rates will be imprecise."), length(few_r))
      }
    }
  }

  new_validation(issues)
}

new_validation <- function(issues) {
  if (length(issues) == 0) {
    tab <- data.frame(level = character(), check = character(), message = character(),
                      n_rows = integer(), stringsAsFactors = FALSE)
  } else {
    tab <- do.call(rbind, issues)
    tab <- tab[order(match(tab$level, c("error", "warning", "info"))), , drop = FALSE]
    rownames(tab) <- NULL
  }
  structure(list(ok = !any(tab$level == "error"), issues = tab), class = "amr_validation")
}

#' @export
print.amr_validation <- function(x, ...) {
  tab <- x$issues
  counts <- table(factor(tab$level, levels = c("error", "warning", "info")))
  cat("<amr_validation> ", if (x$ok) "OK" else "FAILED", " | ",
      counts[["error"]], " error(s), ", counts[["warning"]], " warning(s), ",
      counts[["info"]], " note(s)\n", sep = "")
  if (nrow(tab) > 0) {
    sym <- c(error = "x", warning = "!", info = "i")
    for (i in seq_len(nrow(tab))) {
      n <- if (is.na(tab$n_rows[i])) "" else paste0(" [", tab$n_rows[i], "]")
      cat(" ", sym[[tab$level[i]]], " ", tab$check[i], ": ", tab$message[i], n, "\n", sep = "")
    }
  }
  invisible(x)
}

#' Summarise an AMRvalid prediction table
#'
#' Counts predictions and reference phenotypes per model, dataset, species
#' and antimicrobial.
#'
#' @param x An `amr_data` object.
#' @return A data frame with one row per stratum.
#' @export
#' @examples
#' preds <- data.frame(
#'   isolate_id = c("a", "b", "c", "d"),
#'   species = "K. pneumoniae",
#'   drug = c("meropenem", "meropenem", "ciprofloxacin", "ciprofloxacin"),
#'   truth = c("R", "S", "I", "R"),
#'   prediction = c("R", "S", "R", "R")
#' )
#' amr_data_summary(amr_import_predictions(preds))
amr_data_summary <- function(x) {
  if (!inherits(x, "amr_data")) {
    stop("`x` must be an amr_data object; use amr_import_predictions() first.", call. = FALSE)
  }
  d <- as.data.frame(x)
  cols <- stratum_columns(d)
  d$.key <- interaction(d[, cols], drop = TRUE, lex.order = TRUE)
  truth <- if ("truth_sir" %in% names(d)) as.character(d$truth_sir) else rep(NA_character_, nrow(d))
  rows <- lapply(split(seq_len(nrow(d)), d$.key), function(idx) {
    t <- truth[idx]
    n_s <- sum(t %in% "S"); n_r <- sum(t %in% "R")
    data.frame(
      d[idx[1], cols, drop = FALSE],
      n = length(idx),
      n_isolates = length(unique(d$isolate_id[idx])),
      n_S = n_s,
      n_I_SDD = sum(t %in% c("I", "SDD")),
      n_R = n_r,
      n_missing_truth = sum(is.na(t)),
      prevalence_R = if (n_s + n_r > 0) n_r / (n_s + n_r) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$species_name <- AMR::mo_name(out$species, language = "en")
  out$drug_name <- AMR::ab_name(out$drug, language = "en")
  rownames(out) <- NULL
  out[, c(cols, "species_name", "drug_name", setdiff(names(out), c(cols, "species_name", "drug_name")))]
}
