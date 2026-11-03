#' Record the provenance of an evaluation
#'
#' Collects what was evaluated and how, and labels every item by where it
#' came from: **checked** (derived by AMRvalid from the data or results),
#' **supplied** (declared by you in `model_info` or `dataset_info`) or
#' **not assessable** (needed for reporting but neither available nor
#' declared). Each item is mapped to the TRIPOD+AI checklist item it
#' informs.
#'
#' @param x An `amr_data` object.
#' @param ... Result objects to include: any of `amr_validation`,
#'   `amr_split_check`, `amr_eval`, `amr_bootstrap`, `amr_paired_compare`,
#'   `amr_calibration` and `amr_distance_eval`.
#' @param model_info Named list describing the model(s), e.g.
#'   `list(description = "XGBoost on AMRFinderPlus gene presence",
#'   training_data = "BV-BRC isolates collected 2015-2020",
#'   threshold_selection = "Youden index on the validation split",
#'   prediction_method = "predict() on a saved model object",
#'   class_imbalance = "none", fairness = "...")`.
#' @param dataset_info Named list describing the evaluation data, e.g.
#'   `list(source = "NCBI Pathogen Detection AST browser", setting = "...",
#'   eligibility = "...", data_availability = "...", code_availability = "...")`.
#'
#' @return An object of class `amr_provenance`: a data frame with columns
#'   `section`, `item`, `value`, `status` and `tripod_item`.
#' @references
#' Collins GS, Moons KGM, Dhiman P, et al. (2024). TRIPOD+AI statement:
#' updated guidance for reporting clinical prediction models that use
#' regression or machine learning methods. *BMJ* 385:e078378.
#' @seealso [amr_report()]
#' @export
#' @examples
#' sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10,
#'                     external_lineages = 0, return_distances = FALSE, seed = 1)
#' x <- amr_import_predictions(sim$predictions)
#' e <- amr_species_eval(x, threshold = 0.5)
#' amr_provenance(x, e, model_info = list(threshold_selection = "fixed at 0.5"))
amr_provenance <- function(x, ..., model_info = list(), dataset_info = list()) {
  if (!inherits(x, "amr_data")) {
    stop("`x` must be an amr_data object; use amr_import_predictions() first.", call. = FALSE)
  }
  res <- collect_results(list(...))
  d <- as.data.frame(x)
  rows <- list()
  add <- function(section, item, value, status, tripod) {
    rows[[length(rows) + 1]] <<- data.frame(section = section, item = item, value = paste(value, collapse = "; "),
                                            status = status, tripod_item = tripod, stringsAsFactors = FALSE)
  }
  supplied <- function(section, item, key, info, tripod) {
    v <- info[[key]]
    if (is.null(v) || !nzchar(paste(v, collapse = ""))) {
      add(section, item, "Not declared", "not assessable", tripod)
    } else {
      add(section, item, v, "supplied", tripod)
    }
  }

  # Data -------------------------------------------------------------------
  supplied("Data", "Data source", "source", dataset_info, "5a")
  if ("collection_date" %in% names(d) && any(!is.na(d$collection_date))) {
    r <- range(as.Date(d$collection_date), na.rm = TRUE)
    add("Data", "Collection dates", paste(format(r), collapse = " to "), "checked", "5b")
  } else {
    supplied("Data", "Collection dates", "dates", dataset_info, "5b")
  }
  supplied("Data", "Setting", "setting", dataset_info, "6a")
  supplied("Data", "Eligibility criteria", "eligibility", dataset_info, "6b")
  sp <- unique(suppressMessages(AMR::mo_name(d$species, language = "en")))
  ab <- unique(suppressMessages(AMR::ab_name(d$drug, language = "en")))
  add("Data", "Species", sp[!is.na(sp)], "checked", "6b")
  add("Data", "Antimicrobials", ab[!is.na(ab)], "checked", "8a")
  add("Data", "Predictions, isolates, models, datasets",
      sprintf("%d predictions; %d isolates; %d model(s); %d dataset(s)", nrow(d),
              length(unique(d$isolate_id)), length(unique(d$model_id)), length(unique(d$dataset_id))),
      "checked", "20a")

  # Outcome ------------------------------------------------------------------
  if ("breakpoint_standard" %in% names(d) && any(!is.na(d$breakpoint_standard))) {
    add("Outcome", "Breakpoint standard", unique(stats::na.omit(d$breakpoint_standard)), "checked", "8a")
  } else {
    supplied("Outcome", "Breakpoint standard", "breakpoint_standard", dataset_info, "8a")
  }
  if ("ast_method" %in% names(d) && any(!is.na(d$ast_method))) {
    add("Outcome", "Reference AST method", unique(stats::na.omit(d$ast_method)), "checked", "8a")
  } else {
    supplied("Outcome", "Reference AST method", "ast_method", dataset_info, "8a")
  }
  if ("truth_sir" %in% names(d)) {
    tab <- table(factor(as.character(d$truth_sir), levels = c("S", "SDD", "I", "R", "NI")), useNA = "ifany")
    add("Outcome", "Reference phenotype counts (all rows)",
        paste(names(tab), tab, sep = " = ", collapse = ", "), "checked", "20a")
  }
  pol <- first_setting(res, "intermediate")
  add("Outcome", "Policy for I and SDD results", if (is.null(pol)) "Not recorded in supplied results" else pol,
      if (is.null(pol)) "not assessable" else "checked", "8a")

  # Missing data ---------------------------------------------------------------
  miss_truth <- if ("truth_sir" %in% names(d)) sum(is.na(d$truth_sir)) else NA
  pcol <- switch(d$prediction_type[1], score = "probability", mic = "predicted_mic", "prediction")
  miss_pred <- if (pcol %in% names(d)) sum(is.na(d[[pcol]])) else NA
  add("Missing data", "Rows without reference phenotype / prediction",
      sprintf("%s / %s (excluded from evaluation)", miss_truth, miss_pred), "checked", "11")

  # Validation design ----------------------------------------------------------
  if ("split" %in% names(d)) {
    add("Design", "Partitions present", paste(names(table(d$split)), table(d$split), sep = " = ", collapse = ", "),
        "checked", "12a")
  } else {
    add("Design", "Partitions present", "No split column; all rows treated as test", "checked", "12a")
  }
  sc <- res$amr_split_check
  if (!is.null(sc)) {
    cons <- sc$consistent_designs
    designs <- vapply(seq_len(nrow(cons)), function(i) {
      ok <- names(cons)[-1][vapply(cons[i, -1], isTRUE, logical(1))]
      paste0(cons$model_id[i], ": ", if (length(ok)) paste(ok, collapse = ", ") else "none")
    }, character(1))
    add("Design", "Separation of training and test isolates (consistent designs)", designs, "checked", "12a")
    if (!is.null(sc$claimed_design)) {
      add("Design", "Claimed validation design", paste0(sc$claimed_design, if (sc$ok) " (consistent)" else
        paste0(" (NOT consistent: ", paste(sc$violated, collapse = ", "), ")")), "checked", "12a")
    }
  } else {
    add("Design", "Separation of training and test isolates", "Not checked; run amr_validate_splits()",
        "not assessable", "12a")
  }
  supplied("Model", "Model description", "description", model_info, "12c")
  supplied("Model", "Training data", "training_data", model_info, "12a")
  supplied("Model", "How predictions were produced", "prediction_method", model_info, "12g")
  supplied("Model", "Class imbalance handling", "class_imbalance", model_info, "13")
  supplied("Model", "Fairness considerations", "fairness", model_info, "14")

  # Analysis -------------------------------------------------------------------
  thr <- first_setting(res, "threshold")
  add("Analysis", "Classification threshold", if (is.null(thr)) "None (threshold-free metrics only)" else thr,
      "checked", "12e")
  supplied("Analysis", "How the threshold was chosen", "threshold_selection", model_info, "12e")
  ci <- first_setting(res, "ci"); conf <- first_setting(res, "conf")
  if (!is.null(ci)) {
    add("Analysis", "Interval method for proportions",
        paste0(ci, " (", 100 * conf, "%)"), "checked", "12e")
  }
  vd <- first_setting(res, "vme_denominator")
  if (!is.null(vd)) {
    add("Analysis", "Very major / major error denominator",
        if (vd == "class") "reference-resistant / reference-susceptible isolates" else "all isolates",
        "checked", "12e")
  }
  tp <- first_setting(res, "target_prevalence")
  if (!is.null(tp)) add("Analysis", "Target prevalence for PPV/NPV", tp, "checked", "12e")
  wb <- first_setting(res, "weight_by")
  if (!is.null(wb)) add("Analysis", "Cluster weighting", paste0("each ", wb, " weighted equally"), "checked", "12d")
  measured <- c(if (!is.null(res$amr_eval)) "discrimination and classification metrics",
                if (!is.null(res$amr_calibration)) "calibration",
                if (!is.null(res$amr_distance_eval)) "performance by distance to training",
                if (!is.null(res$amr_paired_compare)) "paired model comparison")
  add("Analysis", "Performance measures assessed", if (length(measured)) measured else "None supplied",
      if (length(measured)) "checked" else "not assessable", "12e")
  if (!is.null(res$amr_bootstrap)) {
    s <- attr(res$amr_bootstrap, "settings")
    add("Analysis", "Group-aware uncertainty",
        sprintf("cluster bootstrap by %s, B = %d, seed = %s", s$group, s$B,
                if (is.null(s$seed)) "not set" else s$seed), "checked", "12d")
  } else {
    add("Analysis", "Group-aware uncertainty", "Not assessed; intervals assume independent isolates",
        "not assessable", "12d")
  }
  if (is.null(res$amr_calibration)) {
    if (identical(d$prediction_type[1], "score")) {
      add("Analysis", "Calibration", "Not assessed; run amr_calibration()", "not assessable", "12e")
    } else {
      add("Analysis", "Calibration", "Not applicable (predictions are categories, not probabilities)",
          "checked", "12e")
    }
  }

  # Availability -----------------------------------------------------------------
  supplied("Open science", "Data availability", "data_availability", dataset_info, "18e")
  supplied("Open science", "Code availability", "code_availability", model_info, "18f")
  add("Open science", "Software", sprintf("AMRvalid %s; AMR %s; %s", utils::packageVersion("AMRvalid"),
                                          utils::packageVersion("AMR"), R.version.string), "checked", "18f")
  add("Open science", "Date of analysis", format(Sys.Date()), "checked", "18f")

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("amr_provenance", "data.frame"))
}

#' @export
print.amr_provenance <- function(x, ...) {
  if (!all(c("section", "item", "value", "status", "tripod_item") %in% names(x))) {
    return(print(as.data.frame(unclass_df(x)), ...))
  }
  counts <- table(factor(x$status, levels = c("checked", "supplied", "not assessable")))
  cat("<amr_provenance> ", counts[["checked"]], " checked | ", counts[["supplied"]], " supplied | ",
      counts[["not assessable"]], " not assessable\n", sep = "")
  sym <- c(checked = "v", supplied = "*", `not assessable` = "?")
  for (sec in unique(x$section)) {
    cat("\n", sec, "\n", sep = "")
    s <- x[x$section == sec, ]
    for (i in seq_len(nrow(s))) {
      val <- if (nchar(s$value[i]) > 70) paste0(substr(s$value[i], 1, 67), "...") else s$value[i]
      cat(" ", sym[[s$status[i]]], " [", s$tripod_item[i], "] ", s$item[i], ": ", val, "\n", sep = "")
    }
  }
  cat("\nv checked by AMRvalid   * supplied by you   ? not assessable (needed for TRIPOD+AI)\n")
  invisible(x)
}

#' Write an evaluation report
#'
#' Writes a self-contained Markdown (or HTML) report: data summary,
#' validation checks, split integrity, performance per stratum with
#' intervals, calibration, distance-to-training results, uncertainty,
#' automatically detected limitations, and a TRIPOD+AI checklist table
#' showing which items the analysis covers.
#'
#' @inheritParams amr_provenance
#' @param file Output path. If `NULL`, the Markdown text is returned
#'   instead of written.
#' @param format `"markdown"` (default) or `"html"`. HTML needs the
#'   rmarkdown package and pandoc.
#' @param title Report title.
#'
#' @return The path to the report (invisibly), or the Markdown text when
#'   `file` is `NULL`.
#' @export
#' @examples
#' sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10,
#'                     external_lineages = 0, return_distances = FALSE, seed = 1)
#' x <- amr_import_predictions(sim$predictions)
#' e <- amr_species_eval(x, threshold = 0.5)
#' report <- amr_report(x, e, file = file.path(tempdir(), "amrvalid-report.md"))
#' cat(readLines(report, n = 20), sep = "\n")
amr_report <- function(x, ..., model_info = list(), dataset_info = list(), file = NULL,
                       format = c("markdown", "html"), title = "AMRvalid evaluation report") {
  format <- match.arg(format)
  objs <- list(...)
  res <- collect_results(objs)
  prov <- amr_provenance(x, ..., model_info = model_info, dataset_info = dataset_info)
  v <- res$amr_validation
  if (is.null(v)) v <- validate_amr_data(x)
  L <- character()
  h <- function(...) L <<- c(L, ...)
  h(paste0("# ", title), "",
    paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), " with AMRvalid ",
           utils::packageVersion("AMRvalid"), "."), "")

  h("## Data", "", md_table(amr_data_summary(x)[, c("model_id", "dataset_id", "species_name", "drug_name",
                                                       "n_isolates", "n_S", "n_I_SDD", "n_R", "prevalence_R")]), "")

  h("## Validation checks", "")
  if (nrow(v$issues) == 0) h("No issues found.", "") else h(md_table(v$issues), "")

  sc <- res$amr_split_check
  if (!is.null(sc)) {
    h("## Split integrity", "", md_table(sc$summary[, c("model_id", "check", "n_test_units", "n_shared",
                                                       "prop_test_isolates_affected", "note")]), "")
    if (length(sc$violated)) h(paste0("**Split is not consistent with:** ", paste(sc$violated, collapse = ", ")), "")
  }

  e <- res$amr_eval
  if (!is.null(e)) {
    key <- c("n", "prevalence_R", "sensitivity", "specificity", "ppv", "npv", "roc_auc", "brier",
             "very_major_error", "major_error", "categorical_agreement", "essential_agreement",
             "ppv_at_target_prevalence", "npv_at_target_prevalence")
    ed <- as.data.frame(e)
    ed <- ed[ed$metric %in% key, , drop = FALSE]
    ed$value <- fmt_ci(ed$estimate, ed$lower, ed$upper)
    h("## Performance by stratum", "", md_table(ed[, c("model_id", "dataset_id", "drug_name", "metric", "value")]), "")
  }
  cal <- res$amr_calibration
  if (!is.null(cal)) {
    s <- cal$summary
    s$intercept_ci <- fmt_ci(s$intercept, s$intercept_lower, s$intercept_upper)
    s$slope_ci <- fmt_ci(s$slope, s$slope_lower, s$slope_upper)
    h("## Calibration", "", md_table(s[, c("model_id", "dataset_id", "drug_name", "n", "oe_ratio",
                                           "intercept_ci", "slope_ci", "ece")]), "",
      "Ideal values: O/E = 1, intercept = 0, slope = 1.", "")
  }
  de <- res$amr_distance_eval
  if (!is.null(de)) {
    tr <- de$trend
    tr$value <- fmt_ci(tr$effect, tr$lower, tr$upper)
    h("## Performance by distance to training", "", md_table(tr[, c("model_id", "dataset_id", "drug_name", "model",
                                                                     "value", "n")]), "")
  }
  b <- res$amr_bootstrap
  if (!is.null(b)) {
    bd <- as.data.frame(b)
    bd$value <- fmt_ci(bd$estimate, bd$lower, bd$upper)
    h(paste0("## Group-aware intervals (resampling unit: ", attr(b, "settings")$group, ")"), "",
      md_table(bd[, c("model_id", "dataset_id", "drug_name", "metric", "value", "n_groups")]), "")
  }
  pc <- res$amr_paired_compare
  if (!is.null(pc)) {
    pd <- as.data.frame(pc)
    pd$difference_ci <- fmt_ci(pd$difference, pd$lower, pd$upper)
    h(paste0("## Paired comparison (", paste(rev(attr(pc, "settings")$models), collapse = " - "), ")"), "",
      md_table(pd[, c("dataset_id", "drug_name", "metric", "difference_ci", "n_matched")]), "")
  }

  lim <- detected_limitations(x, v, res)
  h("## Limitations detected", "")
  if (length(lim)) h(paste0("- ", lim), "") else h("None detected automatically.", "")
  h("AMRvalid checks the supplied data and metadata only. It cannot confirm that feature selection, tuning or",
    "threshold choice used training data alone.", "")

  h("## TRIPOD+AI coverage", "",
    "Status: **checked** = derived by AMRvalid; **supplied** = declared by the analyst; **not assessable** =",
    "needed for reporting but not available.", "",
    md_table(tripod_order(as.data.frame(prov))[, c("tripod_item", "section", "item", "value", "status")]), "",
    "Item numbers refer to the TRIPOD+AI checklist (Collins et al., BMJ 2024;385:e078378).", "")

  text <- paste(L, collapse = "\n")
  if (is.null(file)) return(text)
  if (format == "markdown") {
    writeLines(text, file)
  } else {
    if (!requireNamespace("rmarkdown", quietly = TRUE) || !rmarkdown::pandoc_available()) {
      stop("HTML output needs the rmarkdown package and pandoc; use format = \"markdown\".", call. = FALSE)
    }
    md <- tempfile(fileext = ".md")
    css <- tempfile(fileext = ".html")
    writeLines(sub("^# [^\n]*\n\n", "", text), md)
    writeLines(report_css, css)
    rmarkdown::pandoc_convert(md, to = "html5", output = normalizePath(file, mustWork = FALSE),
                              options = c("--standalone", "--metadata", paste0("title=", title),
                                          "--include-in-header", css))
  }
  invisible(file)
}

# Helpers -------------------------------------------------------------------

report_css <- "<style>
html { color: #2b2b2b; background: #ffffff; }
body { max-width: 1100px; padding: 24px 16px; font-family: system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif;
       font-size: 15px; line-height: 1.5; hyphens: none; }
h1 { font-size: 1.8em; margin-bottom: 0.2em; } h2 { font-size: 1.25em; margin-top: 2em;
     border-bottom: 1px solid #e6e6e3; padding-bottom: 0.2em; }
table { display: block; overflow-x: auto; border-collapse: collapse; font-size: 0.85em; margin: 0.5em 0 1em; }
th, td { padding: 4px 10px; border-bottom: 1px solid #e6e6e3; text-align: left; vertical-align: top; }
th { background: #f5f5f3; white-space: nowrap; }
td { max-width: 420px; }
</style>"

# Sort provenance rows in TRIPOD+AI checklist order (5a, 5b, ..., 12a, ..., 20a).
tripod_order <- function(p) {
  num <- as.numeric(sub("[a-z]+$", "", p$tripod_item))
  let <- sub("^[0-9]+", "", p$tripod_item)
  p[order(num, let, seq_len(nrow(p))), , drop = FALSE]
}

collect_results <- function(objs) {
  known <- c("amr_validation", "amr_split_check", "amr_eval", "amr_bootstrap", "amr_paired_compare",
             "amr_calibration", "amr_distance_eval")
  res <- list()
  for (o in objs) {
    cls <- intersect(class(o), known)
    if (length(cls) == 0) {
      stop("Unsupported object of class '", class(o)[1], "'. Supply AMRvalid result objects.", call. = FALSE)
    }
    res[[cls[1]]] <- o
  }
  res
}

first_setting <- function(res, name) {
  for (o in res) {
    s <- if (is.list(o) && !is.null(o$settings)) o$settings else attr(o, "settings")
    if (!is.null(s[[name]])) return(s[[name]])
  }
  NULL
}

fmt_ci <- function(est, lower, upper, digits = 3) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  whole <- !is.na(est) & is.na(lower) & abs(est - round(est)) < 1e-9 & abs(est) >= 1
  ifelse(is.na(est), "",
         ifelse(whole, format(round(est), trim = TRUE, scientific = FALSE),
                ifelse(is.na(lower), f(est), paste0(f(est), " [", f(lower), ", ", f(upper), "]"))))
}

md_table <- function(df, digits = 3) {
  if (nrow(df) == 0) return("(none)")
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      vapply(col, function(v) {
        if (is.na(v)) "" else if (v == round(v)) format(v, trim = TRUE, scientific = FALSE)
        else formatC(v, digits = digits, format = "f")
      }, character(1))
    } else {
      v <- as.character(col)
      v[is.na(v)] <- ""
      gsub("|", "\\|", gsub("\n", " ", v), fixed = TRUE)
    }
  })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

detected_limitations <- function(x, v, res) {
  out <- character()
  iss <- v$issues
  for (chk in c("few_resistant", "small_stratum", "single_class", "mixed_breakpoints", "approximate_match")) {
    if (chk %in% iss$check) out <- c(out, iss$message[iss$check == chk][1])
  }
  sc <- res$amr_split_check
  if (!is.null(sc) && length(sc$violated)) {
    out <- c(out, paste0("The split is not consistent with: ", paste(sc$violated, collapse = ", "), "."))
  }
  if (is.null(sc)) out <- c(out, "Overlap between training and test isolates was not checked.")
  b <- res$amr_bootstrap
  if (!is.null(b) && any(b$few_groups)) {
    out <- c(out, paste0("Fewer than ", attr(b, "settings")$min_groups,
                         " independent groups in some strata; cluster bootstrap intervals may be too narrow."))
  }
  if (!is.null(b) && any(b$zero_width)) out <- c(out, "Some bootstrap intervals have zero width.")
  if (is.null(b)) out <- c(out, "Intervals assume independent isolates; clonal relatedness was not accounted for.")
  cal <- res$amr_calibration
  if (!is.null(cal) && any(cal$summary$sparse)) out <- c(out, "Calibration estimates are unstable in sparse strata.")
  if (is.null(first_setting(res, "threshold")) && identical(x$prediction_type[1], "score")) {
    out <- c(out, "No classification threshold: only threshold-free metrics were reported.")
  }
  unique(out)
}
