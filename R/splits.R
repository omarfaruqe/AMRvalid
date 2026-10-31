#' Check the integrity of a train/test split
#'
#' Compares the isolates a model was trained on with the isolates it is
#' evaluated on, and reports which validation designs the split is
#' consistent with. AMRvalid cannot see upstream modelling steps, so a
#' clean result shows only that no visible overlap exists in the supplied
#' metadata.
#'
#' @param x An `amr_data` object. If it has a `split` column, rows with
#'   `split == "test"` are the evaluation set and rows with
#'   `split == "train"` are training isolates (unless `train_manifest` is
#'   given). Without a `split` column every row is treated as test.
#' @param train_manifest Optional data frame listing training isolates, with
#'   an `isolate_id` column and optionally `model_id` (if absent, the
#'   manifest applies to every model) and the grouping, date and country
#'   columns to compare.
#' @param group Grouping columns to check for overlap. Defaults to whichever
#'   of `lineage`, `genomic_cluster` and `sequence_type` are available.
#' @param claimed_design Optional design the evaluation claims to follow:
#'   `"isolate"`, `"group"` (held-out groups of the first `group` column),
#'   `"temporal"` or `"geographic"`. A split that contradicts the claim
#'   makes `ok` `FALSE`.
#'
#' @return An object of class `amr_split_check`: a list with `ok`,
#'   `summary` (one row per model and check) and `consistent_designs`
#'   (one row per model).
#' @export
#' @examples
#' preds <- data.frame(
#'   isolate_id = c("a", "b", "c", "d", "e", "f"),
#'   species = "K. pneumoniae",
#'   drug = "meropenem",
#'   truth = c("R", "S", "R", "S", "R", "S"),
#'   probability = c(0.9, 0.1, 0.8, 0.3, 0.7, 0.2),
#'   split = c("train", "train", "train", "test", "test", "test"),
#'   lineage = c("ST258", "ST11", "ST258", "ST15", "ST258", "ST307")
#' )
#' amr_validate_splits(amr_import_predictions(preds), claimed_design = "group")
amr_validate_splits <- function(x, train_manifest = NULL, group = NULL,
                                claimed_design = NULL) {
  if (!inherits(x, "amr_data")) {
    stop("`x` must be an amr_data object; use amr_import_predictions() first.", call. = FALSE)
  }
  d <- as.data.frame(x)
  if (!is.null(claimed_design)) {
    claimed_design <- match.arg(claimed_design, c("isolate", "group", "temporal", "geographic"))
  }

  if ("split" %in% names(d)) {
    test <- d[d$split %in% "test", , drop = FALSE]
  } else {
    test <- d
  }
  if (nrow(test) == 0) stop("No test rows found.", call. = FALSE)

  if (!is.null(train_manifest)) {
    if (!is.data.frame(train_manifest) || !"isolate_id" %in% names(train_manifest)) {
      stop("`train_manifest` must be a data frame with an 'isolate_id' column.", call. = FALSE)
    }
    train <- as.data.frame(train_manifest, stringsAsFactors = FALSE)
    train$isolate_id <- trimws(as.character(train$isolate_id))
  } else if ("split" %in% names(d) && any(d$split %in% "train")) {
    train <- d[d$split %in% "train", , drop = FALSE]
  } else {
    stop("No training isolates: supply `train_manifest` or rows with split = 'train'.", call. = FALSE)
  }

  if (is.null(group)) {
    group <- intersect(c("lineage", "genomic_cluster", "sequence_type"), names(test))
  }
  missing_group <- setdiff(group, names(test))
  if (length(missing_group) > 0) {
    stop("Grouping columns not found in `x`: ", paste(missing_group, collapse = ", "), call. = FALSE)
  }
  if (identical(claimed_design, "group") && length(group) == 0) {
    stop("`claimed_design = \"group\"` needs a grouping column (lineage, genomic_cluster or sequence_type).",
         call. = FALSE)
  }

  models <- unique(test$model_id)
  rows <- list()
  designs <- list()
  for (m in models) {
    te <- test[test$model_id == m, , drop = FALSE]
    tr <- if ("model_id" %in% names(train)) train[train$model_id == m, , drop = FALSE] else train
    te_iso <- unique(te$isolate_id)
    tr_iso <- unique(tr$isolate_id)
    shared_iso <- intersect(te_iso, tr_iso)
    rows[[length(rows) + 1]] <- split_row(m, "isolate", length(tr_iso), length(te_iso),
                                          length(shared_iso), length(shared_iso) / length(te_iso))
    cons <- list(model_id = m, isolate = length(shared_iso) == 0)

    for (g in group) {
      if (!g %in% names(tr)) {
        rows[[length(rows) + 1]] <- split_row(m, g, NA, NA, NA, NA,
                                              "Training metadata lacks this column; overlap not assessable.")
        cons[[g]] <- NA
        next
      }
      te_g <- unique(stats::na.omit(as.character(te[[g]])))
      tr_g <- unique(stats::na.omit(as.character(tr[[g]])))
      shared <- intersect(te_g, tr_g)
      frac_iso <- mean(as.character(te[[g]][!duplicated(te$isolate_id)]) %in% shared)
      note <- if (any(is.na(te[[g]]))) "Some test isolates have no group label." else ""
      rows[[length(rows) + 1]] <- split_row(m, g, length(tr_g), length(te_g), length(shared),
                                            frac_iso, note)
      cons[[g]] <- length(shared) == 0
    }

    if ("collection_date" %in% names(te) && "collection_date" %in% names(tr)) {
      tr_max <- suppressWarnings(max(as.Date(tr$collection_date), na.rm = TRUE))
      te_dates <- as.Date(te$collection_date[!duplicated(te$isolate_id)])
      early <- sum(te_dates <= tr_max, na.rm = TRUE)
      rows[[length(rows) + 1]] <- split_row(
        m, "collection_date", NA, sum(!is.na(te_dates)), early,
        if (any(!is.na(te_dates))) early / sum(!is.na(te_dates)) else NA,
        if (is.finite(tr_max)) paste0("Latest training date: ", format(tr_max), ".") else "No training dates.")
      cons$temporal <- is.finite(tr_max) && early == 0 && all(!is.na(te_dates))
    } else {
      cons$temporal <- NA
    }

    if ("country" %in% names(te) && "country" %in% names(tr)) {
      te_c <- unique(stats::na.omit(as.character(te$country)))
      tr_c <- unique(stats::na.omit(as.character(tr$country)))
      shared_c <- intersect(te_c, tr_c)
      frac <- mean(as.character(te$country[!duplicated(te$isolate_id)]) %in% shared_c)
      rows[[length(rows) + 1]] <- split_row(m, "country", length(tr_c), length(te_c), length(shared_c), frac)
      cons$geographic <- length(shared_c) == 0
    } else {
      cons$geographic <- NA
    }
    designs[[length(designs) + 1]] <- as.data.frame(cons, stringsAsFactors = FALSE)
  }

  summary <- do.call(rbind, rows)
  consistent <- do.call(rbind, lapply(designs, function(z) {
    z[setdiff(c("model_id", "isolate", group, "temporal", "geographic"), names(z))] <- NA
    z[, c("model_id", "isolate", group, "temporal", "geographic"), drop = FALSE]
  }))
  rownames(summary) <- NULL
  rownames(consistent) <- NULL

  ok <- all(consistent$isolate)
  violated <- character()
  if (!all(consistent$isolate)) violated <- c(violated, "isolate")
  if (!is.null(claimed_design)) {
    col <- switch(claimed_design, isolate = "isolate", group = group[1],
                  temporal = "temporal", geographic = "geographic")
    claim_ok <- consistent[[col]]
    if (any(is.na(claim_ok))) {
      ok <- FALSE
      violated <- c(violated, paste0(claimed_design, " (not assessable)"))
    } else if (!all(claim_ok)) {
      ok <- FALSE
      violated <- c(violated, claimed_design)
    }
  }
  structure(list(ok = ok, summary = summary, consistent_designs = consistent,
                 claimed_design = claimed_design, group = group, violated = unique(violated)),
            class = "amr_split_check")
}

split_row <- function(model_id, check, n_train_units, n_test_units, n_shared, prop_test_isolates_affected,
                      note = "") {
  data.frame(model_id = model_id, check = check, n_train_units = n_train_units,
             n_test_units = n_test_units, n_shared = n_shared,
             prop_test_isolates_affected = prop_test_isolates_affected, note = note,
             stringsAsFactors = FALSE)
}

#' @export
print.amr_split_check <- function(x, ...) {
  cat("<amr_split_check> ", if (x$ok) "OK" else "PROBLEMS FOUND",
      if (!is.null(x$claimed_design)) paste0(" | claimed design: ", x$claimed_design), "\n", sep = "")
  if (length(x$violated) > 0) {
    cat(" x Split is not consistent with: ", paste(x$violated, collapse = ", "), "\n", sep = "")
  }
  s <- x$summary
  for (i in seq_len(nrow(s))) {
    label <- if (s$check[i] == "collection_date") {
      paste0(s$n_shared[i], " test isolate(s) collected on or before the latest training date")
    } else if (is.na(s$n_shared[i])) {
      "not assessable"
    } else {
      paste0(s$n_shared[i], " of ", s$n_test_units[i], " test ", s$check[i],
             if (s$check[i] == "isolate") "s" else " value(s)", " also in training")
    }
    pct <- if (is.na(s$prop_test_isolates_affected[i])) "" else
      sprintf(" (%.0f%% of test isolates)", 100 * s$prop_test_isolates_affected[i])
    flag <- if (!is.na(s$n_shared[i]) && s$n_shared[i] > 0) "!" else if (is.na(s$n_shared[i])) "?" else "v"
    cat(" ", flag, " [", s$model_id[i], "] ", s$check[i], ": ", label, pct,
        if (nzchar(s$note[i])) paste0(". ", s$note[i]), "\n", sep = "")
  }
  cat(" i AMRvalid checks supplied metadata only; it cannot certify that upstream feature selection\n",
      "   or tuning used training data alone.\n", sep = "")
  invisible(x)
}
