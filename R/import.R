#' Import predictions into the AMRvalid data contract
#'
#' Renames columns, coerces identifiers to AMR package classes, converts
#' predictions to a common representation and detects the prediction type.
#' The function never drops rows; problems are reported by
#' [validate_amr_data()].
#'
#' @param data A data frame with one row per prediction.
#' @param mapping Optional named character vector mapping contract field
#'   names to column names in `data`, e.g.
#'   `c(isolate_id = "Sample", truth_sir = "phenotype")`.
#' @param prediction_type One of `"binary"`, `"score"`, `"sir"` or `"mic"`.
#'   If `NULL` (default) it is detected from the columns present.
#' @param source Short label for where the predictions came from (recorded
#'   in the `source` attribute).
#' @param strict If `TRUE` (default), species and antimicrobial names must
#'   match an AMR package name, code, abbreviation or synonym exactly
#'   (ignoring case and punctuation). Approximate matches are set to `NA`,
#'   which [validate_amr_data()] reports as errors. The AMR package matches
#'   names fuzzily, so without this check a typo can silently become a
#'   different antimicrobial. If `FALSE`, approximate matches are kept and
#'   reported as warnings.
#'
#' @details
#' A column called `truth` is treated as `truth_sir`. Predicted classes may
#' be given as S/I/R text, `0`/`1` or `FALSE`/`TRUE` (where 1 or `TRUE`
#' means resistant). Missing `model_id` and `dataset_id` columns are filled
#' with `"model_1"` and `"dataset_1"`.
#'
#' The `coercion` attribute of the result lists every distinct input name,
#' the code it was matched to and whether the match was exact.
#'
#' Detection order for `prediction_type`: `predicted_mic` gives `"mic"`;
#' otherwise `probability` gives `"score"`; otherwise a `prediction` column
#' containing I or SDD gives `"sir"`; otherwise `"binary"`.
#'
#' @return A data frame of class `amr_data`.
#' @export
#' @examples
#' preds <- data.frame(
#'   isolate_id = c("a", "b", "c"),
#'   species = "K. pneumoniae",
#'   drug = "meropenem",
#'   truth = c("R", "S", "S"),
#'   probability = c(0.9, 0.2, 0.6)
#' )
#' amr_import_predictions(preds)
amr_import_predictions <- function(data, mapping = NULL, prediction_type = NULL,
                                   source = "custom", strict = TRUE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  x <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!is.null(mapping)) {
    if (is.null(names(mapping)) || any(names(mapping) == "")) {
      stop("`mapping` must be a named character vector: c(field = \"column\").", call. = FALSE)
    }
    missing_cols <- setdiff(unname(mapping), names(x))
    if (length(missing_cols) > 0) {
      stop("Columns named in `mapping` not found in `data`: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    names(x)[match(unname(mapping), names(x))] <- names(mapping)
  }
  if ("truth" %in% names(x) && !"truth_sir" %in% names(x)) {
    names(x)[names(x) == "truth"] <- "truth_sir"
  }

  if ("isolate_id" %in% names(x)) {
    x$isolate_id <- trimws(as.character(x$isolate_id))
  }
  coercion <- list()
  if ("species" %in% names(x)) {
    res <- coerce_identifier(x$species, "species", strict)
    x$species <- res$value
    coercion$species <- res$table
  }
  if ("drug" %in% names(x)) {
    res <- coerce_identifier(x$drug, "drug", strict)
    x$drug <- res$value
    coercion$drug <- res$table
  }
  if ("truth_sir" %in% names(x)) {
    x$truth_sir <- to_sir(x$truth_sir)
  }
  if ("truth_mic" %in% names(x)) {
    x$truth_mic <- suppressWarnings(AMR::as.mic(x$truth_mic))
  }
  if ("prediction" %in% names(x)) {
    x$prediction <- to_sir(x$prediction)
  }
  if ("predicted_mic" %in% names(x)) {
    x$predicted_mic <- suppressWarnings(AMR::as.mic(x$predicted_mic))
  }
  if ("probability" %in% names(x)) {
    x$probability <- suppressWarnings(as.numeric(x$probability))
  }
  if ("split" %in% names(x)) {
    x$split <- tolower(trimws(as.character(x$split)))
  }
  if ("collection_date" %in% names(x) && !inherits(x$collection_date, "Date")) {
    x$collection_date <- suppressWarnings(as.Date(as.character(x$collection_date)))
  }
  for (col in c("fold_id", "lineage", "sequence_type", "genomic_cluster", "country",
                "breakpoint_standard", "ast_method")) {
    if (col %in% names(x)) x[[col]] <- as.character(x[[col]])
  }
  if (!"model_id" %in% names(x)) x$model_id <- "model_1"
  if (!"dataset_id" %in% names(x)) x$dataset_id <- "dataset_1"
  x$model_id <- as.character(x$model_id)
  x$dataset_id <- as.character(x$dataset_id)

  if (is.null(prediction_type)) {
    prediction_type <- detect_prediction_type(x)
  } else {
    prediction_type <- match.arg(prediction_type, valid_prediction_types)
  }
  x$prediction_type <- prediction_type

  front <- intersect(amr_contract()$field, names(x))
  x <- x[, c(front, setdiff(names(x), front)), drop = FALSE]
  rownames(x) <- NULL

  structure(x,
            class = c("amr_data", "data.frame"),
            source = source,
            coercion = do.call(rbind, unname(coercion)),
            strict = strict,
            imported_at = Sys.time(),
            amrvalid_version = as.character(utils::packageVersion("AMRvalid")))
}

detect_prediction_type <- function(x) {
  if ("predicted_mic" %in% names(x) && any(!is.na(x$predicted_mic))) return("mic")
  if ("probability" %in% names(x) && any(!is.na(x$probability))) return("score")
  if ("prediction" %in% names(x)) {
    p <- as.character(x$prediction)
    if (any(p %in% c("I", "SDD"))) return("sir")
    return("binary")
  }
  NA_character_
}

# Convert S/I/R text, 0/1 or logical values to the AMR 'sir' class.
to_sir <- function(v) {
  if (is.logical(v)) {
    v <- ifelse(v, "R", "S")
  } else if (is.numeric(v)) {
    if (!all(v %in% c(0, 1, NA))) {
      stop("Numeric phenotype or prediction columns must contain only 0, 1 or NA.",
           call. = FALSE)
    }
    v <- ifelse(v == 1, "R", "S")
  } else {
    v <- trimws(as.character(v))
    v[v == ""] <- NA
  }
  suppressWarnings(suppressMessages(AMR::as.sir(v)))
}

#' @export
print.amr_data <- function(x, ...) {
  if (!all(c("isolate_id", "model_id", "species", "drug", "prediction_type") %in% names(x))) {
    return(print(as.data.frame(unclass_df(x)), ...))
  }
  n_models <- length(unique(x$model_id))
  n_strata <- nrow(unique(as.data.frame(x)[, c("species", "drug"), drop = FALSE]))
  cat("<amr_data> ", nrow(x), " predictions | ",
      length(unique(x$isolate_id)), " isolates | ",
      n_strata, " species-drug strata | ",
      n_models, " model(s) | type: ", x$prediction_type[1], "\n", sep = "")
  co <- attr(x, "coercion")
  if (!is.null(co)) {
    bad <- co[co$status != "exact", , drop = FALSE]
    for (i in seq_len(min(nrow(bad), 5))) {
      cat("! ", bad$field[i], " '", bad$input[i], "' ",
          if (bad$status[i] == "approximate") paste0("approximately matched to ", bad$matched_name[i])
          else if (bad$status[i] == "rejected") paste0("rejected (closest: ", bad$matched_name[i], ")")
          else "not recognised", "\n", sep = "")
    }
    if (nrow(bad) > 5) cat("! ... and", nrow(bad) - 5, "more; see attr(x, \"coercion\")\n")
  }
  print(utils::head(as.data.frame(x), 6), ...)
  if (nrow(x) > 6) cat("# ... with", nrow(x) - 6, "more rows\n")
  invisible(x)
}

# Normalise names for exact comparison: lower case, punctuation to spaces.
norm_name <- function(v) {
  v <- tolower(as.character(v))
  v <- gsub("[^a-z0-9]+", " ", v)
  trimws(v)
}

# Coerce species or drug labels with the AMR package, then check whether each
# distinct input matched exactly (code, name, abbreviation or synonym).
coerce_identifier <- function(raw, field, strict) {
  raw <- as.character(raw)
  inputs <- unique(raw[!is.na(raw) & trimws(raw) != ""])
  if (field == "species") {
    codes <- suppressWarnings(suppressMessages(AMR::as.mo(inputs, info = FALSE)))
    codes_chr <- as.character(codes)
    codes_chr[codes_chr == "UNKNOWN"] <- NA
    known <- !is.na(codes_chr)
    names_full <- rep(NA_character_, length(inputs))
    names_full[known] <- suppressMessages(AMR::mo_name(codes_chr[known], language = "en"))
    exact <- rep(FALSE, length(inputs))
    if (any(known)) {
      short <- suppressMessages(AMR::mo_shortname(codes_chr[known], language = "en"))
      cand <- cbind(norm_name(codes_chr[known]), norm_name(names_full[known]), norm_name(short))
      exact[known] <- norm_name(inputs[known]) == cand[, 1] |
        norm_name(inputs[known]) == cand[, 2] | norm_name(inputs[known]) == cand[, 3]
    }
  } else {
    codes <- suppressWarnings(suppressMessages(AMR::as.ab(inputs, info = FALSE)))
    codes_chr <- as.character(codes)
    known <- !is.na(codes_chr)
    names_full <- rep(NA_character_, length(inputs))
    names_full[known] <- suppressMessages(AMR::ab_name(codes_chr[known], language = "en"))
    ref <- AMR::antimicrobials
    exact <- rep(FALSE, length(inputs))
    for (i in which(known)) {
      row <- ref[as.character(ref$ab) == codes_chr[i], , drop = FALSE]
      cand <- c(as.character(row$ab), row$name, as.character(row$cid),
                unlist(row$atc), unlist(row$abbreviations), unlist(row$synonyms))
      exact[i] <- norm_name(inputs[i]) %in% norm_name(cand[!is.na(cand)])
    }
  }
  status <- ifelse(!known, "unknown", ifelse(exact, "exact", if (strict) "rejected" else "approximate"))
  keep <- known & (exact | !strict)
  lookup <- ifelse(keep, codes_chr, NA_character_)
  value_chr <- unname(lookup[match(raw, inputs)])
  value <- if (field == "species") AMR::as.mo(value_chr, info = FALSE) else AMR::as.ab(value_chr, info = FALSE)
  table <- data.frame(field = rep(field, length(inputs)), input = inputs,
                      matched_code = codes_chr, matched_name = names_full,
                      status = status, stringsAsFactors = FALSE)
  list(value = value, table = table)
}
