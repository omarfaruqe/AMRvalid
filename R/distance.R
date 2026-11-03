#' Nearest-training distance for each test isolate
#'
#' For every test isolate, finds the genomically closest training isolate
#' using distances computed upstream (for example Mash, SKA, core-genome
#' SNP or cgMLST allele distances). AMRvalid does not compute distances.
#'
#' @param test_ids Character vector of test isolate identifiers.
#' @param train_ids Character vector of training isolate identifiers.
#' @param distances Either a square numeric matrix with isolate identifiers
#'   as row and column names, or a data frame with columns `isolate_1`,
#'   `isolate_2` and `distance` (pairs may appear in either order).
#'
#' @return A data frame with `isolate_id`, `nearest_train_id` and
#'   `distance`. Test isolates with no distance to any training isolate get
#'   `NA`. A test isolate that is itself in the training set gets distance 0.
#' @export
#' @examples
#' m <- matrix(c(0, 2, 9,
#'               2, 0, 7,
#'               9, 7, 0), 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
#' amr_nearest_train_distance(c("b", "c"), "a", m)
amr_nearest_train_distance <- function(test_ids, train_ids, distances) {
  test_ids <- unique(as.character(test_ids))
  train_ids <- unique(as.character(train_ids))
  out <- data.frame(isolate_id = test_ids, nearest_train_id = NA_character_,
                    distance = NA_real_, stringsAsFactors = FALSE)
  if (is.matrix(distances)) {
    if (is.null(rownames(distances)) || is.null(colnames(distances))) {
      stop("A distance matrix needs isolate identifiers as row and column names.", call. = FALSE)
    }
    rows <- intersect(test_ids, rownames(distances))
    cols <- intersect(train_ids, colnames(distances))
    if (length(rows) > 0 && length(cols) > 0) {
      sub <- distances[rows, cols, drop = FALSE]
      best <- apply(sub, 1, function(v) if (all(is.na(v))) NA_integer_ else which.min(v))
      idx <- match(rows, out$isolate_id)
      out$nearest_train_id[idx] <- ifelse(is.na(best), NA, cols[best])
      out$distance[idx] <- vapply(seq_along(rows), function(i) {
        if (is.na(best[i])) NA_real_ else sub[i, best[i]]
      }, numeric(1))
    }
  } else if (is.data.frame(distances)) {
    need <- c("isolate_1", "isolate_2", "distance")
    if (!all(need %in% names(distances))) {
      stop("A distance table needs columns isolate_1, isolate_2 and distance.", call. = FALSE)
    }
    a <- as.character(distances$isolate_1)
    b <- as.character(distances$isolate_2)
    dv <- as.numeric(distances$distance)
    long <- data.frame(test = c(a, b), train = c(b, a), distance = c(dv, dv),
                       stringsAsFactors = FALSE)
    long <- long[long$test %in% test_ids & long$train %in% train_ids & !is.na(long$distance), ]
    if (nrow(long) > 0) {
      long <- long[order(long$test, long$distance), ]
      long <- long[!duplicated(long$test), ]
      idx <- match(long$test, out$isolate_id)
      out$nearest_train_id[idx] <- long$train
      out$distance[idx] <- long$distance
    }
  } else {
    stop("`distances` must be a matrix or a data frame.", call. = FALSE)
  }
  self <- out$isolate_id %in% train_ids
  out$nearest_train_id[self] <- out$isolate_id[self]
  out$distance[self] <- 0
  if (any(out$distance < 0, na.rm = TRUE)) stop("Distances must be non-negative.", call. = FALSE)
  out
}

#' Evaluate performance as a function of distance to the training set
#'
#' Assigns each test isolate the distance to its nearest training isolate,
#' groups test isolates into distance bands and evaluates each band per
#' species-antimicrobial stratum. A trend model then summarises how fast
#' performance decays as isolates become less similar to anything the model
#' was trained on.
#'
#' @inheritParams amr_validate_splits
#' @inheritParams amr_evaluate
#' @param distances Pairwise distances: a matrix or long table, see
#'   [amr_nearest_train_distance()]. Alternatively a data frame with
#'   `isolate_id` and `nearest_distance` (and optionally `model_id`) if
#'   nearest-training distances were computed upstream.
#' @param bands Either a single number of quantile-based bands (default
#'   `3`) or a numeric vector of band boundaries, e.g. `c(0, 10, 100, Inf)`
#'   for SNP distances.
#'
#' @details
#' The trend summary depends on what the predictions allow:
#'
#' * with categories (binary, S/I/R, or scores with a `threshold`): a
#'   logistic regression of correct classification on distance, reported
#'   as the odds ratio per standard deviation of distance;
#' * with scores and no threshold: a linear regression of the squared error
#'   (the per-isolate Brier score) on distance, reported as the change per
#'   standard deviation.
#'
#' For score predictions, each band also reports the calibration intercept
#' and slope (see [amr_calibration()]), showing whether a model becomes
#' over-confident on isolates unlike its training data.
#'
#' An odds ratio below 1 (or a positive Brier slope) means performance
#' worsens with distance. Intervals are Wald intervals that treat isolates
#' as independent; use them as a description, not a test.
#'
#' @return An object of class `amr_distance_eval`: a list with `nearest`
#'   (per-isolate distances and bands), `by_band` (metrics per band, in the
#'   format of [amr_evaluate()]) and `trend` (one row per stratum).
#' @export
#' @examples
#' set.seed(1)
#' ids <- sprintf("iso%02d", 1:40)
#' preds <- data.frame(
#'   isolate_id = ids,
#'   species = "E. coli",
#'   drug = "ciprofloxacin",
#'   truth = rep(c("R", "S"), 20),
#'   probability = runif(40),
#'   split = rep(c("train", "test"), each = 20)
#' )
#' dist <- as.matrix(stats::dist(seq_along(ids)))
#' dimnames(dist) <- list(ids, ids)
#' amr_distance_eval(amr_import_predictions(preds), dist, threshold = 0.5)
amr_distance_eval <- function(x, distances, train_manifest = NULL, bands = 3,
                              threshold = NULL,
                              intermediate = c("exclude", "susceptible", "resistant"),
                              ci = c("wilson", "clopper-pearson"), conf = 0.95) {
  check_amr_data(x)
  intermediate <- match.arg(intermediate)
  ci <- match.arg(ci)
  d <- as.data.frame(x)
  ptype <- d$prediction_type[1]
  test_rows <- if ("split" %in% names(d)) d$split %in% "test" else rep(TRUE, nrow(d))
  if (!any(test_rows)) stop("No test rows found.", call. = FALSE)

  precomputed <- is.data.frame(distances) && all(c("isolate_id", "nearest_distance") %in% names(distances))
  near <- list()
  for (m in unique(d$model_id[test_rows])) {
    te_ids <- unique(d$isolate_id[test_rows & d$model_id == m])
    if (precomputed) {
      pc <- distances
      if ("model_id" %in% names(pc)) pc <- pc[pc$model_id == m, , drop = FALSE]
      nd <- data.frame(isolate_id = te_ids,
                       nearest_train_id = NA_character_,
                       distance = as.numeric(pc$nearest_distance[match(te_ids, as.character(pc$isolate_id))]),
                       stringsAsFactors = FALSE)
    } else {
      tr_ids <- training_ids(d, train_manifest, m)
      nd <- amr_nearest_train_distance(te_ids, tr_ids, distances)
    }
    nd$model_id <- m
    near[[length(near) + 1]] <- nd
  }
  nearest <- do.call(rbind, near)
  if (all(is.na(nearest$distance))) {
    stop("No test isolate has a distance to a training isolate; check identifiers.", call. = FALSE)
  }
  n_missing <- sum(is.na(nearest$distance))
  if (n_missing > 0) {
    warning(n_missing, " test isolate(s) have no distance to any training isolate and are excluded.",
            call. = FALSE)
  }

  breaks <- band_breaks(nearest$distance, bands)
  nearest$band <- cut(nearest$distance, breaks = breaks, include.lowest = TRUE, right = FALSE,
                      dig.lab = 4)
  nearest <- nearest[, c("model_id", "isolate_id", "nearest_train_id", "distance", "band")]

  te <- x[test_rows, , drop = FALSE]
  key <- paste(te$model_id, te$isolate_id, sep = "\r")
  nkey <- paste(nearest$model_id, nearest$isolate_id, sep = "\r")
  te$.distance <- nearest$distance[match(key, nkey)]
  te$.band <- nearest$band[match(key, nkey)]
  te <- te[!is.na(te$.band), , drop = FALSE]
  if ("split" %in% names(te)) te$split <- NULL

  by_band <- list()
  for (b in levels(nearest$band)) {
    sub <- te[te$.band == b, , drop = FALSE]
    if (nrow(sub) == 0) next
    e <- suppressMessages(amr_evaluate(sub, threshold = threshold, intermediate = intermediate,
                                           ci = ci, conf = conf, split = NULL, validate = FALSE))
    e <- as.data.frame(e)
    if (identical(ptype, "score")) {
      cal <- calibration_band_rows(sub, intermediate, conf)
      if (!is.null(cal)) e <- rbind(e, cal[, names(e)])
    }
    e$band <- b
    by_band[[length(by_band) + 1]] <- e
  }
  by_band <- do.call(rbind, by_band)
  by_band$band <- factor(by_band$band, levels = levels(nearest$band))
  front <- c(stratum_columns(by_band), "band")
  by_band <- by_band[, c(front, setdiff(names(by_band), front))]
  rownames(by_band) <- NULL

  trend <- distance_trend(te, threshold, intermediate, conf)

  structure(list(nearest = nearest, by_band = by_band, trend = trend,
                 settings = list(bands = breaks, threshold = threshold, intermediate = intermediate,
                                 ci = ci, conf = conf, prediction_type = ptype)),
            class = "amr_distance_eval")
}

# Training isolate identifiers for one model.
training_ids <- function(d, train_manifest, model) {
  if (!is.null(train_manifest)) {
    if (!is.data.frame(train_manifest) || !"isolate_id" %in% names(train_manifest)) {
      stop("`train_manifest` must be a data frame with an 'isolate_id' column.", call. = FALSE)
    }
    tm <- train_manifest
    if ("model_id" %in% names(tm)) tm <- tm[tm$model_id == model, , drop = FALSE]
    return(unique(trimws(as.character(tm$isolate_id))))
  }
  if ("split" %in% names(d) && any(d$split %in% "train" & d$model_id == model)) {
    return(unique(d$isolate_id[d$split %in% "train" & d$model_id == model]))
  }
  stop("No training isolates for model '", model,
       "': supply `train_manifest` or rows with split = 'train'.", call. = FALSE)
}

band_breaks <- function(dist, bands) {
  dist <- dist[!is.na(dist)]
  if (length(bands) == 1) {
    if (bands < 1 || bands != round(bands)) stop("`bands` must be a positive integer or a vector of boundaries.",
                                                 call. = FALSE)
    br <- unique(stats::quantile(dist, probs = seq(0, 1, length.out = bands + 1), names = FALSE, type = 7))
    if (length(br) < 2) br <- c(min(dist), min(dist) + 1)
    br[length(br)] <- br[length(br)] + 1e-9 * max(1, abs(br[length(br)]))
    return(br)
  }
  br <- sort(unique(as.numeric(bands)))
  if (length(br) < 2) stop("`bands` needs at least two boundaries.", call. = FALSE)
  if (min(dist) < br[1] || max(dist) >= br[length(br)] && is.finite(br[length(br)])) {
    stop("Band boundaries must cover all distances (", min(dist), " to ", max(dist),
         "); use Inf as the last boundary.", call. = FALSE)
  }
  br
}

distance_trend <- function(te, threshold, intermediate, conf) {
  ptype <- te$prediction_type[1]
  cols <- stratum_columns(te)
  key <- interaction(as.data.frame(te)[, cols], drop = TRUE, lex.order = TRUE, sep = "\r")
  z <- stats::qnorm(1 - (1 - conf) / 2)
  rows <- lapply(split(seq_len(nrow(te)), key), function(idx) {
    sub <- te[idx, , drop = FALSE]
    first <- as.data.frame(sub)[1, cols, drop = FALSE]
    first[] <- lapply(first, as.character)
    truth <- collapse_sir(sub$truth_sir, intermediate)
    dist <- sub$.distance
    sdv <- stats::sd(dist[!is.na(truth)])
    res <- data.frame(first, model = NA_character_, effect = NA_real_, lower = NA_real_,
                      upper = NA_real_, n = NA_integer_, distance_sd = sdv,
                      note = "", stringsAsFactors = FALSE)
    if (identical(ptype, "mic")) {
      res$note <- "Trend not implemented for MIC predictions."
      return(res)
    }
    pred <- predicted_class(sub, threshold, intermediate)
    if (!is.null(pred)) {
      ok <- !is.na(truth) & !is.na(pred)
      y <- as.integer(truth[ok] == pred[ok]); dz <- dist[ok] / sdv
      res$model <- "logistic: correct ~ distance (odds ratio per SD)"
      res$n <- sum(ok)
      fit <- trend_fit(y, dz, binomial = TRUE)
      if (!is.null(fit)) {
        res$effect <- exp(fit[1]); res$lower <- exp(fit[1] - z * fit[2]); res$upper <- exp(fit[1] + z * fit[2])
      } else {
        res$note <- "Too little variation in distance or outcome to fit a trend."
      }
    } else {
      ok <- !is.na(truth) & !is.na(sub$probability)
      y <- (sub$probability[ok] - as.integer(truth[ok] == "R"))^2; dz <- dist[ok] / sdv
      res$model <- "linear: squared error ~ distance (change per SD)"
      res$n <- sum(ok)
      fit <- trend_fit(y, dz, binomial = FALSE)
      if (!is.null(fit)) {
        res$effect <- fit[1]; res$lower <- fit[1] - z * fit[2]; res$upper <- fit[1] + z * fit[2]
      } else {
        res$note <- "Too little variation in distance or outcome to fit a trend."
      }
    }
    res
  })
  out <- do.call(rbind, rows)
  out$drug_name <- suppressMessages(AMR::ab_name(out$drug, language = "en"))
  rownames(out) <- NULL
  out
}

# Returns c(slope, se) or NULL when the model cannot be estimated.
trend_fit <- function(y, dz, binomial) {
  if (length(y) < 5 || !is.finite(stats::sd(dz)) || stats::sd(dz) == 0 || length(unique(y)) < 2) {
    return(NULL)
  }
  fit <- tryCatch(
    suppressWarnings(if (binomial) stats::glm(y ~ dz, family = stats::binomial()) else stats::lm(y ~ dz)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients
  if (!"dz" %in% rownames(co) || !is.finite(co["dz", 2]) || co["dz", 2] > 1e3) return(NULL)
  unname(co["dz", 1:2])
}

#' @export
print.amr_distance_eval <- function(x, digits = 3, ...) {
  s <- x$settings
  nr <- x$nearest
  cat("<amr_distance_eval> ", sum(!is.na(nr$distance)), " test isolates | ",
      nlevels(nr$band), " distance bands | type: ", s$prediction_type,
      " | threshold: ", if (is.null(s$threshold)) "none" else s$threshold, "\n", sep = "")
  cat("Isolates per band:\n")
  print(table(model = nr$model_id, band = nr$band))
  key_metrics <- c("sensitivity", "specificity", "very_major_error", "major_error",
                   "roc_auc", "brier", "calibration_slope", "essential_agreement")
  bb <- x$by_band[x$by_band$metric %in% key_metrics, , drop = FALSE]
  if (nrow(bb) > 0) {
    cat("\nKey metrics by band:\n")
    fmt <- function(v) ifelse(is.na(v), "", formatC(v, digits = digits, format = "f"))
    shown <- data.frame(model = bb$model_id, drug = bb$drug_name, band = bb$band, metric = bb$metric,
                        estimate = fmt(bb$estimate),
                        interval = ifelse(is.na(bb$lower), "", paste0("[", fmt(bb$lower), ", ", fmt(bb$upper), "]")),
                        stringsAsFactors = FALSE)
    print(shown, row.names = FALSE, right = FALSE)
  }
  tr <- x$trend
  cat("\nTrend with distance:\n")
  for (i in seq_len(nrow(tr))) {
    if (is.na(tr$effect[i])) {
      cat(" ", tr$model_id[i], " | ", tr$drug_name[i], ": ", tr$note[i], "\n", sep = "")
    } else {
      cat(" ", tr$model_id[i], " | ", tr$drug_name[i], ": ", tr$model[i], " = ",
          formatC(tr$effect[i], digits = digits, format = "f"), " [",
          formatC(tr$lower[i], digits = digits, format = "f"), ", ",
          formatC(tr$upper[i], digits = digits, format = "f"), "], n = ", tr$n[i], "\n", sep = "")
    }
  }
  invisible(x)
}

# Calibration intercept and slope per stratum, in the long metric format of
# amr_evaluate(), for one distance band.
calibration_band_rows <- function(sub, intermediate, conf) {
  cal <- amr_calibration(sub, intermediate = intermediate, conf = conf, split = NULL, min_events = 0)
  s <- cal$summary
  if (nrow(s) == 0) return(NULL)
  cols <- stratum_columns(sub)
  rows <- lapply(seq_len(nrow(s)), function(i) {
    first <- s[i, cols, drop = FALSE]
    data.frame(first[rep(1, 2), , drop = FALSE],
               species_name = suppressMessages(AMR::mo_name(s$species[i], language = "en")),
               drug_name = s$drug_name[i], family = "calibration",
               metric = c("calibration_intercept", "calibration_slope"),
               estimate = c(s$intercept[i], s$slope[i]),
               lower = c(s$intercept_lower[i], s$slope_lower[i]),
               upper = c(s$intercept_upper[i], s$slope_upper[i]),
               numerator = NA_real_, denominator = s$n[i],
               row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
