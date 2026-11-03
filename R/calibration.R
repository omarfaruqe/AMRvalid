#' Calibration of predicted probabilities
#'
#' Assesses whether predicted probabilities of resistance match observed
#' frequencies. Discrimination (AUC) and calibration are different
#' properties: a model can rank isolates well while its probabilities are
#' systematically too high or too low, which matters whenever the
#' probability itself is used, for example to set a threshold or to report
#' a risk.
#'
#' @inheritParams amr_evaluate
#' @param bins Number of bins for the binned reliability curve and the
#'   expected calibration error.
#' @param binning `"equal_width"` (default; bins of width 1 / `bins` on the
#'   probability scale) or `"equal_count"` (quantile bins with similar
#'   numbers of isolates).
#' @param by Optional column name to evaluate calibration separately for
#'   each value, in addition to the stratum (e.g. `"lineage"` or a distance
#'   band column you have added).
#' @param min_events Strata with fewer resistant or fewer susceptible
#'   isolates than this are flagged as too sparse for stable estimates.
#'
#' @details
#' Reported per stratum (model, dataset, species, antimicrobial):
#'
#' * **mean predicted** and **observed** proportion resistant, and their
#'   ratio **O/E** (observed / expected; 1 is ideal);
#' * **calibration intercept** (calibration-in-the-large): the intercept of
#'   a logistic model of the outcome with the logit of the predicted
#'   probability as an offset; 0 is ideal, positive means predictions are
#'   too low on average;
#' * **calibration slope**: the coefficient of a logistic regression of the
#'   outcome on the logit of the predicted probability; 1 is ideal, below 1
#'   means predictions are too extreme (over-confident);
#' * **expected calibration error (ECE)**: the count-weighted mean absolute
#'   difference between mean predicted and observed proportion across bins,
#'   and the **maximum calibration error (MCE)**. Both depend on the binning
#'   and are secondary to the slope and intercept.
#'
#' Probabilities are clipped to \[1e-6, 1 - 1e-6\] before taking logits.
#' Intervals for the intercept and slope are Wald intervals assuming
#' independent isolates.
#'
#' The `curve` element holds a loess-smoothed reliability curve (drawn when
#' there are at least 20 isolates and 10 distinct predicted values) and the
#' `bins` element the binned points with Wilson intervals; [plot()] draws
#' both.
#'
#' @return An object of class `amr_calibration`: a list with `summary` (one
#'   row per stratum), `bins` and `curve`.
#' @references
#' Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW (2019).
#' Calibration: the Achilles heel of predictive analytics. *BMC Medicine*
#' 17:230.
#' @export
#' @examples
#' sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10,
#'                     external_lineages = 20, return_distances = FALSE, seed = 1)
#' x <- amr_import_predictions(sim$predictions)
#' cal <- amr_calibration(x[x$dataset_id == "external", ])
#' cal
#' plot(cal)
amr_calibration <- function(x, bins = 10, binning = c("equal_width", "equal_count"), by = NULL,
                            intermediate = c("exclude", "susceptible", "resistant"),
                            conf = 0.95, split = "test", min_events = 10) {
  check_amr_data(x)
  binning <- match.arg(binning)
  intermediate <- match.arg(intermediate)
  if (!identical(x$prediction_type[1], "score")) {
    stop("Calibration needs predicted probabilities (prediction_type 'score').", call. = FALSE)
  }
  if (!is.numeric(bins) || length(bins) != 1 || bins < 2 || bins != round(bins)) {
    stop("`bins` must be a whole number of at least 2.", call. = FALSE)
  }
  if (!is.null(by) && !by %in% names(x)) stop("Column '", by, "' not found.", call. = FALSE)
  x <- select_split(x, split)

  cols <- c(stratum_columns(x), by)
  key <- interaction(as.data.frame(x)[, cols, drop = FALSE], drop = TRUE, lex.order = TRUE, sep = "\r")
  summ <- list(); bin_rows <- list(); curve_rows <- list()
  for (idx in split(seq_len(nrow(x)), key)) {
    sub <- x[idx, , drop = FALSE]
    truth <- collapse_sir(sub$truth_sir, intermediate)
    ok <- !is.na(truth) & !is.na(sub$probability)
    y <- as.integer(truth[ok] == "R"); p <- sub$probability[ok]
    first <- as.data.frame(sub)[1, cols, drop = FALSE]
    first[] <- lapply(first, as.character)
    cs <- calibration_stats(y, p, bins, binning, conf)
    sparse <- min(sum(y == 1), sum(y == 0)) < min_events
    summ[[length(summ) + 1]] <- data.frame(first, cs$summary, sparse = sparse,
                                           row.names = NULL, stringsAsFactors = FALSE)
    if (nrow(cs$bins) > 0) bin_rows[[length(bin_rows) + 1]] <- data.frame(first[rep(1, nrow(cs$bins)), , drop = FALSE],
                                                                          cs$bins, row.names = NULL)
    if (nrow(cs$curve) > 0) curve_rows[[length(curve_rows) + 1]] <- data.frame(first[rep(1, nrow(cs$curve)), , drop = FALSE],
                                                                              cs$curve, row.names = NULL)
  }
  bind <- function(l) if (length(l)) { r <- do.call(rbind, l); rownames(r) <- NULL; r } else data.frame()
  summary <- bind(summ)
  summary$drug_name <- suppressMessages(AMR::ab_name(summary$drug, language = "en"))
  structure(list(summary = summary, bins = bind(bin_rows), curve = bind(curve_rows),
                 settings = list(bins = bins, binning = binning, by = by, intermediate = intermediate,
                                 conf = conf, split = split, min_events = min_events, strata = cols)),
            class = "amr_calibration")
}

# Calibration statistics for one set of outcomes and probabilities.
calibration_stats <- function(y, p, bins, binning, conf) {
  n <- length(y)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  out <- data.frame(n = n, n_resistant = sum(y), mean_predicted = NA_real_, observed = NA_real_,
                    oe_ratio = NA_real_, intercept = NA_real_, intercept_lower = NA_real_,
                    intercept_upper = NA_real_, slope = NA_real_, slope_lower = NA_real_,
                    slope_upper = NA_real_, ece = NA_real_, mce = NA_real_, note = "",
                    stringsAsFactors = FALSE)
  empty <- list(summary = out, bins = data.frame(), curve = data.frame())
  if (n == 0) { out$note <- "No evaluable isolates."; empty$summary <- out; return(empty) }
  out$mean_predicted <- mean(p); out$observed <- mean(y)
  out$oe_ratio <- if (mean(p) > 0) mean(y) / mean(p) else NA_real_

  # bins
  if (binning == "equal_width") {
    br <- seq(0, 1, length.out = bins + 1)
  } else {
    br <- unique(stats::quantile(p, probs = seq(0, 1, length.out = bins + 1), names = FALSE))
    if (length(br) < 2) br <- c(0, 1)
    br[1] <- 0; br[length(br)] <- 1
  }
  b <- cut(p, br, include.lowest = TRUE)
  bt <- do.call(rbind, lapply(split(seq_len(n), b, drop = TRUE), function(i) {
    ci <- amr_prop_ci(sum(y[i]), length(i), conf = conf)
    data.frame(bin = as.character(b[i[1]]), n = length(i), mean_predicted = mean(p[i]),
               observed = mean(y[i]), lower = ci$lower, upper = ci$upper, stringsAsFactors = FALSE)
  }))
  rownames(bt) <- NULL
  gap <- abs(bt$mean_predicted - bt$observed)
  out$ece <- sum(bt$n * gap) / n
  out$mce <- max(gap)

  if (length(unique(y)) < 2) {
    out$note <- "Only one phenotype class: intercept and slope undefined."
    return(list(summary = out, bins = bt, curve = data.frame()))
  }
  lp <- stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  fit_int <- tryCatch(suppressWarnings(stats::glm(y ~ 1 + offset(lp), family = stats::binomial())),
                      error = function(e) NULL)
  if (!is.null(fit_int)) {
    co <- summary(fit_int)$coefficients
    out$intercept <- co[1, 1]
    out$intercept_lower <- co[1, 1] - z * co[1, 2]; out$intercept_upper <- co[1, 1] + z * co[1, 2]
  }
  if (stats::sd(lp) > 0) {
    fit_slope <- tryCatch(suppressWarnings(stats::glm(y ~ lp, family = stats::binomial())),
                          error = function(e) NULL)
    if (!is.null(fit_slope) && fit_slope$converged) {
      co <- summary(fit_slope)$coefficients
      if ("lp" %in% rownames(co) && co["lp", 2] < 1e3) {
        out$slope <- co["lp", 1]
        out$slope_lower <- co["lp", 1] - z * co["lp", 2]; out$slope_upper <- co["lp", 1] + z * co["lp", 2]
      }
    }
  } else {
    out$note <- "All predicted probabilities are equal: slope undefined."
  }

  curve <- data.frame()
  if (n >= 20 && length(unique(p)) >= 10) {
    fit <- tryCatch(suppressWarnings(stats::loess(y ~ p, span = 0.75, degree = 2)), error = function(e) NULL)
    if (!is.null(fit)) {
      grid <- seq(min(p), max(p), length.out = 50)
      curve <- data.frame(predicted = grid,
                          observed = pmin(1, pmax(0, stats::predict(fit, data.frame(p = grid)))))
    }
  }
  list(summary = out, bins = bt, curve = curve)
}

#' @export
print.amr_calibration <- function(x, digits = 3, ...) {
  s <- x$summary
  cat("<amr_calibration> ", nrow(s), " stratum/strata | ", x$settings$bins, " ",
      sub("_", "-", x$settings$binning), " bins | ", x$settings$conf * 100, "% Wald intervals\n", sep = "")
  fmt <- function(v) ifelse(is.na(v), "", formatC(v, digits = digits, format = "f"))
  ci <- function(e, l, u) ifelse(is.na(e), "", paste0(fmt(e), " [", fmt(l), ", ", fmt(u), "]"))
  shown <- data.frame(model = s$model_id, drug = s$drug_name, n = s$n,
                      `O/E` = fmt(s$oe_ratio), intercept = ci(s$intercept, s$intercept_lower, s$intercept_upper),
                      slope = ci(s$slope, s$slope_lower, s$slope_upper), ECE = fmt(s$ece),
                      check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(x$settings$by)) shown <- cbind(shown[, 1:2], by = s[[x$settings$by]], shown[, -(1:2)])
  print(shown, row.names = FALSE, right = FALSE)
  cat("Ideal: O/E = 1, intercept = 0, slope = 1. Slope < 1: predictions too extreme.\n")
  if (any(s$sparse)) cat("! Some strata have fewer than ", x$settings$min_events,
                         " resistant or susceptible isolates; estimates are unstable.\n", sep = "")
  notes <- unique(s$note[nzchar(s$note)])
  for (n in notes) cat("i ", n, "\n", sep = "")
  invisible(x)
}

#' Plot a reliability diagram
#'
#' Draws one panel per stratum (model, dataset, species and antimicrobial,
#' plus the `by` column if used) with the binned observed proportions, their
#' Wilson intervals and the loess-smoothed calibration curve.
#'
#' @param x An `amr_calibration` object.
#' @param max_panels Maximum number of panels to draw.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
#' @examples
#' sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10,
#'                     external_lineages = 20, return_distances = FALSE, seed = 1)
#' x <- amr_import_predictions(sim$predictions)
#' plot(amr_calibration(x[x$dataset_id == "external" & x$model_id == "marker/lineage", ]))
plot.amr_calibration <- function(x, max_panels = 12, ...) {
  cols <- x$settings$strata
  s <- x$summary
  keys <- do.call(paste, c(s[, cols, drop = FALSE], sep = "\r"))
  if (length(keys) > max_panels) {
    message("Drawing the first ", max_panels, " of ", length(keys), " panels; set `max_panels` to show more.")
    keys <- keys[seq_len(max_panels)]
  }
  key_of <- function(d) if (nrow(d)) do.call(paste, c(d[, cols, drop = FALSE], sep = "\r")) else character()
  bkey <- key_of(x$bins); ckey <- key_of(x$curve)
  n <- length(keys)
  nc <- min(n, 3); nr <- ceiling(n / nc)
  old <- graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1), las = 1, col.axis = "#6b6b6b",
                       fg = "#6b6b6b", col.lab = "#2b2b2b", mgp = c(2.4, 0.6, 0), tcl = -0.25)
  on.exit(graphics::par(old), add = TRUE)
  col_curve <- "#2a78d6"; col_pts <- "#eb6834"
  for (k in keys) {
    i <- match(k, keys)
    row <- s[match(k, do.call(paste, c(s[, cols, drop = FALSE], sep = "\r"))), ]
    graphics::plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Predicted probability of resistance",
                   ylab = "Observed proportion resistant", bty = "n")
    graphics::abline(h = seq(0, 1, 0.2), v = seq(0, 1, 0.2), col = "#e6e6e3")
    graphics::abline(0, 1, lty = 2, col = "#6b6b6b")
    b <- x$bins[bkey == k, , drop = FALSE]
    if (nrow(b)) {
      graphics::arrows(b$mean_predicted, b$lower, b$mean_predicted, b$upper, angle = 90, code = 3,
                       length = 0.03, col = col_pts)
      graphics::points(b$mean_predicted, b$observed, pch = 21, bg = col_pts, col = "white",
                       cex = 0.8 + 1.6 * sqrt(b$n / max(b$n)))
    }
    cv <- x$curve[ckey == k, , drop = FALSE]
    if (nrow(cv)) graphics::lines(cv$predicted, cv$observed, col = col_curve, lwd = 2)
    title_txt <- paste(c(row$model_id, if (length(unique(s$dataset_id)) > 1) row$dataset_id,
                         row$drug_name, if (!is.null(x$settings$by)) row[[x$settings$by]]), collapse = " | ")
    graphics::title(main = title_txt, cex.main = 0.9, font.main = 1, adj = 0, col.main = "#2b2b2b")
    sub_txt <- if (is.na(row$slope)) "slope n/a" else sprintf("slope %.2f, intercept %.2f", round(row$slope, 2) + 0, round(row$intercept, 2) + 0)
    graphics::mtext(sub_txt, side = 3, line = 0.1, adj = 0, cex = 0.7, col = "#6b6b6b")
    if (i == 1) {
      graphics::legend("topleft", bty = "n", cex = 0.75, text.col = "#2b2b2b",
                       legend = c("binned observed (95% CI)", "smoothed curve", "perfect calibration"),
                       pch = c(21, NA, NA), pt.bg = c(col_pts, NA, NA), pt.cex = 1.4,
                       col = c("white", col_curve, "#6b6b6b"),
                       lty = c(NA, 1, 2), lwd = c(NA, 2, 1))
    }
  }
  invisible(x)
}
