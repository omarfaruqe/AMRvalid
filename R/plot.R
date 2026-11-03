# Plot methods. Base graphics only, so plotting adds no dependencies.
# Colours: a colourblind-safe categorical palette, assigned to models in a
# fixed order (sorted model identifiers), never cycled.

amr_palette <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
ink <- "#2b2b2b"
ink_muted <- "#6b6b6b"
grid_col <- "#e6e6e3"

model_colours <- function(models) {
  m <- sort(unique(models))
  if (length(m) > length(amr_palette)) {
    stop("Plots show at most ", length(amr_palette), " models; subset the results first.", call. = FALSE)
  }
  stats::setNames(amr_palette[seq_along(m)], m)
}

with_plot_par <- function(code, mfrow = c(1, 1), mar = c(4, 12, 3, 1), widths = NULL, oma = c(0, 0, 0, 0)) {
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = mar, las = 1, col.axis = ink_muted, fg = ink_muted, oma = oma,
                col.lab = ink, col.main = ink, mgp = c(2.4, 0.6, 0), tcl = -0.25)
  if (!is.null(widths)) {
    graphics::layout(matrix(seq_along(widths), nrow = 1), widths = widths)
    graphics::par(cex = 0.9)
  } else {
    graphics::par(mfrow = mfrow)
  }
  force(code)
}

# Panel widths so that every panel gets the same plotting width when the
# first one also carries the row labels. Widths are in text lines, computed
# from the current device width (0.9 text size, as set in with_plot_par()).
panel_widths <- function(n, label_margin = 13) {
  total <- graphics::par("din")[1] / (0.2 * 0.9)
  first_margin <- label_margin + 1
  other_margin <- 2
  plot_w <- max(1, (total - first_margin - other_margin * (n - 1)) / n)
  c(plot_w + first_margin, rep(plot_w + other_margin, n - 1))
}

# One dot-and-interval panel. Rows run top to bottom. The first panel of a
# row of panels carries the row labels and a wide left margin.
forest_panel <- function(labels, est, lower, upper, cols, main, xlim = NULL, ref = NULL, show_labels = TRUE,
                         label_margin = 13) {
  graphics::par(mar = c(3, if (show_labels) label_margin else 1, 3, 1))
  n <- length(est)
  y <- rev(seq_len(n))
  if (is.null(xlim)) {
    rng <- range(c(est, lower, upper, ref), na.rm = TRUE)
    if (!all(is.finite(rng))) rng <- c(0, 1)
    pad <- max(diff(rng) * 0.08, 0.02)
    xlim <- c(rng[1] - pad, rng[2] + pad)
  }
  graphics::plot(NA, xlim = xlim, ylim = c(0.5, n + 0.5), yaxt = "n", xlab = "", ylab = "", bty = "n")
  graphics::abline(v = pretty(xlim), col = grid_col)
  if (!is.null(ref)) graphics::abline(v = ref, lty = 2, col = ink_muted)
  ok <- !is.na(lower) & !is.na(upper)
  graphics::segments(lower[ok], y[ok], upper[ok], y[ok], col = cols[ok], lwd = 2)
  graphics::points(est, y, pch = 21, bg = cols, col = "white", cex = 1.6, lwd = 1.5)
  if (show_labels) graphics::axis(2, at = y, labels = labels, lwd = 0, cex.axis = 0.8, col.axis = ink)
  graphics::title(main = main, adj = 0, cex.main = 0.95, font.main = 1)
}

model_legend <- function(cols) {
  if (length(cols) < 2) return(invisible(NULL))
  graphics::legend("bottomright", legend = names(cols), pt.bg = cols, pch = 21, col = "white", pt.cex = 1.4,
                   bty = "n", cex = 0.75, text.col = ink, inset = c(0, 0.01))
}

row_labels <- function(d, extra = NULL) {
  parts <- list(d$drug_name)
  if ("dataset_id" %in% names(d) && length(unique(d$dataset_id)) > 1) parts <- c(parts, list(d$dataset_id))
  for (e in extra) parts <- c(parts, list(as.character(d[[e]])))
  lab <- do.call(paste, c(parts, sep = " | "))
  if ("model_id" %in% names(d) && length(unique(d$model_id)) > 1) lab <- paste0(lab, "\n", d$model_id)
  lab
}

#' Plot evaluation results
#'
#' Dot-and-interval plots of selected metrics, one panel per metric and one
#' row per stratum. Rows are labelled with the antimicrobial (and model,
#' dataset and `by` values where there are several); colour repeats the
#' model.
#'
#' @param x An `amr_eval` object from [amr_evaluate()].
#' @param metrics Metrics to draw. By default sensitivity and specificity
#'   when available, otherwise ROC-AUC (threshold-free scores) or essential
#'   agreement (MIC predictions).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
#' @examples
#' x <- amr_import_predictions(amr_example_data())
#' plot(amr_evaluate(x), metrics = c("very_major_error", "major_error"))
plot.amr_eval <- function(x, metrics = NULL, ...) {
  d <- as.data.frame(x)
  if (is.null(metrics)) {
    metrics <- if (all(c("sensitivity", "specificity") %in% d$metric)) c("sensitivity", "specificity")
               else if ("roc_auc" %in% d$metric) "roc_auc" else if ("essential_agreement" %in% d$metric)
                 "essential_agreement" else stop("No plottable metrics found.", call. = FALSE)
  }
  bad <- setdiff(metrics, d$metric)
  if (length(bad)) stop("Metrics not in the results: ", paste(bad, collapse = ", "), call. = FALSE)
  by <- attr(x, "settings")$by
  cols <- model_colours(d$model_id)
  first <- d[d$metric == metrics[1], , drop = FALSE]
  labels <- row_labels(first, by)
  with_plot_par({
    for (i in seq_along(metrics)) {
      m <- d[d$metric == metrics[i], , drop = FALSE]
      m <- m[match(paste(first$model_id, first$dataset_id, first$drug, do.call(paste, c(first[by], ""))),
                   paste(m$model_id, m$dataset_id, m$drug, do.call(paste, c(m[by], "")))), ]
      forest_panel(labels, m$estimate, m$lower, m$upper, cols[m$model_id], main = metrics[i],
                   xlim = if (metrics[i] %in% c("roc_auc", "brier", "log_loss", "mcc")) NULL else c(0, 1),
                   ref = if (metrics[i] == "roc_auc") 0.5 else NULL, show_labels = i == 1)
    }
  }, widths = panel_widths(length(metrics)))
  invisible(x)
}

#' Plot performance by distance to the training set
#'
#' Draws one metric across the distance bands, with intervals: one line per
#' stratum, coloured by model, with line type distinguishing datasets.
#'
#' @param x An `amr_distance_eval` object from [amr_distance_eval()].
#' @param metric Metric to draw. Defaults to ROC-AUC for scores without a
#'   threshold and accuracy otherwise.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
#' @examples
#' sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10,
#'                     external_lineages = 20, seed = 1)
#' x <- amr_import_predictions(sim$predictions)
#' d <- amr_distance_eval(x[x$model_id == "knn/random", ], sim$distances,
#'                        train_manifest = sim$train_manifest)
#' plot(d)
plot.amr_distance_eval <- function(x, metric = NULL, ...) {
  bb <- x$by_band
  if (is.null(metric)) metric <- if ("accuracy" %in% bb$metric) "accuracy" else "roc_auc"
  bb <- bb[bb$metric == metric, , drop = FALSE]
  if (nrow(bb) == 0) stop("Metric '", metric, "' is not in the results.", call. = FALSE)
  bands <- levels(bb$band)
  cols <- model_colours(bb$model_id)
  datasets <- sort(unique(bb$dataset_id))
  ltys <- stats::setNames(rep_len(c(1, 2, 3, 4), length(datasets)), datasets)
  key <- paste(bb$model_id, bb$dataset_id, bb$drug, sep = " | ")
  keys <- unique(key)
  offset <- if (length(keys) > 1) seq(-0.12, 0.12, length.out = length(keys)) else 0
  with_plot_par({
    graphics::plot(NA, xlim = c(0.6, length(bands) + 0.4), ylim = c(0, 1), xaxt = "n", bty = "n",
                   xlab = "Distance to nearest training isolate", ylab = metric)
    graphics::abline(h = seq(0, 1, 0.2), col = grid_col)
    if (metric == "roc_auc") graphics::abline(h = 0.5, lty = 2, col = ink_muted)
    graphics::axis(1, at = seq_along(bands), labels = bands, lwd = 0, cex.axis = 0.8)
    leg <- list(label = character(), col = character(), lty = numeric())
    for (k in seq_along(keys)) {
      s <- bb[key == keys[k], , drop = FALSE]
      xs <- match(as.character(s$band), bands) + offset[k]
      col <- cols[[s$model_id[1]]]
      lty <- ltys[[s$dataset_id[1]]]
      graphics::lines(xs, s$estimate, col = col, lwd = 2, lty = lty)
      ok <- !is.na(s$lower)
      graphics::segments(xs[ok], s$lower[ok], xs[ok], s$upper[ok], col = col, lwd = 2)
      graphics::points(xs, s$estimate, pch = 21, bg = col, col = "white", cex = 1.6, lwd = 1.5)
      parts <- c(s$model_id[1], if (length(datasets) > 1) s$dataset_id[1],
                 if (length(unique(bb$drug)) > 1) s$drug_name[1])
      leg$label <- c(leg$label, paste(parts, collapse = " | "))
      leg$col <- c(leg$col, col); leg$lty <- c(leg$lty, lty)
    }
    if (length(keys) > 1) {
      graphics::legend("bottomleft", legend = leg$label, col = leg$col, lty = leg$lty, lwd = 2, pch = 21,
                       pt.bg = leg$col, pt.cex = 1.2, bty = "n", cex = 0.75, text.col = ink, seg.len = 3)
    }
    graphics::title(main = paste(metric, "by distance band"), adj = 0, cex.main = 0.95, font.main = 1)
  }, mar = c(4, 4, 3, 1))
  invisible(x)
}

#' Plot bootstrap intervals
#'
#' Dot-and-interval plots of bootstrap estimates, one panel per metric.
#'
#' @param x An `amr_bootstrap` object from [amr_bootstrap()].
#' @param metrics Metrics to draw; by default sensitivity and specificity
#'   when available, otherwise ROC-AUC.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
#' @examples
#' x <- amr_import_predictions(amr_example_data())
#' b <- amr_bootstrap(x, resample_by = "genomic_cluster", B = 100, seed = 1)
#' plot(b)
plot.amr_bootstrap <- function(x, metrics = NULL, ...) {
  d <- as.data.frame(x)
  if (is.null(metrics)) {
    metrics <- if (all(c("sensitivity", "specificity") %in% d$metric)) c("sensitivity", "specificity") else "roc_auc"
  }
  bad <- setdiff(metrics, d$metric)
  if (length(bad)) stop("Metrics not in the results: ", paste(bad, collapse = ", "), call. = FALSE)
  cols <- model_colours(d$model_id)
  first <- d[d$metric == metrics[1], , drop = FALSE]
  labels <- row_labels(first)
  with_plot_par({
    for (i in seq_along(metrics)) {
      m <- d[d$metric == metrics[i], , drop = FALSE]
      m <- m[match(paste(first$model_id, first$dataset_id, first$drug), paste(m$model_id, m$dataset_id, m$drug)), ]
      forest_panel(labels, m$estimate, m$lower, m$upper, cols[m$model_id], main = metrics[i],
                   xlim = if (metrics[i] %in% c("brier", "log_loss", "mcc")) NULL else c(0, 1),
                   ref = if (metrics[i] == "roc_auc") 0.5 else NULL, show_labels = i == 1)
    }
    graphics::mtext(paste0(100 * attr(x, "settings")$conf, "% intervals from a bootstrap resampling ",
                           attr(x, "settings")$resample_by, " (B = ", attr(x, "settings")$B, ")"),
                    side = 3, outer = TRUE, adj = 0.01, line = 0.3, cex = 0.8, col = ink_muted)
  }, widths = panel_widths(length(metrics)), oma = c(0, 0, 1.5, 0))
  invisible(x)
}

#' Plot a paired model comparison
#'
#' Differences between the two models (second minus first) with bootstrap
#' intervals and a reference line at zero.
#'
#' @param x An `amr_paired_compare` object from [amr_paired_compare()].
#' @param metrics Metrics to draw; by default every metric in the results.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
#' @examples
#' x <- amr_import_predictions(amr_example_data())
#' mero <- x[x$drug == "MEM", ]
#' p <- amr_paired_compare(mero, models = c("carbapenemase genes",
#'                         "carbapenemase or porin + ESBL/AmpC"), B = 100, seed = 1)
#' plot(p)
plot.amr_paired_compare <- function(x, metrics = NULL, ...) {
  d <- as.data.frame(x)
  if (is.null(metrics)) metrics <- unique(d$metric)
  d <- d[d$metric %in% metrics, , drop = FALSE]
  if (nrow(d) == 0) stop("None of the requested metrics are in the results.", call. = FALSE)
  s <- attr(x, "settings")
  lab <- if (length(unique(d$drug_name)) > 1) paste(d$drug_name, d$metric, sep = " | ") else d$metric
  col <- ifelse(is.na(d$lower), ink_muted, ifelse(d$lower > 0 | d$upper < 0, amr_palette[1], ink_muted))
  with_plot_par({
    forest_panel(lab, d$difference, d$lower, d$upper, col, ref = 0,
                 main = paste0(s$models[2], " minus ", s$models[1]))
    graphics::mtext("Blue: interval excludes zero", side = 1, line = 2, adj = 1, cex = 0.7, col = ink_muted)
  }, mar = c(4, 12, 3, 1))
  invisible(x)
}
