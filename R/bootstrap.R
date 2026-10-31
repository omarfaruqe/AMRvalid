#' Group-aware bootstrap intervals
#'
#' Resamples whole groups (isolates, lineages or genomic clusters) with
#' replacement and recomputes metrics in every replicate. Resampling rows
#' as if isolates were independent understates uncertainty when related
#' isolates share outcomes; resampling clusters respects that dependence.
#'
#' @inheritParams amr_species_eval
#' @param group Column defining the resampling unit, e.g. `"isolate_id"`
#'   (default), `"lineage"` or `"genomic_cluster"`. Rows with a missing
#'   group label are treated as their own single-isolate groups.
#' @param B Number of bootstrap replicates.
#' @param seed Optional random seed for reproducibility. The global random
#'   state is restored afterwards.
#' @param min_groups Strata with fewer independent groups than this are
#'   flagged: cluster bootstrap intervals are unreliable with few clusters.
#'
#' @details
#' Metrics are computed on S/R categories after applying the `intermediate`
#' policy, so the very major error rate equals 1 - sensitivity and the
#' major error rate equals 1 - specificity. Minor errors are not
#' bootstrapped. Intervals are percentile intervals; `n_valid` counts the
#' replicates in which a metric was defined (for example, sensitivity is
#' undefined in a replicate that drew no resistant isolates). Treat
#' intervals with many invalid replicates with caution.
#'
#' When every replicate gives the same value (typically a sensitivity of 1
#' or a very major error rate of 0), the percentile interval has zero width
#' and cannot express uncertainty. Such rows are marked `zero_width`; use
#' the Wilson or Clopper-Pearson intervals from [amr_species_eval()] for
#' them instead.
#'
#' @return An object of class `amr_bootstrap`: a data frame with the
#'   stratum columns, `metric`, `estimate`, `lower`, `upper`, `n_groups`,
#'   `B`, `n_valid`, `few_groups` and `zero_width`.
#' @export
#' @examples
#' set.seed(2)
#' preds <- data.frame(
#'   isolate_id = sprintf("iso%02d", 1:30),
#'   species = "K. pneumoniae",
#'   drug = "meropenem",
#'   truth = rep(c("R", "S", "S"), 10),
#'   probability = runif(30),
#'   lineage = rep(sprintf("ST%d", 1:10), each = 3)
#' )
#' amr_bootstrap(amr_import_predictions(preds), group = "lineage", B = 200,
#'               threshold = 0.5, seed = 1)
amr_bootstrap <- function(x, group = "isolate_id", B = 2000, threshold = NULL,
                          intermediate = c("exclude", "susceptible", "resistant"),
                          conf = 0.95, seed = NULL, split = "test", min_groups = 20) {
  check_amr_data(x)
  intermediate <- match.arg(intermediate)
  check_B(B)
  if (!group %in% names(x)) stop("Group column '", group, "' not found.", call. = FALSE)
  ptype <- x$prediction_type[1]
  if (identical(ptype, "mic")) stop("Bootstrap of MIC metrics is not implemented yet.", call. = FALSE)
  check_threshold(threshold, ptype)
  x <- select_split(x, split)
  restore <- local_seed(seed)
  on.exit(restore(), add = TRUE)

  cols <- stratum_columns(x)
  key <- interaction(as.data.frame(x)[, cols], drop = TRUE, lex.order = TRUE, sep = "\r")
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  pieces <- lapply(split(seq_len(nrow(x)), key), function(idx) {
    sub <- x[idx, , drop = FALSE]
    truth <- collapse_sir(sub$truth_sir, intermediate)
    pred <- predicted_class(sub, threshold, intermediate)
    prob <- if (identical(ptype, "score")) sub$probability else NULL
    ok <- !is.na(truth) & (if (is.null(pred)) TRUE else !is.na(pred)) &
      (if (is.null(prob)) TRUE else !is.na(prob))
    y <- as.integer(truth[ok] == "R")
    yhat <- if (is.null(pred)) NULL else as.integer(pred[ok] == "R")
    pp <- if (is.null(prob)) NULL else prob[ok]
    g <- group_labels(sub, group)[ok]

    est <- fast_metrics(y, yhat, pp)
    groups <- split(seq_along(y), g)
    reps <- matrix(NA_real_, nrow = B, ncol = length(est), dimnames = list(NULL, names(est)))
    for (b in seq_len(B)) {
      draw <- unlist(groups[sample.int(length(groups), length(groups), replace = TRUE)], use.names = FALSE)
      reps[b, ] <- fast_metrics(y[draw], if (is.null(yhat)) NULL else yhat[draw],
                                if (is.null(pp)) NULL else pp[draw])
    }
    qs <- apply(reps, 2, function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 2) c(NA_real_, NA_real_) else stats::quantile(v, probs, names = FALSE, type = 7)
    })
    first <- as.data.frame(sub)[1, cols, drop = FALSE]
    first[] <- lapply(first, as.character)
    data.frame(first[rep(1, length(est)), , drop = FALSE], metric = names(est), estimate = unname(est),
               lower = qs[1, ], upper = qs[2, ], n_groups = length(groups), B = B,
               n_valid = colSums(is.finite(reps)), few_groups = length(groups) < min_groups,
               zero_width = !is.na(qs[1, ]) & abs(qs[2, ] - qs[1, ]) < 1e-10,
               row.names = NULL, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, pieces)
  out$drug_name <- suppressMessages(AMR::ab_name(out$drug, language = "en"))
  rownames(out) <- NULL
  structure(out, class = c("amr_bootstrap", "data.frame"),
            settings = list(group = group, B = B, threshold = threshold, intermediate = intermediate,
                            conf = conf, seed = seed, split = split, min_groups = min_groups))
}

#' @export
print.amr_bootstrap <- function(x, digits = 3, ...) {
  s <- attr(x, "settings")
  cat("<amr_bootstrap> resampling unit: ", s$group, " | B = ", s$B, " | ", s$conf * 100,
      "% percentile intervals | threshold: ", if (is.null(s$threshold)) "none" else s$threshold, "\n", sep = "")
  d <- as.data.frame(x)
  fmt <- function(v) ifelse(is.na(v), "", formatC(v, digits = digits, format = "f"))
  shown <- data.frame(model = d$model_id, drug = d$drug_name, metric = d$metric, estimate = fmt(d$estimate),
                      interval = ifelse(is.na(d$lower), "", paste0("[", fmt(d$lower), ", ", fmt(d$upper), "]")),
                      groups = d$n_groups, valid = paste0(d$n_valid, "/", d$B), stringsAsFactors = FALSE)
  print(shown, row.names = FALSE, right = FALSE)
  if (any(d$few_groups)) {
    cat("! Fewer than ", s$min_groups, " independent groups in some strata; cluster bootstrap intervals",
        " may be too narrow.\n", sep = "")
  }
  if (any(d$zero_width)) {
    cat("! ", sum(d$zero_width), " interval(s) have zero width (every replicate gave the same value);",
        " use the Wilson intervals from amr_species_eval() for these.\n", sep = "")
  }
  invisible(x)
}

#' Paired comparison of two models on the same isolates
#'
#' Restricts both models to the isolates they share in each stratum and
#' bootstraps the difference in each metric, resampling the same groups for
#' both models so the pairing is kept in every replicate.
#'
#' @inheritParams amr_bootstrap
#' @param models Character vector of two `model_id` values. Differences are
#'   reported as `models[2] - models[1]`. If `NULL`, the data must contain
#'   exactly two models.
#' @param threshold Probability threshold for scores: a single number, or a
#'   named numeric vector with one threshold per model.
#'
#' @details
#' Besides the bootstrap intervals, an exact McNemar test compares the
#' discordant pairs of correct and incorrect classifications. It assumes
#' independent isolates, so prefer the group-aware intervals when isolates
#' are clustered.
#'
#' @return An object of class `amr_paired_compare`: a data frame with
#'   `dataset_id`, `species`, `drug`, `metric`, the estimate for each model
#'   (`model_a`, `model_b`), `difference`, `lower`, `upper`, `n_matched`,
#'   `n_groups`, `n_valid` and `mcnemar_p` (on the `accuracy` row).
#' @export
#' @examples
#' set.seed(3)
#' ids <- sprintf("iso%02d", 1:40)
#' truth <- rep(c("R", "S"), 20)
#' preds <- rbind(
#'   data.frame(isolate_id = ids, truth = truth, model_id = "old",
#'              probability = ifelse(truth == "R", 0.6, 0.4) + rnorm(40, 0, 0.2)),
#'   data.frame(isolate_id = ids, truth = truth, model_id = "new",
#'              probability = ifelse(truth == "R", 0.8, 0.2) + rnorm(40, 0, 0.1))
#' )
#' preds$probability <- pmin(pmax(preds$probability, 0), 1)
#' preds$species <- "E. coli"
#' preds$drug <- "ciprofloxacin"
#' amr_paired_compare(amr_import_predictions(preds), models = c("old", "new"),
#'                    threshold = 0.5, B = 200, seed = 1)
amr_paired_compare <- function(x, models = NULL, group = "isolate_id", B = 2000, threshold = NULL,
                               intermediate = c("exclude", "susceptible", "resistant"),
                               conf = 0.95, seed = NULL, split = "test", min_groups = 20) {
  check_amr_data(x)
  intermediate <- match.arg(intermediate)
  check_B(B)
  if (!group %in% names(x)) stop("Group column '", group, "' not found.", call. = FALSE)
  ptype <- x$prediction_type[1]
  if (identical(ptype, "mic")) stop("Paired comparison of MIC predictions is not implemented yet.", call. = FALSE)
  if (is.null(models)) models <- unique(x$model_id)
  if (length(models) != 2) stop("Exactly two models are needed; set `models`.", call. = FALSE)
  if (!all(models %in% x$model_id)) {
    stop("Models not found: ", paste(setdiff(models, x$model_id), collapse = ", "), call. = FALSE)
  }
  thr <- function(m) {
    if (is.null(threshold) || length(threshold) == 1) return(threshold)
    if (is.null(names(threshold)) || !m %in% names(threshold)) {
      stop("A threshold vector must be named with both model identifiers.", call. = FALSE)
    }
    threshold[[m]]
  }
  for (m in models) check_threshold(thr(m), ptype)
  x <- select_split(x, split)
  restore <- local_seed(seed)
  on.exit(restore(), add = TRUE)

  cols <- intersect(c("dataset_id", "species", "drug"), names(x))
  key <- interaction(as.data.frame(x)[, cols], drop = TRUE, lex.order = TRUE, sep = "\r")
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  pieces <- lapply(split(seq_len(nrow(x)), key), function(idx) {
    sub <- x[idx, , drop = FALSE]
    prep <- lapply(models, function(m) {
      s <- sub[sub$model_id == m, , drop = FALSE]
      truth <- collapse_sir(s$truth_sir, intermediate)
      pred <- predicted_class(s, thr(m), intermediate)
      prob <- if (identical(ptype, "score")) s$probability else NULL
      ok <- !is.na(truth) & (if (is.null(pred)) TRUE else !is.na(pred)) &
        (if (is.null(prob)) TRUE else !is.na(prob))
      list(id = s$isolate_id[ok], y = as.integer(truth[ok] == "R"),
           yhat = if (is.null(pred)) NULL else as.integer(pred[ok] == "R"),
           p = if (is.null(prob)) NULL else prob[ok], g = group_labels(s, group)[ok])
    })
    shared <- intersect(prep[[1]]$id, prep[[2]]$id)
    first <- as.data.frame(sub)[1, cols, drop = FALSE]
    first[] <- lapply(first, as.character)
    if (length(shared) == 0) return(NULL)
    a <- lapply(prep[[1]], function(v) if (is.null(v)) NULL else v[match(shared, prep[[1]]$id)])
    b <- lapply(prep[[2]], function(v) if (is.null(v)) NULL else v[match(shared, prep[[2]]$id)])
    if (any(a$y != b$y)) stop("Reference phenotypes differ between models for the same isolates.", call. = FALSE)
    est_a <- fast_metrics(a$y, a$yhat, a$p)
    est_b <- fast_metrics(b$y, b$yhat, b$p)
    groups <- split(seq_along(shared), a$g)
    reps <- matrix(NA_real_, nrow = B, ncol = length(est_a), dimnames = list(NULL, names(est_a)))
    for (r in seq_len(B)) {
      draw <- unlist(groups[sample.int(length(groups), length(groups), replace = TRUE)], use.names = FALSE)
      reps[r, ] <- fast_metrics(b$y[draw], if (is.null(b$yhat)) NULL else b$yhat[draw],
                                if (is.null(b$p)) NULL else b$p[draw]) -
        fast_metrics(a$y[draw], if (is.null(a$yhat)) NULL else a$yhat[draw],
                     if (is.null(a$p)) NULL else a$p[draw])
    }
    qs <- apply(reps, 2, function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 2) c(NA_real_, NA_real_) else stats::quantile(v, probs, names = FALSE, type = 7)
    })
    p_mcnemar <- rep(NA_real_, length(est_a))
    if (!is.null(a$yhat)) {
      ca <- a$yhat == a$y; cb <- b$yhat == b$y
      n01 <- sum(ca & !cb); n10 <- sum(!ca & cb)
      if (n01 + n10 > 0) {
        p_mcnemar[names(est_a) == "accuracy"] <- stats::binom.test(n10, n01 + n10, 0.5)$p.value
      }
    }
    data.frame(first[rep(1, length(est_a)), , drop = FALSE], metric = names(est_a),
               model_a = unname(est_a), model_b = unname(est_b), difference = unname(est_b - est_a),
               lower = qs[1, ], upper = qs[2, ], n_matched = length(shared),
               n_unmatched = length(union(prep[[1]]$id, prep[[2]]$id)) - length(shared),
               n_groups = length(groups), n_valid = colSums(is.finite(reps)),
               few_groups = length(groups) < min_groups, mcnemar_p = p_mcnemar,
               row.names = NULL, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, pieces)
  if (is.null(out)) stop("The two models share no evaluable isolates.", call. = FALSE)
  out$drug_name <- suppressMessages(AMR::ab_name(out$drug, language = "en"))
  rownames(out) <- NULL
  structure(out, class = c("amr_paired_compare", "data.frame"),
            settings = list(models = models, group = group, B = B, threshold = threshold,
                            intermediate = intermediate, conf = conf, seed = seed, split = split,
                            min_groups = min_groups))
}

#' @export
print.amr_paired_compare <- function(x, digits = 3, ...) {
  s <- attr(x, "settings")
  cat("<amr_paired_compare> ", s$models[2], " - ", s$models[1], " | resampling unit: ", s$group,
      " | B = ", s$B, " | ", s$conf * 100, "% percentile intervals\n", sep = "")
  d <- as.data.frame(x)
  fmt <- function(v) ifelse(is.na(v), "", formatC(v, digits = digits, format = "f"))
  shown <- data.frame(drug = d$drug_name, metric = d$metric, a = fmt(d$model_a), b = fmt(d$model_b),
                      difference = fmt(d$difference),
                      interval = ifelse(is.na(d$lower), "", paste0("[", fmt(d$lower), ", ", fmt(d$upper), "]")),
                      matched = d$n_matched, stringsAsFactors = FALSE)
  names(shown)[3:4] <- s$models
  print(shown, row.names = FALSE, right = FALSE)
  if (any(d$n_unmatched > 0)) {
    cat("i ", sum(unique(d[, c("drug", "n_unmatched")])$n_unmatched),
        " isolate(s) predicted by only one model were excluded.\n", sep = "")
  }
  p <- d$mcnemar_p[!is.na(d$mcnemar_p)]
  if (length(p) > 0) {
    cat("McNemar exact p (accuracy, assumes independent isolates): ",
        paste(paste0(d$drug_name[!is.na(d$mcnemar_p)], " ", formatC(p, digits = 3, format = "g")), collapse = "; "),
        "\n", sep = "")
  }
  if (any(d$few_groups)) {
    cat("! Fewer than ", s$min_groups, " independent groups in some strata; intervals may be too narrow.\n",
        sep = "")
  }
  invisible(x)
}

# Helpers -------------------------------------------------------------------

# Metrics on 0/1 vectors; y = reference (1 = R), yhat = predicted class,
# p = predicted probability. Returns a named numeric vector.
fast_metrics <- function(y, yhat = NULL, p = NULL) {
  out <- c()
  if (!is.null(yhat)) {
    tp <- sum(y == 1 & yhat == 1); fn <- sum(y == 1 & yhat == 0)
    tn <- sum(y == 0 & yhat == 0); fp <- sum(y == 0 & yhat == 1)
    sens <- safe_div(tp, tp + fn); spec <- safe_div(tn, tn + fp)
    ppv <- safe_div(tp, tp + fp); npv <- safe_div(tn, tn + fn)
    mcc_den <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    out <- c(sensitivity = sens, specificity = spec, ppv = ppv, npv = npv,
             accuracy = safe_div(tp + tn, length(y)),
             balanced_accuracy = (sens + spec) / 2,
             f1 = safe_div(2 * ppv * sens, ppv + sens),
             mcc = safe_div(as.numeric(tp) * tn - as.numeric(fp) * fn, mcc_den),
             very_major_error = safe_div(fn, tp + fn),
             major_error = safe_div(fp, tn + fp))
  }
  if (!is.null(p)) {
    m <- sum(y == 1); n <- sum(y == 0)
    auc <- if (m > 0 && n > 0) (sum(rank(p)[y == 1]) - m * (m + 1) / 2) / (m * n) else NA_real_
    pc <- pmin(pmax(p, 1e-15), 1 - 1e-15)
    out <- c(out, roc_auc = auc,
             brier = if (length(y)) mean((p - y)^2) else NA_real_,
             log_loss = if (length(y)) -mean(y * log(pc) + (1 - y) * log(1 - pc)) else NA_real_)
  }
  out
}

group_labels <- function(d, group) {
  g <- as.character(d[[group]])
  miss <- is.na(g) | g == ""
  g[miss] <- paste0(".isolate:", d$isolate_id[miss])
  g
}

select_split <- function(x, split) {
  if (!is.null(split) && "split" %in% names(x)) {
    x <- x[x$split %in% split, , drop = FALSE]
    if (nrow(x) == 0) stop("No rows with split = '", paste(split, collapse = "', '"), "'.", call. = FALSE)
  }
  x
}

check_B <- function(B) {
  if (!is.numeric(B) || length(B) != 1 || B < 10 || B != round(B)) {
    stop("`B` must be a whole number of at least 10.", call. = FALSE)
  }
  invisible(TRUE)
}

# Set a seed if given and return a function that restores the previous state.
local_seed <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  set.seed(seed)
  function() {
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
    invisible(NULL)
  }
}
