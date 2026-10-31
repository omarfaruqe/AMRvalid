#' Simulate AMR predictions with known population structure
#'
#' Generates isolates in lineages, a resistance phenotype driven by a
#' causal marker plus lineage-inherited background, genomic distances, and
#' predictions from two fixed predictors evaluated under a random and a
#' lineage-held-out split. It exists to teach and test evaluation methods:
#' because the generating process is known, it shows how much each
#' validation design overstates performance on new lineages.
#'
#' @param n_lineages Number of lineages in the development dataset.
#' @param isolates_per_lineage Isolates per lineage.
#' @param external_lineages Number of fresh lineages in an external dataset
#'   (`dataset_id = "external"`), never seen in training. Large values give
#'   a close approximation of true performance on new lineages.
#' @param marker_effect Log-odds effect of the causal marker.
#' @param lineage_sd Standard deviation of the lineage effect on the
#'   log-odds of resistance (resistance inherited clonally but not
#'   explained by the marker).
#' @param intercept Baseline log-odds of resistance.
#' @param within_sd Spread of isolates around their lineage centre, in
#'   distance units (lineage centres are spread over 0 to 1000).
#' @param test_fraction Fraction of isolates (random split) or lineages
#'   (lineage-held-out split) placed in the test set.
#' @param k Number of neighbours for the nearest-neighbour predictor.
#' @param return_distances If `TRUE`, include the full pairwise distance
#'   matrix (memory grows with the square of the number of isolates).
#' @param seed Optional random seed; the global random state is restored.
#'
#' @details
#' Each lineage has a centre in a two-dimensional space; isolate distances
#' are Euclidean distances between isolates, standing in for genomic
#' distances. Marker frequency varies by lineage (Beta(0.5, 0.5)), so the
#' marker is partly confounded with lineage, as resistance genes often are.
#'
#' The two predictors are fixed rules, not trained machine-learning models:
#'
#' * `knn`: the proportion of resistant isolates among the `k` nearest
#'   training isolates. It memorises lineages.
#' * `marker`: the proportion resistant among training isolates with the
#'   same marker status (with a +1/+2 smoothing). It learns the mechanism.
#'
#' Each predictor is applied under two designs, giving model identifiers
#' `knn/random`, `knn/lineage`, `marker/random` and `marker/lineage`.
#'
#' @return A list with `predictions` (a data frame ready for
#'   [amr_import_predictions()]; test isolates only), `train_manifest`
#'   (training isolates per model), `isolates` (all isolates with lineage,
#'   marker, phenotype and coordinates), `distances` (a matrix or `NULL`)
#'   and `parameters`.
#' @export
#' @examples
#' sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10,
#'                     external_lineages = 10, seed = 1)
#' x <- amr_import_predictions(sim$predictions)
#' e <- amr_species_eval(x, threshold = 0.5, validate = FALSE)
#' e[e$metric == "accuracy", c("model_id", "dataset_id", "estimate")]
amr_simulate <- function(n_lineages = 40, isolates_per_lineage = 25, external_lineages = 40,
                         marker_effect = 3, lineage_sd = 1.5, intercept = -1.5, within_sd = 5,
                         test_fraction = 0.2, k = 5, return_distances = TRUE, seed = NULL) {
  stopifnot(n_lineages >= 5, isolates_per_lineage >= 1, external_lineages >= 0,
            test_fraction > 0, test_fraction < 1, k >= 1)
  restore <- local_seed(seed)
  on.exit(restore(), add = TRUE)

  make_lineages <- function(n_lin, prefix) {
    centre <- matrix(stats::runif(n_lin * 2, 0, 1000), ncol = 2)
    freq <- stats::rbeta(n_lin, 0.5, 0.5)
    effect <- stats::rnorm(n_lin, 0, lineage_sd)
    lin <- rep(seq_len(n_lin), each = isolates_per_lineage)
    n <- length(lin)
    marker <- stats::rbinom(n, 1, freq[lin])
    p_true <- stats::plogis(intercept + marker_effect * marker + effect[lin])
    data.frame(
      isolate_id = sprintf("%s%04d", prefix, seq_len(n)),
      lineage = sprintf("%sL%03d", prefix, lin),
      marker = marker,
      p_true = p_true,
      truth = ifelse(stats::rbinom(n, 1, p_true) == 1, "R", "S"),
      x = centre[lin, 1] + stats::rnorm(n, 0, within_sd),
      y = centre[lin, 2] + stats::rnorm(n, 0, within_sd),
      stringsAsFactors = FALSE
    )
  }
  dev <- make_lineages(n_lineages, "D")
  dev$dataset_id <- "internal"
  ext <- if (external_lineages > 0) make_lineages(external_lineages, "E") else dev[0, ]
  if (nrow(ext) > 0) ext$dataset_id <- "external"

  n_test_iso <- max(1, round(test_fraction * nrow(dev)))
  test_random <- dev$isolate_id[sample.int(nrow(dev), n_test_iso)]
  lins <- unique(dev$lineage)
  test_lins <- lins[sample.int(length(lins), max(1, round(test_fraction * length(lins))))]
  test_lineage <- dev$isolate_id[dev$lineage %in% test_lins]
  designs <- list(random = test_random, lineage = test_lineage)

  preds <- list()
  manifest <- list()
  for (design in names(designs)) {
    tr <- dev[!dev$isolate_id %in% designs[[design]], , drop = FALSE]
    te <- rbind(dev[dev$isolate_id %in% designs[[design]], , drop = FALSE], ext)
    tr_r <- as.integer(tr$truth == "R")
    p_knn <- knn_resistance(te, tr, tr_r, k)
    p_marker <- ifelse(te$marker == 1,
                       (sum(tr_r[tr$marker == 1]) + 1) / (sum(tr$marker == 1) + 2),
                       (sum(tr_r[tr$marker == 0]) + 1) / (sum(tr$marker == 0) + 2))
    for (model in c("knn", "marker")) {
      id <- paste0(model, "/", design)
      preds[[id]] <- data.frame(
        isolate_id = te$isolate_id, species = "Klebsiella pneumoniae", drug = "meropenem",
        truth = te$truth, probability = if (model == "knn") p_knn else p_marker,
        model_id = id, dataset_id = te$dataset_id, lineage = te$lineage, split = "test",
        stringsAsFactors = FALSE
      )
      manifest[[id]] <- data.frame(isolate_id = tr$isolate_id, model_id = id, lineage = tr$lineage,
                                   stringsAsFactors = FALSE)
    }
  }
  isolates <- rbind(dev, ext)
  distances <- NULL
  if (return_distances) {
    distances <- as.matrix(stats::dist(isolates[, c("x", "y")]))
    dimnames(distances) <- list(isolates$isolate_id, isolates$isolate_id)
  }
  out_preds <- do.call(rbind, preds)
  rownames(out_preds) <- NULL
  out_manifest <- do.call(rbind, manifest)
  rownames(out_manifest) <- NULL
  list(predictions = out_preds, train_manifest = out_manifest, isolates = isolates,
       distances = distances,
       parameters = list(n_lineages = n_lineages, isolates_per_lineage = isolates_per_lineage,
                         external_lineages = external_lineages, marker_effect = marker_effect,
                         lineage_sd = lineage_sd, intercept = intercept, within_sd = within_sd,
                         test_fraction = test_fraction, k = k, seed = seed))
}

# Proportion resistant among the k nearest training isolates, in blocks to
# keep memory bounded.
knn_resistance <- function(te, tr, tr_r, k) {
  k <- min(k, nrow(tr))
  out <- numeric(nrow(te))
  block <- 500
  for (start in seq(1, nrow(te), by = block)) {
    idx <- start:min(nrow(te), start + block - 1)
    dx <- outer(te$x[idx], tr$x, "-")
    dy <- outer(te$y[idx], tr$y, "-")
    d2 <- dx^2 + dy^2
    out[idx] <- apply(d2, 1, function(v) mean(tr_r[order(v)[seq_len(k)]]))
  }
  out
}
