# AMRvalid simulation study
#
# Question 1 (optimism): how much does each validation design overstate
#   performance on new lineages?
# Question 2 (coverage): how often do isolate-level and lineage-level
#   bootstrap intervals contain the true performance on new lineages?
#
# "True" performance is the ROC-AUC on a large external set of lineages that
# were never used for training.
#
# Usage (from the package root):
#   Rscript simulation/run_study.R [replicates] [cores] [B]
# Defaults: 200 replicates, all but one core, B = 200 bootstrap replicates.
# Outputs go to simulation/results/.

args <- commandArgs(trailingOnly = TRUE)
replicates <- if (length(args) >= 1) as.integer(args[1]) else 200L
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)
B <- if (length(args) >= 3) as.integer(args[3]) else 200L

if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  library(AMRvalid)
}

out_dir <- file.path("simulation", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

scenarios <- data.frame(
  scenario = c("weak lineage effect", "moderate lineage effect", "strong lineage effect"),
  lineage_sd = c(0.5, 1.5, 2.5),
  stringsAsFactors = FALSE
)
settings <- list(n_lineages = 40, isolates_per_lineage = 25, external_lineages = 200,
                 marker_effect = 3, intercept = -1.5, test_fraction = 0.2, k = 5)

auc_of <- function(x) {
  m <- amr_metrics(x)
  m$estimate[m$metric == "roc_auc"]
}

# Each bootstrap gets its own seed, so results do not depend on the number of
# cores or on the order in which replicates run.
boot_ci <- function(x, group, seed) {
  b <- as.data.frame(amr_bootstrap(x, group = group, B = B, split = NULL, seed = seed))
  b <- b[b$metric == "roc_auc", ]
  c(lower = b$lower, upper = b$upper, n_groups = b$n_groups)
}

one_replicate <- function(r, sc) {
  sim <- amr_simulate(n_lineages = settings$n_lineages,
                      isolates_per_lineage = settings$isolates_per_lineage,
                      external_lineages = settings$external_lineages,
                      marker_effect = settings$marker_effect, lineage_sd = sc$lineage_sd,
                      intercept = settings$intercept, test_fraction = settings$test_fraction,
                      k = settings$k, return_distances = FALSE,
                      seed = 100000 * match(sc$scenario, scenarios$scenario) + r)
  base_seed <- 1e7 * match(sc$scenario, scenarios$scenario) + 1000 * r
  x <- amr_import_predictions(sim$predictions)
  rows <- list()
  for (model in sort(unique(x$model_id))) {
    model_seed <- base_seed + 10 * match(model, sort(unique(x$model_id)))
    internal <- x[x$model_id == model & x$dataset_id == "internal", ]
    external <- x[x$model_id == model & x$dataset_id == "external", ]
    truth <- auc_of(external)
    est <- auc_of(internal)
    ci_iso <- boot_ci(internal, "isolate_id", model_seed + 1)
    ci_lin <- boot_ci(internal, "lineage", model_seed + 2)
    rows[[length(rows) + 1]] <- data.frame(
      scenario = sc$scenario, lineage_sd = sc$lineage_sd, replicate = r,
      predictor = sub("/.*", "", model), design = sub(".*/", "", model),
      internal_auc = est, true_auc = truth, optimism = est - truth,
      iso_lower = ci_iso[["lower"]], iso_upper = ci_iso[["upper"]],
      lin_lower = ci_lin[["lower"]], lin_upper = ci_lin[["upper"]],
      n_test_isolates = nrow(internal), n_test_lineages = ci_lin[["n_groups"]],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

started <- Sys.time()
cat("Running", replicates, "replicates x", nrow(scenarios), "scenarios on", cores, "core(s), B =", B, "\n")
results <- list()
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  run <- function(r) tryCatch(one_replicate(r, sc), error = function(e) {
    message("replicate ", r, " failed: ", conditionMessage(e)); NULL
  })
  res <- if (cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(seq_len(replicates), run, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(replicates), run)
  }
  results[[i]] <- do.call(rbind, res)
  cat(sprintf("  %-24s done (%d rows)\n", sc$scenario, nrow(results[[i]])))
}
reps <- do.call(rbind, results)
reps$iso_covers <- reps$iso_lower <= reps$true_auc & reps$true_auc <= reps$iso_upper
reps$lin_covers <- reps$lin_lower <= reps$true_auc & reps$true_auc <= reps$lin_upper
reps$iso_width <- reps$iso_upper - reps$iso_lower
reps$lin_width <- reps$lin_upper - reps$lin_lower
utils::write.csv(reps, file.path(out_dir, "replicates.csv"), row.names = FALSE)

# Summary -------------------------------------------------------------------
key <- interaction(reps$scenario, reps$predictor, reps$design, drop = TRUE, lex.order = TRUE)
summ <- do.call(rbind, lapply(split(reps, key), function(d) {
  data.frame(
    scenario = d$scenario[1], lineage_sd = d$lineage_sd[1], predictor = d$predictor[1], design = d$design[1],
    replicates = nrow(d),
    mean_internal_auc = mean(d$internal_auc), mean_true_auc = mean(d$true_auc),
    mean_optimism = mean(d$optimism),
    optimism_q05 = unname(stats::quantile(d$optimism, 0.05)), optimism_q95 = unname(stats::quantile(d$optimism, 0.95)),
    coverage_isolate_bootstrap = mean(d$iso_covers, na.rm = TRUE),
    coverage_lineage_bootstrap = mean(d$lin_covers, na.rm = TRUE),
    mean_width_isolate = mean(d$iso_width, na.rm = TRUE),
    mean_width_lineage = mean(d$lin_width, na.rm = TRUE),
    test_lineages = stats::median(d$n_test_lineages),
    stringsAsFactors = FALSE
  )
}))
summ <- summ[order(summ$lineage_sd, summ$predictor, summ$design), ]
rownames(summ) <- NULL
utils::write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)

# Figures -------------------------------------------------------------------
source(file.path("simulation", "make_figures.R"))

info <- c(
  paste("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  paste("Duration (minutes):", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)),
  paste("Replicates per scenario:", replicates), paste("Bootstrap replicates:", B),
  paste("Settings:", paste(names(settings), unlist(settings), sep = " = ", collapse = ", ")),
  paste("AMRvalid version:", as.character(utils::packageVersion("AMRvalid"))),
  paste("R version:", R.version.string)
)
writeLines(info, file.path(out_dir, "run_info.txt"))
print(summ[, c("scenario", "predictor", "design", "mean_optimism", "coverage_isolate_bootstrap",
               "coverage_lineage_bootstrap", "test_lineages")], digits = 3, row.names = FALSE)
