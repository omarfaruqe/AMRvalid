report_inputs <- function() {
  sim <- amr_simulate(n_lineages = 20, isolates_per_lineage = 10, external_lineages = 10, seed = 4)
  x <- amr_import_predictions(sim$predictions)
  list(sim = sim, x = x)
}

test_that("provenance labels items as checked, supplied or not assessable", {
  r <- report_inputs()
  e <- amr_evaluate(r$x, threshold = 0.5)
  p <- amr_provenance(r$x, e, model_info = list(threshold_selection = "fixed at 0.5"))
  expect_s3_class(p, "amr_provenance")
  expect_true(all(p$status %in% c("checked", "supplied", "not assessable")))
  expect_equal(p$status[p$item == "How the threshold was chosen"], "supplied")
  expect_equal(p$value[p$item == "Classification threshold"], "0.5")
  expect_equal(p$status[p$item == "Data source"], "not assessable")
  expect_equal(p$status[p$item == "Group-aware uncertainty"], "not assessable")
  expect_true(all(c("5a", "8a", "12a", "12e", "18f") %in% p$tripod_item))
  expect_output(print(p), "not assessable")
})

test_that("provenance records split checks, bootstrap and calibration", {
  r <- report_inputs()
  x <- r$x[r$x$dataset_id == "internal", ]
  sc <- amr_validate_splits(x, train_manifest = r$sim$train_manifest, group = "lineage",
                            claimed_design = "group")
  b <- amr_bootstrap(x[x$model_id == "marker/lineage", ], resample_by = "lineage", B = 20, seed = 1)
  cal <- amr_calibration(x)
  p <- amr_provenance(x, sc, b, cal)
  expect_match(p$value[p$item == "Claimed validation design"], "NOT consistent")
  expect_match(p$value[p$item == "Group-aware uncertainty"], "lineage, B = 20, seed = 1")
  expect_match(p$value[p$item == "Performance measures assessed"], "calibration")
  expect_error(amr_provenance(x, data.frame(a = 1)), "Unsupported")
})

test_that("report contains every supplied section and the TRIPOD+AI table", {
  r <- report_inputs()
  x <- r$x[r$x$dataset_id == "internal", ]
  e <- amr_evaluate(x, threshold = 0.5)
  sc <- amr_validate_splits(x, train_manifest = r$sim$train_manifest, group = "lineage")
  cal <- amr_calibration(x)
  d <- amr_distance_eval(x, r$sim$distances, train_manifest = r$sim$train_manifest)
  b <- amr_bootstrap(x, resample_by = "lineage", B = 20, threshold = 0.5, seed = 1)
  txt <- amr_report(x, e, sc, cal, d, b, dataset_info = list(source = "simulation"))
  for (heading in c("## Data", "## Validation checks", "## Split integrity", "## Performance by stratum",
                    "## Calibration", "## Performance by distance to training", "## Group-aware intervals",
                    "## Limitations detected", "## TRIPOD+AI coverage")) {
    expect_true(grepl(heading, txt, fixed = TRUE), info = heading)
  }
  expect_match(txt, "| 5a | Data | Data source | simulation | supplied |", fixed = TRUE)
  expect_match(txt, "Fewer than 20 independent groups")
})

test_that("report writes Markdown and HTML files", {
  r <- report_inputs()
  e <- amr_evaluate(r$x, threshold = 0.5)
  md <- amr_report(r$x, e, file = tempfile(fileext = ".md"))
  expect_true(file.exists(md))
  expect_match(readLines(md, n = 1), "^# AMRvalid evaluation report")
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  html <- amr_report(r$x, e, file = tempfile(fileext = ".html"), format = "html")
  expect_true(any(grepl("<table", readLines(html))))
})

test_that("markdown tables escape pipes and format numbers", {
  t <- AMRvalid:::md_table(data.frame(a = c("x|y", NA), b = c(1, 0.12345)))
  expect_equal(t[3], "| x\\|y | 1 |")
  expect_equal(t[4], "|  | 0.123 |")
})

test_that("TRIPOD+AI rows are sorted in checklist order", {
  p <- data.frame(tripod_item = c("12e", "5a", "20a", "12a", "5b"))
  expect_equal(AMRvalid:::tripod_order(p)$tripod_item, c("5a", "5b", "12a", "12e", "20a"))
})
