#' The AMRvalid data contract
#'
#' Returns the field specification that every prediction table must follow
#' after import. One row is one prediction for one isolate, one antimicrobial,
#' one model and one evaluation context (dataset, split and fold).
#'
#' @return A data frame with columns `field`, `requirement`, `type` and
#'   `description`.
#' @export
#' @examples
#' amr_contract()
amr_contract <- function() {
  fields <- list(
    c("isolate_id",          "required",    "character", "Stable isolate identifier."),
    c("species",             "required",    "mo",        "Organism, coerced to the AMR package 'mo' class."),
    c("drug",                "required",    "ab",        "Antimicrobial, coerced to the AMR package 'ab' class."),
    c("truth_sir",           "conditional", "sir",       "Reference category (S, SDD, I, R, NI). Required unless truth_mic is supplied."),
    c("truth_mic",           "conditional", "mic",       "Reference MIC. Required unless truth_sir is supplied."),
    c("prediction_type",     "required",    "character", "One of 'binary', 'score', 'sir' or 'mic'; detected at import if absent."),
    c("probability",         "conditional", "numeric",   "Predicted probability of resistance in [0, 1]."),
    c("prediction",          "conditional", "sir",       "Predicted category (binary S/R or S/I/R)."),
    c("predicted_mic",       "conditional", "mic",       "Predicted MIC, for MIC-regression models."),
    c("model_id",            "recommended", "character", "External model or software identifier and version."),
    c("dataset_id",          "recommended", "character", "Dataset or cohort identifier."),
    c("split",               "recommended", "character", "One of 'train', 'validation', 'calibration' or 'test'."),
    c("fold_id",             "optional",    "character", "Cross-validation fold identifier."),
    c("lineage",             "optional",    "character", "Lineage label."),
    c("sequence_type",       "optional",    "character", "MLST sequence type."),
    c("genomic_cluster",     "optional",    "character", "Population cluster identifier (e.g. from PopPUNK)."),
    c("country",             "optional",    "character", "Geographic origin."),
    c("collection_date",     "optional",    "Date",      "Isolate collection date."),
    c("breakpoint_standard", "recommended", "character", "Breakpoint authority and version, e.g. 'EUCAST 2025' or 'CLSI 2025'."),
    c("ast_method",          "optional",    "character", "Reference AST method, e.g. 'broth microdilution'.")
  )
  out <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE)
  names(out) <- c("field", "requirement", "type", "description")
  out
}

#' Path to the JSON Schema of the data contract
#'
#' The schema lets workflows in other languages (for example Python) export
#' prediction tables that AMRvalid can import without changes.
#'
#' @return Path to `amrvalid-contract.json` inside the installed package.
#' @export
#' @examples
#' amr_contract_schema()
amr_contract_schema <- function() {
  system.file("schema", "amrvalid-contract.json", package = "AMRvalid", mustWork = TRUE)
}

# Columns that, together, identify one prediction record.
key_columns <- function(x) {
  intersect(c("isolate_id", "species", "drug", "model_id", "dataset_id", "split", "fold_id"),
            names(x))
}

# Columns that identify one evaluation stratum.
stratum_columns <- function(x) {
  intersect(c("model_id", "dataset_id", "species", "drug"), names(x))
}

valid_splits <- c("train", "validation", "calibration", "test")
valid_prediction_types <- c("binary", "score", "sir", "mic")
