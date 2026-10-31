#' AMRvalid: population-aware evaluation of AMR predictions
#'
#' AMRvalid evaluates antimicrobial resistance (AMR) phenotype predictions
#' produced by external models and tools. It does not train models, detect
#' resistance genes or interpret raw MIC values into categories.
#'
#' The core workflow is:
#' 1. [amr_import_predictions()] maps a table onto the AMRvalid data contract.
#' 2. [validate_amr_data()] reports blocking errors and warnings.
#' 3. [amr_metrics()], [amr_ast_agreement()] and [amr_species_eval()]
#'    compute performance per species-antimicrobial stratum.
#'
#' @keywords internal
"_PACKAGE"
