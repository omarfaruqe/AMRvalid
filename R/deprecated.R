#' Deprecated functions
#'
#' These functions still work but will be removed in a future release. Use
#' the replacements listed below.
#'
#' * `validate_amr_data()`: use [amr_validate()].
#' * `amr_species_eval()`: use [amr_evaluate()].
#'
#' @param ... Arguments passed to the replacement function.
#' @return The result of the replacement function.
#' @name AMRvalid-deprecated
#' @keywords internal
NULL

#' @rdname AMRvalid-deprecated
#' @export
validate_amr_data <- function(...) {
  .Deprecated("amr_validate", package = "AMRvalid")
  amr_validate(...)
}

#' @rdname AMRvalid-deprecated
#' @export
amr_species_eval <- function(...) {
  .Deprecated("amr_evaluate", package = "AMRvalid")
  amr_evaluate(...)
}
