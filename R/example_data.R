#' Example data: rule-based predictions for Klebsiella pneumoniae
#'
#' Real, public data for trying AMRvalid: laboratory susceptibility results
#' for meropenem and ciprofloxacin in *Klebsiella pneumoniae* isolates from
#' NCBI Pathogen Detection, with predictions from two transparent rule-based
#' predictors per antimicrobial applied to the AMRFinderPlus genotypes.
#'
#' @return A data frame ready for [amr_import_predictions()], with one row
#'   per isolate, antimicrobial and predictor:
#' \describe{
#'   \item{isolate_id}{NCBI BioSample accession.}
#'   \item{species}{Always "Klebsiella pneumoniae".}
#'   \item{drug}{"meropenem" or "ciprofloxacin".}
#'   \item{truth}{Reference category (S, I or R) as submitted to BioSample.}
#'   \item{prediction}{Predicted category (S or R) from the rule.}
#'   \item{model_id}{The rule (see Details).}
#'   \item{dataset_id}{Source and NCBI analysis version.}
#'   \item{genomic_cluster}{NCBI SNP cluster (PDS accession); empty when the
#'     isolate is not in a cluster.}
#'   \item{bioproject}{BioProject accession (the submitting study).}
#'   \item{collection_year, collection_date}{Year, and full date when
#'     available.}
#'   \item{country}{Country of isolation.}
#'   \item{amrfinder_version}{AMRFinderPlus version used by NCBI.}
#' }
#'
#' @details
#' Source: NCBI Pathogen Detection, *Klebsiella* analysis
#' PDG000000012.2568 (metadata and SNP cluster files, downloaded
#' 2 November 2026). Isolates were kept if NCBI labelled them
#' *K. pneumoniae* and they had a meropenem or ciprofloxacin phenotype.
#' Breakpoint standards are those of the submitting laboratories and are not
#' recorded consistently, so the data are suitable for demonstrating
#' methods, not for drawing clinical conclusions.
#'
#' Rules applied to the AMRFinderPlus calls (calls flagged MISTRANSLATION
#' ignored):
#'
#' * **carbapenemase genes** (meropenem): any KPC, NDM, VIM, IMP or
#'   OXA-48-like gene;
#' * **carbapenemase or porin + ESBL/AmpC** (meropenem): the above, or an
#'   OmpK35/OmpK36 alteration together with a CTX-M, DHA, CMY, ACT or FOX
#'   gene;
#' * **QRDR mutations** (ciprofloxacin): gyrA S83/D87 or parC S80/E84
#'   substitutions;
#' * **QRDR mutations or PMQR genes** (ciprofloxacin): the above, or a qnr,
#'   aac(6')-Ib-cr or qepA gene.
#'
#' These are illustrative textbook rules, not recommendations. The script
#' that builds the data, `data-raw/ncbi_kpneumoniae.R`, is in the package
#' source repository.
#'
#' @source NCBI Pathogen Detection,
#'   <https://www.ncbi.nlm.nih.gov/pathogens/>. Data produced by the U.S.
#'   National Center for Biotechnology Information are in the public domain.
#' @export
#' @examples
#' d <- amr_example_data()
#' head(d)
#' x <- amr_import_predictions(d)
#' x
amr_example_data <- function() {
  path <- system.file("extdata", "ncbi_kpneumoniae.csv.gz", package = "AMRvalid", mustWork = TRUE)
  utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "",
                  colClasses = c(collection_year = "integer", collection_date = "character",
                                 genomic_cluster = "character", amrfinder_version = "character"))
}
