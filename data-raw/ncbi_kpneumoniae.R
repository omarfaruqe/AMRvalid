# Build the example dataset inst/extdata/ncbi_kpneumoniae.csv.gz
#
# Source: NCBI Pathogen Detection, Klebsiella, PDG000000012.2568
#   Metadata: AST phenotypes (as submitted to BioSample) and AMRFinderPlus
#             genotypes for each isolate
#   Clusters: SNP cluster (PDS) membership
#
# The script downloads the two files into data-raw/cache/ if they are not
# there yet, keeps Klebsiella pneumoniae isolates with a meropenem or
# ciprofloxacin phenotype, applies two transparent rule-based predictors per
# antimicrobial to the AMRFinderPlus calls, and writes a small table in the
# AMRvalid data contract.
#
# Run from the package root:
#   Rscript data-raw/ncbi_kpneumoniae.R
#
# The rules are illustrative textbook rules for demonstrating evaluation
# methods. They are not recommendations for genotype-based susceptibility
# prediction.

pdg <- "PDG000000012.2568"
base_url <- "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/latest_snps"
cache <- file.path("data-raw", "cache")
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
meta_file <- Sys.getenv("AMRVALID_NCBI_METADATA", file.path(cache, "kpneumoniae.metadata.tsv"))
clus_file <- Sys.getenv("AMRVALID_NCBI_CLUSTERS", file.path(cache, "kpneumoniae.cluster_list.tsv"))
if (!file.exists(meta_file)) {
  options(timeout = 3600)
  utils::download.file(file.path(base_url, "Metadata", paste0(pdg, ".metadata.tsv")), meta_file, mode = "wb")
}
if (!file.exists(clus_file)) {
  utils::download.file(file.path(base_url, "Clusters", paste0(pdg, ".reference_target.cluster_list.tsv")),
                       clus_file, mode = "wb")
}

# Read --------------------------------------------------------------------------
cols <- c("target_acc", "biosample_acc", "bioproject_acc", "scientific_name", "collection_date",
          "geo_loc_name", "AST_phenotypes", "AMR_genotypes", "amrfinder_version", "refgene_db_version")
meta <- utils::read.delim(meta_file, quote = "", comment.char = "", na.strings = c("", "NULL"),
                          colClasses = "character", check.names = FALSE)
names(meta) <- sub("^#", "", names(meta))
meta <- meta[!is.na(meta$AST_phenotypes) &
               meta$scientific_name %in% c("Klebsiella pneumoniae", "Klebsiella pneumoniae subsp. pneumoniae"),
             cols]
clusters <- utils::read.delim(clus_file, colClasses = "character")

# Phenotypes -------------------------------------------------------------------------
# AST_phenotypes looks like "amikacin=S,ciprofloxacin=R,..."
phenotype <- function(ast, drug) {
  vapply(strsplit(ast, ",", fixed = TRUE), function(items) {
    hit <- items[startsWith(items, paste0(drug, "="))]
    if (length(hit) == 0) return(NA_character_)
    v <- sub("^[^=]+=", "", hit[1])
    if (v %in% c("S", "I", "R", "SDD")) v else NA_character_
  }, character(1))
}

# Genotype rules ---------------------------------------------------------------------
# AMR_genotypes looks like "blaKPC-2,gyrA_S83I=POINT,ompK35_E42RfsTer47=POINT,...".
# Calls flagged MISTRANSLATION are ignored; PARTIAL calls count as present.
tokens <- strsplit(ifelse(is.na(meta$AMR_genotypes), "", meta$AMR_genotypes), ",", fixed = TRUE)
tokens <- lapply(tokens, function(t) t[!grepl("=MISTRANSLATION$", t)])
has <- function(pattern) vapply(tokens, function(t) any(grepl(pattern, t, perl = TRUE)), logical(1))

carbapenemase <- has("^bla(KPC|NDM|VIM|IMP)-[0-9]+|^blaOXA-(48|162|163|181|204|232|244|245|370|484|515|519|535|547)(\\b|=|$)")
porin <- has("^ompK3[56]_.*=POINT$")
esbl_ampc <- has("^bla(CTX-M|DHA|CMY|ACT|FOX)-[0-9]+")
qrdr <- has("^(gyrA_(S83|D87)|parC_(S80|E84))[A-Z]+=POINT$")
pmqr <- has("^(qnr[A-Z]+[0-9]*(-[0-9]+)?|aac\\(6'\\)-Ib-cr[0-9]*|qepA[0-9]*)")

rules <- list(
  meropenem = list(
    "carbapenemase genes" = carbapenemase,
    "carbapenemase or porin + ESBL/AmpC" = carbapenemase | (porin & esbl_ampc)
  ),
  ciprofloxacin = list(
    "QRDR mutations" = qrdr,
    "QRDR mutations or PMQR genes" = qrdr | pmqr
  )
)

# Assemble ---------------------------------------------------------------------------
year <- suppressWarnings(as.integer(substr(meta$collection_date, 1, 4)))
full_date <- ifelse(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", meta$collection_date), meta$collection_date, NA)
country <- sub(":.*$", "", meta$geo_loc_name)
cluster <- clusters$PDS_acc[match(meta$target_acc, clusters$target_acc)]

out <- list()
for (drug in names(rules)) {
  truth <- phenotype(meta$AST_phenotypes, drug)
  keep <- !is.na(truth)
  for (model in names(rules[[drug]])) {
    out[[length(out) + 1]] <- data.frame(
      isolate_id = meta$biosample_acc[keep],
      species = "Klebsiella pneumoniae",
      drug = drug,
      truth = truth[keep],
      prediction = ifelse(rules[[drug]][[model]][keep], "R", "S"),
      model_id = model,
      dataset_id = paste("NCBI Pathogen Detection", pdg),
      genomic_cluster = cluster[keep],
      bioproject = meta$bioproject_acc[keep],
      collection_year = year[keep],
      collection_date = full_date[keep],
      country = country[keep],
      amrfinder_version = meta$amrfinder_version[keep],
      stringsAsFactors = FALSE
    )
  }
}
ex <- do.call(rbind, out)
ex <- ex[order(ex$drug, ex$model_id, ex$isolate_id), ]
rownames(ex) <- NULL

dir.create(file.path("inst", "extdata"), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(file.path("inst", "extdata", "ncbi_kpneumoniae.csv.gz"), "w")
utils::write.csv(ex, con, row.names = FALSE, na = "")
close(con)

cat("Isolates:", length(unique(ex$isolate_id)), "| rows:", nrow(ex), "\n")
print(table(ex$drug, ex$truth))
cat("SNP clusters:", length(unique(stats::na.omit(ex$genomic_cluster))),
    "| isolates not in a cluster:", sum(is.na(ex$genomic_cluster[!duplicated(ex$isolate_id)])), "\n")
