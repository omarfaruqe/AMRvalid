toy_scores <- function() {
  data.frame(
    isolate_id  = sprintf("iso%02d", 1:10),
    species     = "Klebsiella pneumoniae",
    drug        = "meropenem",
    truth       = c("R", "R", "R", "R", "S", "S", "S", "S", "S", "S"),
    probability = c(0.95, 0.80, 0.60, 0.30, 0.70, 0.40, 0.20, 0.10, 0.05, 0.02),
    stringsAsFactors = FALSE
  )
}
