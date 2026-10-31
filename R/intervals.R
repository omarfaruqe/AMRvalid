#' Confidence intervals for a proportion
#'
#' Wilson score and Clopper-Pearson (exact) intervals. Both behave well
#' when the denominator is small or the proportion is near 0 or 1, as is
#' common for very major error rates where few isolates are resistant.
#'
#' @param x Number of events (numerator).
#' @param n Number of trials (denominator).
#' @param method `"wilson"` (default) or `"clopper-pearson"`.
#' @param conf Confidence level.
#'
#' @return A data frame with `estimate`, `lower` and `upper`. Rows with
#'   `n == 0` give `NA`.
#' @references
#' Wilson EB (1927). Probable inference, the law of succession, and
#' statistical inference. *Journal of the American Statistical Association*
#' 22:209-212.
#'
#' Clopper CJ, Pearson ES (1934). The use of confidence or fiducial limits
#' illustrated in the case of the binomial. *Biometrika* 26:404-413.
#' @export
#' @examples
#' # 1 very major error among 12 resistant isolates
#' amr_prop_ci(1, 12)
#' amr_prop_ci(1, 12, method = "clopper-pearson")
amr_prop_ci <- function(x, n, method = c("wilson", "clopper-pearson"), conf = 0.95) {
  method <- match.arg(method)
  if (length(conf) != 1 || conf <= 0 || conf >= 1) {
    stop("`conf` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (any(x < 0, na.rm = TRUE) || any(x > n, na.rm = TRUE)) {
    stop("`x` must lie between 0 and `n`.", call. = FALSE)
  }
  len <- max(length(x), length(n))
  x <- rep_len(as.numeric(x), len)
  n <- rep_len(as.numeric(n), len)
  est <- ifelse(n > 0, x / n, NA_real_)
  alpha <- 1 - conf
  if (method == "wilson") {
    z <- stats::qnorm(1 - alpha / 2)
    denom <- 1 + z^2 / n
    centre <- (est + z^2 / (2 * n)) / denom
    half <- z * sqrt(est * (1 - est) / n + z^2 / (4 * n^2)) / denom
    lower <- pmax(0, centre - half)
    upper <- pmin(1, centre + half)
  } else {
    lower <- ifelse(x == 0, 0, stats::qbeta(alpha / 2, x, n - x + 1))
    upper <- ifelse(x == n, 1, stats::qbeta(1 - alpha / 2, x + 1, n - x))
  }
  bad <- is.na(n) | n == 0
  lower[bad] <- NA_real_
  upper[bad] <- NA_real_
  data.frame(estimate = est, lower = lower, upper = upper)
}
