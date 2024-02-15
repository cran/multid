ci_to_ddsc_ml <- function(results, level) {
  ci.lower <- results[, "estimate"] +
    stats::qt(p = (1 - level) / 2, df = results[, "df"]) * results[, "SE"]
  ci.upper <- results[, "estimate"] +
    stats::qt(p = 1 - (1 - level) / 2, df = results[, "df"]) * results[, "SE"]
  CIs <- cbind(ci.lower, ci.upper)
  return(CIs)
}
