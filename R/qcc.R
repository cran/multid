#' Quantile correlation coefficient
#'
#' For computation of tail dependence as correlations estimated at different variable quantiles (Choi & Shin, 2022; Lee et al., 2022) summarized across two quantile regression models where x and y switch roles as independent/dependent variables.
#'
#' @param data Data frame.
#' @param x Name of x variable. Character string.
#' @param y Name of y variable. Character string.
#' @param tau The quantile(s) to be estimated. A vector of values between 0 and 1, default c(.1,.5,.9). @seealso \code{\link[quantreg]{rq}}
#' @param method The algorithmic method used to compute the fit (default "br"). @seealso \code{\link[quantreg]{rq}}
#' @param boot_n Number of bootstrap redraws (default NULL = no bootstrap inference).
#' @param ci_level Level for percentile bootstrap confidence interval. Numeric values between 0 and 1. Default .95.
#'
#' @return
#' \item{r}{Pearson's correlation estimate for comparison.}
#' \item{rho_tau}{Correlations at different tau values (quantiles).}
#' \item{r_boot_est}{Pearson's correlation bootstrap estimates.}
#' \item{rho_tau_boot_est}{Bootstrap estimates for correlations at different tau values (quantiles).}
#' @references Choi, J.-E., & Shin, D. W. (2022). Quantile correlation coefficient: A new tail dependence measure. Statistical Papers, 63(4), 1075–1104. https://doi.org/10.1007/s00362-021-01268-7
#' @references Lee, J. A., Bardi, A., Gerrans, P., Sneddon, J., van Herk, H., Evers, U., & Schwartz, S. (2022). Are value–behavior relations stronger than previously thought? It depends on value importance. European Journal of Personality, 36(2), 133–148. https://doi.org/10.1177/08902070211002965
#' @export
#'
#' @examples
#' set.seed(2321)
#' d <- data.frame(x = rnorm(2000))
#' d$y <- 0.10 * d$x + (0.20) * d$x^2 + 0.40 * d$x^3 + (-0.20) * d$x^4 + rnorm(2000)
#' qcc_boot <- qcc(x = "x", y = "y", data = d, tau = 1:9 / 10, boot_n = 50)
#' qcc_boot$rho_tau
qcc <- function(x, y, tau = c(.1, .5, .9), data, method = "br", boot_n = NULL, ci_level = .95) {
  rq1 <- quantreg::rq(stats::as.formula(paste(y, "~", x)),
    data = data, tau = tau, method = method
  )

  b1 <- rq1$coefficients[rownames(rq1$coefficients) == x, ]

  rq2 <- quantreg::rq(stats::as.formula(paste(x, "~", y)),
    data = data, tau = tau, method = method
  )

  b2 <- rq2$coefficients[rownames(rq2$coefficients) == y, ]

  rho_tau <- sign(b1) * sqrt(b1 * b2)

  r <- stats::cor(data[, x], data[, y], method = "pearson")

  output <-
    list(
      r = r,
      rho_tau = rho_tau
    )

  if (is.numeric(boot_n)) {
    boot.list <- list()

    for (i in 1:boot_n) {
      # sample with replacement
      temp.d <- dplyr::sample_n(data, nrow(data), replace = TRUE)

      # save results of the reanalysis
      boot.list[[i]] <-
        qcc(
          data = temp.d,
          x = x,
          y = y,
          tau = tau,
          method = method
        )
    }

    r_boots <-
      data.frame(r = do.call(rbind, lapply(boot.list, "[[", 1)))

    r_pct_CI <-
      t(sapply(
        data.frame(r_boots),
        function(x) {
          stats::quantile(
            x,
            c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2),
            na.rm = TRUE
          )
        }
      ))

    rho_tau_boots <-
      data.frame(do.call(rbind, lapply(boot.list, "[[", 2)))

    rho_tau_pct_CI <-
      t(sapply(
        data.frame(rho_tau_boots),
        function(x) {
          stats::quantile(
            x,
            c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2),
            na.rm = TRUE
          )
        }
      ))

    output <-
      list(
        r = cbind(r, r_pct_CI),
        rho_tau = cbind(rho = rho_tau, rho_tau_pct_CI),
        r_boots = r_boots,
        rho_tau_boots = rho_tau_boots
      )
  }

  return(output)
}
