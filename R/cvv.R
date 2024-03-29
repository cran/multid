#' Coefficient of variance variation
#'
#' Calculates three different indices for variation between two or more variance estimates.
#' VR = Variance ratio between the largest and the smallest variance.
#' CVV = Coefficient of variance variation (Box, 1954).
#' SVH = Standardized variance heterogeneity (Ruscio & Roche, 2012).
#'
#' @param data Data frame of two or more columns or list of two or more variables.
#'
#' @return A vector including VR, CVV, and SVH.
#' @references Box, G. E. P. (1954). Some Theorems on Quadratic Forms Applied in the Study of Analysis of Variance Problems, I. Effect of Inequality of Variance in the One-Way Classification. The Annals of Mathematical Statistics, 25(2), 290–302.
#' @references Ruscio, J., & Roche, B. (2012). Variance Heterogeneity in Published Psychological Research: A Review and a New Index. Methodology, 8(1), 1–11. https://doi.org/10.1027/1614-2241/a000034
#' @export
#'
#' @examples
#' d <- list(
#'   X1 = rnorm(10, sd = 10),
#'   X2 = rnorm(100, sd = 7.34),
#'   X3 = rnorm(1000, sd = 6.02),
#'   X4 = rnorm(100, sd = 5.17),
#'   X5 = rnorm(10, sd = 4.56)
#' )
#' cvv(d)
cvv <- function(data) {

  # change to list
  data <- as.list(data)

  # obtain basic parameters
  k <- length(data)
  n <- sapply(data, length)
  df <- n - 1

  # obtain sample variances
  vars <- sapply(data, stats::var)

  # calculate CVV
  s2_p <- sum(df * vars) / (sum(df))

  numerator <- sum(df * (vars - s2_p)^2)
  denominator <- (sum(n) - k)
  CVV <- sqrt(numerator / denominator) / s2_p

  # calculate SVH
  adj.vars <- k * vars / sum(vars)
  mean.adj.vars <- mean(adj.vars)
  devs.adj.var <- (adj.vars - mean.adj.vars)^2
  adj.vars.sd <- sqrt(sum(devs.adj.var) / k)
  SVH <- adj.vars.sd / sqrt(k - 1)

  # calculate VR
  VR <- max(vars) / min(vars)

  output <- c(VR = VR, CVV = CVV, SVH = SVH)
  return(output)
}
