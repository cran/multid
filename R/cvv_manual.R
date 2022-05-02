#' Coefficient of variance variation from manual input sample sizes and variance estimates
#'
#' Calculates three different indices for variation between two or more variance estimates.
#' VR = Variance ratio between the largest and the smallest variance.
#' CVV = Coefficient of variance variation (Box, 1954).
#' SVH = Standardized variance heterogeneity (Ruscio & Roche, 2012).
#'
#' @param sample_sizes Numeric vector of length > 1. Sample sizes used for each variance estimate.
#' @param variances Numeric vector of length > 1. Variance estimates.
#'
#' @return A vector including VR, CVV, and SVH.
#' @references Box, G. E. P. (1954). Some Theorems on Quadratic Forms Applied in the Study of Analysis of Variance Problems, I. Effect of Inequality of Variance in the One-Way Classification. The Annals of Mathematical Statistics, 25(2), 290–302.
#' @references Ruscio, J., & Roche, B. (2012). Variance Heterogeneity in Published Psychological Research: A Review and a New Index. Methodology, 8(1), 1–11. https://doi.org/10.1027/1614-2241/a000034
#' @export
#'
#' @examples
#' cvv_manual(sample_sizes=c(10,100,1000,75,3),
#' variances=c(1.5,2,2.5,3,3.5))
cvv_manual <- function(sample_sizes,variances) {

  if(length(sample_sizes)!=length(variances)){
    stop("Unequal number of sample sizes and variance estimates")
  }

  # obtain basic parameters
  n <- sample_sizes
  k <- length(n)

  df <- n - 1

  # obtain sample variances
  vars <- variances


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
