#' Variance partition coefficient calculated at different level-1 values
#'
#' Calculates variance estimates (level-2 Intercept variance) and variance partition coefficients (i.e., intra-class correlation) at selected values of predictor values in two-level linear models with random effects (intercept, slope, and their covariation).
#'
#' @param model Two-level model fitted with lme4. Must include random intercept, slope, and their covariation.
#' @param lvl1.var Character string. Level 1 variable name to which random slope is also estimated.
#' @param lvl1.values Level 1 variable values.
#'
#' @return Data frame of level 2 variance and std.dev. estimates at level 1 variable values, respective VPCs (ICC1s) and group-mean reliabilities (ICC2s) (Bliese, 2000).
#' @references Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning Variation in Multilevel Models. Understanding Statistics, 1(4), 223–231. https://doi.org/10.1207/S15328031US0104_02
#' @references Bliese, P. D. (2000). Within-group agreement, non-independence, and reliability: Implications for data aggregation and analysis. In K. J. Klein & S. W. J. Kozlowski (Eds.), Multilevel theory, research, and methods in organizations: Foundations, extensions, and new directions (pp. 349–381). Jossey-Bass.
#' @export
#'
#' @examples
#' fit <- lme4::lmer(Sepal.Length ~ Petal.Length +
#'   (Petal.Length | Species),
#' data = iris
#' )
#'
#' lvl1.values <-
#'   c(
#'     mean(iris$Petal.Length) - stats::sd(iris$Petal.Length),
#'     mean(iris$Petal.Length),
#'     mean(iris$Petal.Length) + stats::sd(iris$Petal.Length)
#'   )
#'
#' vpc_at(
#'   model = fit,
#'   lvl1.var = "Petal.Length",
#'   lvl1.values = lvl1.values
#' )
vpc_at <- function(model, lvl1.var, lvl1.values) {
  VC <- as.data.frame(lme4::VarCorr(model))

  VC.frame <-
    cbind(VC[, c(1:3)],
      sdcor = VC[, 5],
      vcov = VC[, 4]
    )

  VC.frame[is.na(VC.frame)] <- "empty"

  # obtain the covariance estimates from the VC.frame

  # intercept variance
  var_u0 <-
    VC.frame[
      VC.frame[, "var1"] ==
        "(Intercept)" & VC.frame[, "var2"] == "empty",
      "vcov"
    ]

  # covariance between intercept and slope
  cov_u01 <-
    VC.frame[VC.frame[, "var1"] == "(Intercept)" &
      VC.frame[, "var2"] == lvl1.var, "vcov"]

  # if not available, code as zer0
  if (length(cov_u01) == 0) {
    cov_u01 <- 0
  } else {
    cov_u01 <- cov_u01
  }

  # variance of slope
  var_u1 <- VC.frame[VC.frame[, "var1"] ==
    lvl1.var & VC.frame[, "var2"] ==
    "empty", "vcov"]

  # placeholder vector for conditional lvl2 variances
  cond.lvl2.values <- rep(NA, length(lvl1.values))

  for (i in 1:length(lvl1.values)) {
    cond.lvl2.values[i] <- var_u0 +
      2 * cov_u01 * lvl1.values[i] +
      var_u1 * (lvl1.values[i])^2
  }

  output <-
    cbind.data.frame(
      lvl1.value = lvl1.values,
      Intercept.var = cond.lvl2.values,
      Intercept.sd = sqrt(cond.lvl2.values)
    )

  output$Residual.var <-
    VC.frame[VC.frame[
      ,
      "grp"
    ] == "Residual", "vcov"]
  output$Total.var <-
    output$Intercept.var + output$Residual.var

  output$VPC <-
    output$Intercept.var / output$Total.var

  # extract the data used in the model

  data <- model@frame

  # only calculate ICC2 if there are at most 5 different values for level-1 predictor
  # and these values are represented in the data many times
  if (all(lvl1.values %in% names(table(data[, lvl1.var]))) &
    length(table(data[, lvl1.var])) < 6) {
    n.groups <- unname(lme4::getME(model, "l_i"))
    n.obs <- table(data[, lvl1.var])

    # limit to the requested "at" values
    n.obs <- n.obs[lvl1.values %in% names(n.obs)]

    mean.n.obs <- n.obs / n.groups

    output$mean.n.obs <-
      as.numeric(unname(mean.n.obs))

    output$ICC2 <-
      (output$VPC * output$mean.n.obs) / (output$VPC * output$mean.n.obs +
        output$Residual.var)
    output$ICC2_SP <-
      (output$VPC * output$mean.n.obs) /
        (1 + (output$mean.n.obs - 1) * output$VPC)
  }

  return(output)
}
