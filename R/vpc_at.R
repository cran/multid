#' Variance partition coefficient calculated at different level-1 values
#'
#' Calculates variance estimates (level-2 Intercept variance) and variance partition coefficients (i.e., intra-class correlation) at selected values of predictor values in two-level linear models with random effects (intercept, slope, and their covariation).
#'
#' @param model Two-level model fitted with lme4. Must include random intercept, slope, and their covariation.
#' @param lvl1.var Character string. Level 1 variable name to which random slope is also estimated.
#' @param lvl1.values Level 1 variable values.
#'
#' @return Data frame of level 2 variance and std.dev. estimates at level 1 variable values and respective VPCs.
#' @references Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning Variation in Multilevel Models. Understanding Statistics, 1(4), 223–231. https://doi.org/10.1207/S15328031US0104_02
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

  # rename <NA> in the frame as "empty
  VC.frame[is.na(VC.frame)] <- "empty"

  # placeholder vector for variances at each value
  cond.lvl2.values <- rep(NA, length(lvl1.values))

  for (i in 1:length(lvl1.values)) {
    cond.lvl2.values[i] <-
      VC.frame[VC.frame[, "var1"] == "(Intercept)" &
        VC.frame[, "var2"] == "empty", "vcov"] +
      2 * VC.frame[VC.frame[, "var1"] == "(Intercept)" &
        VC.frame[, "var2"] == lvl1.var, "vcov"] * (lvl1.values[i]) +
      VC.frame[VC.frame[, "var1"] == lvl1.var &
        VC.frame[, "var2"] == "empty", "vcov"] * (lvl1.values[i])^2
  }

  output <- cbind.data.frame(
    lvl1.value = lvl1.values,
    Intercept.var = cond.lvl2.values,
    Intercept.sd = sqrt(cond.lvl2.values)
  )
  output$Total.var <- output$Intercept.var +
    VC.frame[VC.frame[, "grp"] == "Residual", "vcov"]
  output$VPC <- output$Intercept.var / output$Total.var

  return(output)
}
