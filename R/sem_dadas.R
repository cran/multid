#' Predicting algebraic difference scores in structural equation model
#'
#' @param data A data frame.
#' @param var1 Character string. Variable name of first component score of difference score (Y_1).
#' @param var2 Character string. Variable name of second component score of difference score (Y_2).
#' @param predictor Character string. Variable name of independent variable predicting difference score.
#' @param covariates Character string or vector. Variable names of covariates (Default NULL).
#' @param estimator Character string. Estimator used in SEM (Default "MLR").
#' @param scale Logical. Are var1 and var2 scaled with their pooled sd? (Default FALSE)
#' @param center Logical. Are var1 and var2 centered around their grand mean? (Default FALSE)
#' @param level Numeric. The confidence level required for the result output (Default .95)
#' @param sampling.weights Character string. Name of sampling weights variable.
#'
#' @return
#' \item{descriptives}{Means, standard deviations, and intercorrelations.}
#' \item{parameter_estimates}{Parameter estimates from the structural equation model.}
#' \item{variance_test}{Variances and covariances of component scores.}
#' \item{transformed_data}{Data frame with variables used in SEM.}
#' \item{dadas}{One sided dadas-test for positivity of abs(b_11-b_21)-abs(b_11+b_21).}
#' \item{results}{Summary of key results.}
#' @references Edwards, J. R. (1995). Alternatives to Difference Scores as Dependent Variables in the Study of Congruence in Organizational Research. Organizational Behavior and Human Decision Processes, 64(3), 307–324.
#'
#' @export
#'
#' @examples set.seed(342356)
#' d <- data.frame(
#'   var1 = rnorm(50),
#'   var2 = rnorm(50),
#'   x = rnorm(50)
#' )
#' sem_dadas(
#'   data = d, var1 = "var1", var2 = "var2",
#'   predictor = "x", center = TRUE, scale = TRUE
#' )$results
sem_dadas <- function(data,
                      var1,
                      var2,
                      center = FALSE,
                      scale = FALSE,
                      predictor,
                      covariates = NULL,
                      estimator = "MLR",
                      level = .95,
                      sampling.weights = NULL) {
  if (center) {
    pooled_mean <-
      mean(c(mean(data[, var1]), mean(data[, var2])))

    data[, var1] <- data[, var1] - pooled_mean
    data[, var2] <- data[, var2] - pooled_mean
  }

  if (scale) {
    pooled_sd <-
      sqrt(((nrow(data) - 1) * stats::sd(data[, var1])^2 +
        (nrow(data) - 1) * stats::sd(data[, var2])^2) /
        (nrow(data) * 2 - 2))

    data[, var1] <- data[, var1] / pooled_sd
    data[, var2] <- data[, var2] / pooled_sd
    data[, predictor] <-
      (data[, predictor] - mean(data[, predictor])) /
        stats::sd(data[, predictor])
  }


  data[, "diff"] <- data[, var1] - data[, var2]

  descriptives <-
    cbind(
      rbind(
        c(mean(data[, var1]), stats::sd(data[, var1])),
        c(mean(data[, var2]), stats::sd(data[, var2])),
        c(mean(data[, "diff"]), stats::sd(data[, "diff"])),
        c(mean(data[, predictor]), stats::sd(data[, predictor]))
      ),
      stats::cor(data[, c(var1, var2, "diff", predictor)])
    )

  colnames(descriptives) <-
    c("M", "SD", var1, var2, "diff", predictor)

  model <- paste0(
    paste0(var1, "~b_11*", predictor), "\n",
    paste0(var2, "~b_21*", predictor), "\n",
    paste0(var1, "~b_10*1"), "\n",
    paste0(var2, "~b_20*1"), "\n",
    paste0(var1, "~~rescov_12*", var2), "\n",
    paste0(var1, "~~e_1*", var1), "\n",
    paste0(var2, "~~e_2*", var2), "\n",
    paste0(predictor, "~pred_int*1"), "\n",
    paste0(predictor, "~~pred_var*", predictor), "\n",
    paste0("coef_diff:=b_11-b_21"), "\n",
    paste0(
      "coef_diff_std:=((b_11-b_21)*sqrt(pred_var))/",
      as.character(round(descriptives["diff", "SD"], 8))
    ), "\n",
    paste0("coef_sum:=b_11+b_21"), "\n",
    paste0("diff_abs_magnitude:=sqrt(b_11^2)-sqrt(b_21^2)"), "\n",
    paste0("abs_coef_diff:=sqrt((b_11-b_21)^2)"), "\n",
    paste0("abs_coef_sum:=sqrt((b_11+b_21)^2)"), "\n",
    paste0("dadas_two_sided:=sqrt((b_11-b_21)^2)-sqrt((b_11+b_21)^2)")
  )

  # include covariates to the model


  if (!is.null(covariates)) {
    model <- paste0(model, "\n")
    model <- paste0(
      model,
      paste0(var1, "~", covariates, collapse = "\n")
    )
    model <- paste0(model, "\n")
    model <- paste0(
      model,
      paste0(var2, "~", covariates, collapse = "\n")
    )
    model <- paste0(model, "\n")
    model <- paste0(
      model,
      paste0(predictor, "~~", covariates, collapse = "\n")
    )
  }

  # fit model

  fit <-
    lavaan::sem(
      model = model,
      data = data,
      estimator = estimator,
      sampling.weights = sampling.weights
    )

  output.data <-
    data[, c(var1, var2, "diff", predictor)]

  pars <- data.frame(lavaan::parameterestimates(fit,
    rsquare = TRUE,
    level = level
  ))

  # one-sided tests for absolute parameters

  labels <- c("abs_coef_diff", "abs_coef_sum", "dadas_two_sided")

  osts <- pars[pars$label %in% labels, ]
  osts$p.pos <- stats::pnorm(osts[, "z"], lower.tail = F)

  res.pars <- c(
    "b_11", "b_21", "b_10", "b_20", "rescov_12",
    "coef_diff", "coef_diff_std",
    "coef_sum", "diff_abs_magnitude",
    "abs_coef_diff", "abs_coef_sum",
    "dadas_two_sided"
  )

  results <- pars[pars$label %in% res.pars, 4:ncol(pars)]
  rownames(results) <- results$label
  results <- results[, 2:ncol(results)]

  # replace two-sided tests with one-sided for absolute parameters

  results[labels, "pvalue"] <- osts[, "p.pos"]

  rownames(results) <- c(
    rownames(results)[1:(nrow(results) - 1)],
    "dadas"
  )

  rsquared <- pars[pars[, "op"] == "r2", c(1, 5)]
  rownames(rsquared) <- NULL
  colnames(rsquared) <- c("score", "r2")

  b_11 <- results["b_11", "est"]
  b_21 <- results["b_21", "est"]
  sd_predictor <- descriptives[predictor, "SD"]
  cov_C1C2 <- stats::cov(output.data[, c(var1, var2)])

  r_diff <- (b_11 - b_21) * sd_predictor /
    sqrt(cov_C1C2[var1, var1] + cov_C1C2[var2, var2]
      - 2 * cov_C1C2[var1, var2])
  r2_diff <- r_diff^2

  r2_diff_row <-
    cbind.data.frame(
      score = "difference",
      r2 = r2_diff
    )
  rsquared <-
    rbind(rsquared, r2_diff_row)

  # variance test in a separate model

  var_model <- paste0(
    paste0(var1, "~~cov_12*", var2), "\n",
    paste0(var1, "~~var_1*", var1), "\n",
    paste0(var2, "~~var_2*", var2), "\n",
    paste0("var_diff:=var_1-var_2"), "\n",
    paste0("var_ratio:=var_1/var_2"), "\n",
    paste0("cor_12:=cov_12/(sqrt(var_1)*sqrt(var_2))")
  )

  var_fit <-
    lavaan::sem(
      model = var_model,
      data = data,
      estimator = estimator,
      sampling.weights = sampling.weights
    )

  var_pars <- data.frame(lavaan::parameterestimates(var_fit,
    level = level
  ))

  var_results <- var_pars[, 4:ncol(var_pars)]
  rownames(var_results) <- var_results$label
  var_results <- var_results[, 2:ncol(var_results)]

  output <- list(
    variance_test = var_results,
    descriptives = descriptives,
    parameter_estimates = pars,
    transformed_data = output.data,
    rsquared = rsquared,
    results = results
  )

  return(output)
}
