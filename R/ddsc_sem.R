#' Deconstructing difference score correlation with structural equation modeling
#'
#' @param data A data frame.
#' @param y1 Character string. Variable name of first component score of difference score.
#' @param y2 Character string. Variable name of second component score of difference score.
#' @param x Character string. Variable name of independent variable.
#' @param covariates Character string or vector. Variable names of covariates (Default NULL).
#' @param estimator Character string. Estimator used in SEM (Default "ML").
#' @param center_yvars Logical. Should y1 and y2 be centered around their grand mean? (Default FALSE)
#' @param level Numeric. The confidence level required for the result output (Default .95)
#' @param sampling.weights Character string. Name of sampling weights variable.
#' @param q_sesoi Numeric. The smallest effect size of interest for Cohen's q estimates (Default 0; See Lakens et al. 2018).
#' @param min_cross_over_point_location Numeric. Z-score for the minimal slope cross-over point of interest (Default 0).
#'
#' @return
#' \item{descriptives}{Means, standard deviations, and intercorrelations.}
#' \item{parameter_estimates}{Parameter estimates from the structural equation model.}
#' \item{variance_test}{Variances and covariances of component scores.}
#' \item{data}{Data frame with original and scaled variables used in SEM.}
#' \item{results}{Summary of key results.}
#' @references Edwards, J. R. (1995). Alternatives to Difference Scores as Dependent Variables in the Study of Congruence in Organizational Research. Organizational Behavior and Human Decision Processes, 64(3), 307–324.
#' @references Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence Testing for Psychological Research: A Tutorial. Advances in Methods and Practices in Psychological Science, 1(2), 259–269. https://doi.org/10.1177/2515245918770963
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(342356)
#' d <- data.frame(
#'   y1 = rnorm(50),
#'   y2 = rnorm(50),
#'   x = rnorm(50)
#' )
#' ddsc_sem(
#'   data = d, y1 = "y1", y2 = "y2",
#'   x = "x",
#'   q_sesoi = 0.20,
#'   min_cross_over_point_location = 1
#' )$results
#' }
ddsc_sem <- function(data,
                     x,
                     y1,
                     y2,
                     center_yvars = FALSE,
                     covariates = NULL,
                     estimator = "ML",
                     level = .95,
                     sampling.weights = NULL,
                     q_sesoi = 0,
                     min_cross_over_point_location = 0) {
  if (center_yvars) {
    pooled_mean <-
      mean(c(mean(data[, y1]), mean(data[, y2])))

    data[, y1] <- data[, y1] - pooled_mean
    data[, y2] <- data[, y2] - pooled_mean
  }

  # pooled sd for the components as the arithmetic mean
  pooled_sd_y1y2 <- sqrt(mean(c(
    stats::sd(data[, y1])^2,
    stats::sd(data[, y2])^2
  )))
  y1_scaled <- paste0(y1, "_scaled")
  y2_scaled <- paste0(y2, "_scaled")
  x_scaled <- paste0(x, "_scaled")

  data[, y1_scaled] <- data[, y1] / pooled_sd_y1y2
  data[, y2_scaled] <- data[, y2] / pooled_sd_y1y2
  data[, x_scaled] <- (data[, x] - mean(data[, x])) /
    stats::sd(data[, x])

  data[, "diff_score"] <- data[, y1] - data[, y2]
  data[, "diff_score_scaled"] <- data[, y1_scaled] - data[, y2_scaled]

  # descriptive statistics
  descriptives <-
    cbind(
      rbind(
        c(mean(data[, y1]), stats::sd(data[, y1])),
        c(mean(data[, y1_scaled]), stats::sd(data[, y1_scaled])),
        c(mean(data[, y2]), stats::sd(data[, y2])),
        c(mean(data[, y2_scaled]), stats::sd(data[, y2_scaled])),
        c(mean(data[, x]), stats::sd(data[, x])),
        c(mean(data[, x_scaled]), stats::sd(data[, x_scaled])),
        c(mean(data[, "diff_score"]), stats::sd(data[, "diff_score"])),
        c(mean(data[, "diff_score_scaled"]), stats::sd(data[, "diff_score_scaled"]))
      ),
      stats::cor(data[, c(
        y1, y1_scaled,
        y2, y2_scaled,
        x, x_scaled,
        "diff_score", "diff_score_scaled"
      )])
    )

  colnames(descriptives) <-
    c(
      "M", "SD", y1, y1_scaled, y2, y2_scaled, x, x_scaled,
      "diff_score", "diff_score_scaled"
    )

  # divergence of standardized slopes
  Cohens_q <- atanh(descriptives[x, y1]) - atanh(descriptives[x, y2])

  # divergence of common scale covariances with standardized X
  cov_scaled_x_y1 <- stats::cov(data[, x_scaled], data[, y1_scaled])
  cov_scaled_x_y2 <- stats::cov(data[, x_scaled], data[, y2_scaled])

  q_b <- atanh(cov_scaled_x_y1) - atanh(cov_scaled_x_y2)

  # Association table

  x_with_y1 <- c(
    cor = descriptives[x, y1],
    cov_scaled = cov_scaled_x_y1
  )

  x_with_y2 <- c(
    cor = descriptives[x, y2],
    cov_scaled = cov_scaled_x_y2
  )

  x_with_diff_score <- c(
    cor = descriptives[x, "diff_score"],
    q = Cohens_q,
    q_b = q_b
  )

  # structural equation model

  model <- paste0(
    paste0(y1_scaled, "~b_11*", x_scaled), "\n",
    paste0(y2_scaled, "~b_21*", x_scaled), "\n",
    paste0(y1_scaled, "~b_10*1"), "\n",
    paste0(y2_scaled, "~b_20*1"), "\n",
    paste0(y1_scaled, "~~res_cov_y1_y2*", y2_scaled), "\n",
    paste0(y1_scaled, "~~e_1*", y1_scaled), "\n",
    paste0(y2_scaled, "~~e_2*", y2_scaled), "\n",
    paste0(x_scaled, "~pred_int*1"), "\n",
    paste0(x_scaled, "~~pred_var*", x_scaled), "\n",
    paste0("diff_b11_b21:=b_11-b_21"), "\n",
    paste0("diff_b10_b20:=b_10-b_20"), "\n",
    paste0("q_b11_b21:=atanh(b_11)-atanh(b_21)"), "\n",
    paste0(
      "r_xy1:=b_11/",
      as.character(round(descriptives[y1_scaled, "SD"], 10))
    ), "\n",
    paste0(
      "r_xy2:=b_21/",
      as.character(round(descriptives[y2_scaled, "SD"], 10))
    ), "\n",
    paste0(
      "r_xy1_y2:=",
      paste0(
        "(r_xy1*",
        as.character(round(descriptives[y1, "SD"], 10))
      ), "-",
      paste0(
        "r_xy2*",
        as.character(round(descriptives[y2, "SD"], 10))
      ), ")/",
      as.character(round(descriptives["diff_score", "SD"], 10))
    ), "\n",
    paste0(
      "diff_rxy1_rxy2:=r_xy1-r_xy2"
    ), "\n",
    paste0("q_rxy1_rxy2:=atanh(r_xy1)-atanh(r_xy2)"), "\n",
    paste0("diff_y1_y2:=b_10-b_20"), "\n",
    paste0("cross_over_point:=(-1)*(b_10-b_20)/(b_11-b_21)"), "\n",
    paste0("sum_b11_b21:=b_11+b_21"), "\n",
    paste0("main_effect:=(b_11+b_21)/2"), "\n",
    paste0("interaction_vs_main_effect:=sqrt((b_11-b_21)^2)-sqrt(((b_11+b_21)/2)^2)"), "\n",
    paste0("diff_abs_b11_abs_b21:=sqrt(b_11^2)-sqrt(b_21^2)"), "\n",
    paste0("abs_diff_b11_b21:=sqrt((b_11-b_21)^2)"), "\n",
    paste0("abs_sum_b11_b21:=sqrt((b_11+b_21)^2)"), "\n",
    paste0("dadas_two_sided:=sqrt((b_11-b_21)^2)-sqrt((b_11+b_21)^2)")
  )

  # include covariates into the model

  if (!is.null(covariates)) {
    model <- paste0(model, "\n")
    model <- paste0(
      model,
      paste0(y1_scaled, "~", covariates, collapse = "\n")
    )
    model <- paste0(model, "\n")
    model <- paste0(
      model,
      paste0(y2_scaled, "~", covariates, collapse = "\n")
    )
    model <- paste0(model, "\n")
    model <- paste0(
      model,
      paste0(x_scaled, "~~", covariates, collapse = "\n")
    )
  }

  # include equivalence test for q and minimal-effects test for cross-over point

  if (!is.null(q_sesoi)) {
    model <- paste0(model, "\n")
    model <- paste0(
      model,
      paste0(
        "q_b_equivalence:=sqrt((q_b11_b21)^2)-",
        q_sesoi, "\n"
      ),
      paste0(
        "q_r_equivalence:=sqrt((q_rxy1_rxy2)^2)-",
        q_sesoi, "\n"
      ),
      paste0(
        "cross_over_point_equivalence:=sqrt((cross_over_point)^2)-",
        sqrt(min_cross_over_point_location^2)
      )
    )
  }

  if (!is.null(min_cross_over_point_location)) {
    model <- paste0(model, "\n")
    model <- paste0(
      model,
      paste0(
        "cross_over_point_minimal_effect:=sqrt((cross_over_point)^2)-",
        sqrt(min_cross_over_point_location^2)
      )
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
    data[, c(
      y1, y2, y1_scaled, y2_scaled,
      x, x_scaled, "diff_score", "diff_score_scaled"
    )]

  pars <- data.frame(lavaan::parameterestimates(fit,
    rsquare = TRUE,
    level = level
  ))

  # one-sided tests for absolute parameters

  labels <- c(
    "abs_diff_b11_b21", "abs_sum_b11_b21",
    "dadas_two_sided", "cross_over_point_minimal_effect"
  )

  osts <- pars[pars$label %in% labels, ]
  osts$p.pos <- stats::pnorm(osts[, "z"], lower.tail = F)

  # collect a summary of key results
  res.pars <- c(
    "r_xy1_y2", "r_xy1", "r_xy2",
    "b_11", "b_21", "b_10", "b_20", "res_cov_y1_y2",
    "diff_b10_b20",
    "diff_b11_b21", "diff_rxy1_rxy2",
    "q_b11_b21", "q_rxy1_rxy2",
    "cross_over_point", "sum_b11_b21",
    "main_effect", "interaction_vs_main_effect",
    "diff_abs_b11_abs_b21", "abs_diff_b11_b21", "abs_sum_b11_b21",
    "dadas_two_sided", "q_r_equivalence", "q_b_equivalence",
    "cross_over_point_equivalence", "cross_over_point_minimal_effect"
  )

  results <- pars[pars$label %in% res.pars, 4:ncol(pars)]
  rownames(results) <- results$label
  results <- results[, 2:ncol(results)]

  # retain the defined order of the parameters
  res_order <-
    sapply(res.pars, function(x) {
      which(rownames(results) == x)
    })

  results <- results[res_order, ]

  # replace two-sided tests with one-sided for absolute parameters

  results[labels, "pvalue"] <- osts[, "p.pos"]

  # replace two_sided_dadas name with dadas
  rownames(results)[rownames(results) == "dadas_two_sided"] <- "dadas"

  # add equivalence test to results with two one-sided tests
  results["q_r_equivalence", "pvalue"] <-
    stats::pnorm(results["q_r_equivalence", "z"], lower.tail = TRUE)

  results["q_r_equivalence", "ci.lower"] <- NA
  results["q_r_equivalence", "ci.upper"] <- NA

  results["q_b_equivalence", "pvalue"] <-
    stats::pnorm(results["q_b_equivalence", "z"], lower.tail = TRUE)

  results["q_b_equivalence", "ci.lower"] <- NA
  results["q_b_equivalence", "ci.upper"] <- NA

  results["cross_over_point_equivalence", "pvalue"] <-
    stats::pnorm(results["cross_over_point_equivalence", "z"], lower.tail = TRUE)

  results["cross_over_point_equivalence", "ci.lower"] <- NA
  results["cross_over_point_equivalence", "ci.upper"] <- NA

  results["cross_over_point_minimal_effect", "pvalue"] <-
    stats::pnorm(results["cross_over_point_minimal_effect", "z"], lower.tail = FALSE)

  results["cross_over_point_minimal_effect", "ci.lower"] <- NA
  results["cross_over_point_minimal_effect", "ci.upper"] <- NA

  # variance test in a separate model

  var_model <- paste0(
    paste0(y1_scaled, "~~cov_12*", y2_scaled), "\n",
    paste0(y1_scaled, "~~var_1*", y1_scaled), "\n",
    paste0(y2_scaled, "~~var_2*", y2_scaled), "\n",
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
    data = output.data,
    results = results
  )

  return(output)
}
