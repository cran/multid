#' Predicting algebraic difference scores in multilevel model
#'
#' @param model Multilevel model fitted with lmerTest.
#' @param predictor Character string. Variable name of independent variable predicting difference score.
#' @param diff_var Character string. A variable indicative of difference score components (two groups).
#' @param diff_var_values Vector. Values of the component score groups in diff_var.
#'
#' @return A data frame including regression coefficients for component scores and dadas (One sided dadas-test for positivity of abs(b_11-b_21)-abs(b_11+b_21)).
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(95332)
#' n1 <- 10 # groups
#' n2 <- 10 # observations per group
#'
#' dat <- data.frame(
#'   group = rep(c(LETTERS[1:n1]), each = n2),
#'   x = sample(c(-0.5, 0.5), n1 * n2, replace = TRUE),
#'   w = rep(sample(1:5, n1, replace = TRUE), each = n2),
#'   y = sample(1:5, n1 * n2, replace = TRUE)
#' )
#' library(lmerTest)
#' fit <- lmerTest::lmer(y ~ x * w + (x | group),
#'   data = dat
#' )
#'
#' round(ml_dadas(fit,
#'   predictor = "w",
#'   diff_var = "x",
#'   diff_var_values = c(0.5, -0.5)
#' ), 3)
#' }
ml_dadas <- function(model,
                     predictor,
                     diff_var,
                     diff_var_values) {

  # reorder diff_var_values
  diff_var_values <-
    diff_var_values[order(diff_var_values)]

  # format list for contrast values
  at.list <- list(diff_var_values)
  names(at.list) <- diff_var

  trends <- emmeans::emtrends(
    object = model,
    specs = diff_var,
    var = predictor,
    lmerTest.limit = nrow(model@frame),
    disable.pbkrtest = TRUE,
    infer = c(FALSE, TRUE),
    at = at.list
  )
  # obtain trend signs for correct contrasts

  df.trends <- data.frame(trends)

  trend.diff <- df.trends[1, 2] - df.trends[2, 2]
  trend.sum <- df.trends[1, 2] + df.trends[2, 2]

  trend.strength <- df.trends[1, 2] - df.trends[2, 2]

  # define correct contrasts based on signs

  if (trend.diff > 0 & trend.sum < 0) {
    mlist <- list(
      abs_diff = c(1, -1),
      abs_sum = c(-1, -1)
    )
  } else if (trend.diff > 0 & trend.sum > 0) {
    mlist <- list(
      abs_diff = c(1, -1),
      abs_sum = c(1, 1)
    )
  } else if (trend.diff < 0 & trend.sum < 0) {
    mlist <- list(
      abs_diff = c(-1, 1),
      abs_sum = c(-1, -1)
    )
  } else if (trend.diff < 0 & trend.sum > 0) {
    mlist <- list(
      abs_diff = c(-1, 1),
      abs_sum = c(1, 1)
    )
  }

  temp.cont <-
    emmeans::contrast(trends,
      method = mlist,
      side = ">"
    )

  ml_abstest <-
    data.frame(emmeans::contrast(temp.cont,
      method = list(abs_test = c(1, -1)),
      side = ">"
    ))

  trends.df <- data.frame(trends)
  colnames(trends.df) <- colnames(data.frame(temp.cont))


  output <- rbind(
    trends.df,
    data.frame(temp.cont),
    ml_abstest
  )
  rownames(output) <- output$contrast
  output <- output[, 2:ncol(output)]

  rownames(output) <- c(
    rownames(output)[1:(nrow(output) - 1)],
    "dadas"
  )

  return(output)
}
