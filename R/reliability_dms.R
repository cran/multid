#' Reliability calculation for difference score variable that is a difference between two mean variables calculated over upper-level units (e.g., sex differences across countries)
#'
#' Calculates reliability of difference score (Johns, 1981) based on two separate ICC2 values (Bliese, 2000), standard deviations of mean values over upper-level units, and correlations between the mean values across upper-level units.
#'
#' @param model Multilevel model fitted with lmer (default NULL)
#' @param data Long format data frame (default NULL)
#' @param var Character string. Name of the dependent variable or variable of which mean values are calculated.
#' @param diff_var Character string. A variable indicative of difference score components (two groups).
#' @param diff_var_values Vector. Values of the component score groups in diff_var.
#' @param group_var Character string. Upper-level clustering unit.
#'
#' @return A vector including ICC2s (r11 and r22), SDs (sd1, sd2, and sd_d12), means (m1, m2, and m_d12), correlation between means (r12), and reliability of the mean difference variable.
#' @references Bliese, P. D. (2000). Within-group agreement, non-independence, and reliability: Implications for data aggregation and analysis. In K. J. Klein & S. W. J. Kozlowski (Eds.), Multilevel theory, research, and methods in organizations: Foundations, extensions, and new directions (pp. 349–381). Jossey-Bass.
#' @references Johns, G. (1981). Difference score measures of organizational behavior variables: A critique. Organizational Behavior and Human Performance, 27(3), 443–463. https://doi.org/10.1016/0030-5073(81)90033-7
#' @export
#'
#' @examples set.seed(4317)
#' n2 <- 20
#' n1 <- 200
#' ri <- rnorm(n2, m = 0.5, sd = 0.2)
#' rs <- 0.5 * ri + rnorm(n2, m = 0.3, sd = 0.15)
#' d.list <- list()
#' for (i in 1:n2) {
#'   x <- rep(c(-0.5, 0.5), each = n1 / 2)
#'   y <- ri[i] + rs[i] * x + rnorm(n1)
#'   d.list[[i]] <- cbind(x, y, i)
#' }
#'
#' d <- data.frame(do.call(rbind, d.list))
#' names(d) <- c("x", "y", "cntry")
#' reliability_dms(
#'   data = d, diff_var = "x",
#'   diff_var_values = c(-0.5, 0.5), var = "y", group_var = "cntry"
#' )
reliability_dms <-
  function(model = NULL, data = NULL, diff_var, diff_var_values, var, group_var) {

    if (!is.null(data)){
      d.mod <- data
    } else (
      d.mod <- model@frame
    )

    d.mod[, group_var] <- as.factor(d.mod[, group_var])

    d.1 <- d.mod[d.mod[, diff_var] == diff_var_values[1], ]
    d.2 <- d.mod[d.mod[, diff_var] == diff_var_values[2], ]

    aov.1 <-
      stats::aov(stats::as.formula(paste0(var, "~", group_var)), d.1)
    r11.obs <- (summary(aov.1)[[1]][1, 3] -
      summary(aov.1)[[1]][2, 3]) /
      summary(aov.1)[[1]][1, 3]

    aov.2 <-
      stats::aov(stats::as.formula(paste0(var, "~", group_var)), d.2)
    r22.obs <- (summary(aov.2)[[1]][1, 3] -
      summary(aov.2)[[1]][2, 3]) /
      summary(aov.2)[[1]][1, 3]

    dm.1 <- stats::aggregate(d.1[, var],
                             mean, by = list(d.1[, group_var]))
    dm.2 <- stats::aggregate(d.2[, var],
                             mean, by = list(d.2[, group_var]))

    sd1.obs <- stats::sd(dm.1[, 2])
    sd2.obs <- stats::sd(dm.2[, 2])
    r12.obs <- stats::cor(dm.1[, 2], dm.2[, 2])

    sd_d12.obs <- stats::sd(dm.1[, 2]-dm.2[, 2])

    mean1.obs <- mean(dm.1[, 2])
    mean2.obs <- mean(dm.2[, 2])

    mean_d12.obs <- mean(dm.1[, 2]-dm.2[, 2])

    obs.ds.rel <-
      (r11.obs * sd1.obs^2 + r22.obs * sd2.obs^2 - 2 * r12.obs * sd1.obs * sd2.obs) /
        (sd1.obs^2 + sd2.obs^2 - 2 * r12.obs * sd1.obs * sd2.obs)

    output <-
      c(
        r11 = r11.obs, r22 = r22.obs,
        r12 = r12.obs, sd1 = sd1.obs,
        sd2 = sd2.obs, sd_d12 = sd_d12.obs,
        m1 = mean1.obs,
        m2 = mean2.obs,
        m_d12 = mean_d12.obs,
        reliability_dmsa = obs.ds.rel
      )
    return(output)
  }
