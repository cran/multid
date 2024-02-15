#' Predicting algebraic difference scores in multilevel model
#'
#' Decomposes difference score predictions to predictions of difference score components by probing simple effects at the levels of the binary moderator.
#'
#' @param model Multilevel model fitted with lmerTest.
#' @param predictor Character string. Variable name of independent variable predicting difference score.
#' @param diff_var Character string. A variable indicative of difference score components (two groups).
#' @param diff_var_values Vector. Values of the component score groups in diff_var.
#' @param scaled_estimates Logical. Are scaled estimates obtained? Does fit a reduced model for correct standard deviations. (Default FALSE)
#' @param re_cov_test Logical. Significance test for random effect covariation? Does fit a reduced model without the correlation. (Default FALSE)
#' @param var_boot_test Logical. Compare variance by lower-level groups at the upper-level in a reduced model with bootstrap? (Default FALSE)
#' @param nsim Numeric. Number of bootstrap simulations.
#' @param seed Numeric. Seed number for bootstrap simulations.
#' @param level Numeric. The confidence level required for the var_boot_test output (Default .95)
#' @param abs_diff_test Numeric. A value against which absolute difference between component score predictions is tested (Default 0).
#'
#' @return
#' \item{dadas}{A data frame including main effect, interaction, regression coefficients for component scores, dadas, and comparison between interaction and main effect.}
#' \item{scaled_estimates}{Scaled regression coefficients for difference score components and difference score.}
#' \item{vpc_at_reduced}{Variance partition coefficients in the model without the predictor and interactions.}
#' \item{re_cov_test}{Likelihood ratio significance test for random effect covariation.}
#' \item{boot_var_diffs}{List of different variance bootstrap tests.}
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
#'   w = sample(c(-0.5, 0.5), n1 * n2, replace = TRUE),
#'   x = rep(sample(1:5, n1, replace = TRUE), each = n2),
#'   y = sample(1:5, n1 * n2, replace = TRUE)
#' )
#' library(lmerTest)
#' fit <- lmerTest::lmer(y ~ x * w + (w | group),
#'   data = dat
#' )
#'
#' round(ml_dadas(fit,
#'   predictor = "x",
#'   diff_var = "w",
#'   diff_var_values = c(0.5, -0.5)
#' )$dadas, 3)
#' }
ml_dadas <- function(model,
                     predictor,
                     diff_var,
                     diff_var_values,
                     scaled_estimates = FALSE,
                     re_cov_test = FALSE,
                     var_boot_test = FALSE,
                     nsim = NULL,
                     level = .95,
                     seed = NULL,
                     abs_diff_test = 0) {
  .Deprecated("ddsc_ml")

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

  # one-way dadas test

  ml_abstest <-
    data.frame(emmeans::contrast(temp.cont,
      method = list(dadas = c(1, -1)),
      side = ">"
    ))

  # test of magnitude difference between interaction and main effect

  interaction_vs_main <-
    data.frame(emmeans::contrast(temp.cont,
      method = list(
        interaction_vs_main =
          c(1, -1 / 2)
      )
    ))

  # obtain main effect and interaction for the output
  model.coefs <- summary(model)$coefficients
  main_effect <- model.coefs[predictor, ]
  moderator_effect <- model.coefs[diff_var, ]
  interaction.term <-
    ifelse(paste0(predictor, ":", diff_var) %in% rownames(model.coefs),
      paste0(predictor, ":", diff_var),
      paste0(diff_var, ":", predictor)
    )
  interaction <- model.coefs[interaction.term, ]

  mi.coefs <-
    rbind(
      main_effect,
      moderator_effect,
      interaction
    )
  colnames(mi.coefs) <-
    c("estimate", "SE", "df", "t.ratio", "p.value")

  # simple slopes

  trends.df <- data.frame(trends)
  colnames(trends.df) <- colnames(data.frame(temp.cont))

  # absolute difference test

  adt <-
    emmeans::contrast(trends,
      method = mlist,
      side = ">",
      null = abs_diff_test
    )
  adt <- data.frame(adt)[1, ]
  adt$contrast <-
    paste0(adt$contrast, "_test_null", adt$null)
  adt <-
    adt[, -which(colnames(adt) == "null")]


  # combine to same output

  dadas <- rbind(
    trends.df,
    data.frame(temp.cont),
    ml_abstest,
    interaction_vs_main,
    adt
  )

  rownames(dadas) <- dadas$contrast
  dadas <- dadas[, 2:ncol(dadas)]

  dadas <- rbind(mi.coefs, dadas)

  output <- list(dadas = dadas)

  if (scaled_estimates | re_cov_test | var_boot_test) {
    # refit a reduced model to obtain scaled estimates
    # obtain fixed effects as character vector

    FEs <- attributes(stats::terms(model))$term.labels

    # obtain random effects and DV as a character vector

    DV <- as.character(stats::formula(model,
      random.only = T
    ))[2]
    RE <- as.character(stats::formula(model,
      random.only = T
    ))[3]

    # define the interaction term correctly

    if (paste0(predictor, ":", diff_var) %in% FEs) {
      int.term <- paste0(predictor, ":", diff_var)
    } else {
      int.term <- paste0(diff_var, ":", predictor)
    }

    # reformulate

    new.formula <-
      stats::reformulate(c(FEs[-c(
        which(FEs == predictor),
        which(FEs == int.term)
      )], RE), response = DV)

    # update to reduced model

    reduced_model <-
      stats::update(model, new.formula)

    # get component standard deviations

    vpc_at_reduced <-
      vpc_at(
        model = reduced_model,
        lvl1.var = diff_var,
        lvl1.values = diff_var_values
      )

    # get slope sd

    vc <-
      as.data.frame(lme4::VarCorr(reduced_model))

    slope_sd_reduced <-
      vc[vc$var1 == diff_var &
        !is.na(vc$var1) &
        is.na(vc$var2), "sdcor"]

    # compile to same frame

    scaled_estimates_df <-
      rbind(
        c(
          est = dadas[
            rownames(dadas) == diff_var_values[1],
            "estimate"
          ],
          scaling_SD =
            vpc_at_reduced[
              vpc_at_reduced[
                ,
                "lvl1.value"
              ] == diff_var_values[1],
              "Intercept.sd"
            ],
          scaled_est = dadas[
            rownames(dadas) == diff_var_values[1],
            "estimate"
          ] /
            vpc_at_reduced[
              vpc_at_reduced[
                ,
                "lvl1.value"
              ] == diff_var_values[1],
              "Intercept.sd"
            ]
        ),
        c(
          est = dadas[
            rownames(dadas) == diff_var_values[2],
            "estimate"
          ],
          scaling_SD =
            vpc_at_reduced[
              vpc_at_reduced[
                ,
                "lvl1.value"
              ] == diff_var_values[2],
              "Intercept.sd"
            ],
          scaled_est = dadas[
            rownames(dadas) == diff_var_values[2],
            "estimate"
          ] /
            vpc_at_reduced[
              vpc_at_reduced[
                ,
                "lvl1.value"
              ] == diff_var_values[2],
              "Intercept.sd"
            ]
        ),
        c(
          est = dadas[rownames(dadas) == diff_var_values[1], "estimate"] -
            dadas[rownames(dadas) == diff_var_values[2], "estimate"],
          scaling_SD = slope_sd_reduced,
          scaled_est = (dadas[rownames(dadas) == diff_var_values[1], "estimate"] -
            dadas[rownames(dadas) == diff_var_values[2], "estimate"]) /
            slope_sd_reduced
        )
      )

    coef_diff <-
      c(
        est = scaled_estimates_df[1, "est"] -
          scaled_estimates_df[2, "est"],
        scaling_SD = NA,
        scaled_est = scaled_estimates_df[1, "scaled_est"] -
          scaled_estimates_df[2, "scaled_est"]
      )

    component_correlation <-
      c(
        est = NA,
        scaling_SD = NA,
        scaled_est =
          (scaled_estimates_df[1, "scaling_SD"]^2 +
            scaled_estimates_df[2, "scaling_SD"]^2 -
            scaled_estimates_df[3, "scaling_SD"]^2) /
            (2 * scaled_estimates_df[1, "scaling_SD"] *
              scaled_estimates_df[2, "scaling_SD"])
      )

    rownames(scaled_estimates_df) <-
      c(diff_var_values, "difference")

    scaled_estimates_df <-
      rbind(
        scaled_estimates_df,
        coef_diff,
        component_correlation
      )
  }

  if (re_cov_test) {

    # drop the random effect correlation
    RE.no.cov <-
      gsub(
        x = RE, pattern = " | ",
        replacement = " || ",
        fixed = TRUE
      )

    # reformulate
    no.cov.formula <-
      stats::reformulate(c(FEs[-c(
        which(FEs == predictor),
        which(FEs == int.term)
      )], RE.no.cov), response = DV)

    # update to model without the covariance

    no.cov.mod <-
      stats::update(model, no.cov.formula)

    # test against the less reduced model

    re.cov.test <-
      stats::anova(
        no.cov.mod,
        reduced_model
      )

    # obtain important numbers

    re_cov_reduced <-
      vc[vc$var1 == "(Intercept)" &
        vc$var2 == diff_var &
        !is.na(vc$var1) &
        !is.na(vc$var2), c("vcov", "sdcor")]

    # compile to a frame

    re_cov_test_df <-
      c(
        RE_cov = re_cov_reduced[1, 1],
        RE_cor = re_cov_reduced[1, 2],
        Chisq = re.cov.test$Chisq[2],
        Df = re.cov.test$Df[2],
        p = re.cov.test$`Pr(>Chisq)`[2]
      )
  }

  if (var_boot_test) {

    # obtain covariance matrix in variance metric

    ranefmat <- function(.) {
      c(
        ranefcovmat = unlist(lme4::VarCorr(.)),
        sigma2 = stats::sigma(.)^2
      )
    }

    # parametric bootstrap for the covariance matrix

    boot.fit <-
      lme4::bootMer(reduced_model,
        FUN = ranefmat,
        nsim = nsim,
        seed = seed,
        use.u = FALSE,
        type = c("parametric"),
        verbose = FALSE
      )

    # obtain the variance estimates at lower-level values across bootstraps

    intvar1 <-
      boot.fit$t[, 1] +
      2 * boot.fit$t[, 2] * diff_var_values[1] +
      boot.fit$t[, 4] * diff_var_values[1]^2

    intvar2 <- boot.fit$t[, 1] +
      2 * boot.fit$t[, 2] * diff_var_values[2] +
      boot.fit$t[, 4] * diff_var_values[2]^2

    # obtain the values in the estimated model

    intvar1.t0 <-
      unname(boot.fit$t0[1] +
        2 * boot.fit$t0[2] * diff_var_values[1] +
        boot.fit$t0[4] * diff_var_values[1]^2)

    intvar2.t0 <- unname(boot.fit$t0[1] +
      2 * boot.fit$t0[2] * diff_var_values[2] +
      boot.fit$t0[4] * diff_var_values[2]^2)

    # estimate for the difference in variance and in variance ratio

    est.diff.var <- intvar1.t0 - intvar2.t0
    sd.diff.var.boot <- stats::sd(intvar1 - intvar2)

    est.ratio.var <- intvar1.t0 / intvar2.t0
    sd.ratio.var.boot <- stats::sd(intvar1 / intvar2)

    diff.var.norm <-
      c(
        est = est.diff.var,
        sd.boot = sd.diff.var.boot,
        p = 2 * (1 - stats::pnorm(abs(est.diff.var / sd.diff.var.boot))),
        LL = est.diff.var + stats::qnorm((1 - level) / 2) * sd.diff.var.boot,
        UL = est.diff.var + stats::qnorm(1 - (1 - level) / 2) * sd.diff.var.boot
      )

    ratio.var.norm <-
      c(
        est = est.ratio.var,
        sd.boot = sd.ratio.var.boot,
        p = 2 * (1 - stats::pnorm(abs((est.ratio.var - 1) / sd.ratio.var.boot))),
        LL = est.ratio.var + stats::qnorm((1 - level) / 2) * sd.ratio.var.boot,
        UL = est.ratio.var + stats::qnorm(1 - (1 - level) / 2) * sd.ratio.var.boot
      )

    # simple percentile interval

    diff.var.perc <-
      c(
        est = est.diff.var,
        LL = unname(stats::quantile(
          intvar1 - intvar2,
          stats::pnorm(stats::qnorm((1 - level) / 2))
        )),
        UL = unname(stats::quantile(
          intvar1 - intvar2,
          stats::pnorm(stats::qnorm(1 - (1 - level) / 2))
        ))
      )

    ratio.var.perc <-
      c(
        est = est.ratio.var,
        LL = unname(stats::quantile(
          intvar1 / intvar2,
          stats::pnorm(stats::qnorm((1 - level) / 2))
        )),
        UL = unname(stats::quantile(
          intvar1 / intvar2,
          stats::pnorm(stats::qnorm(1 - (1 - level) / 2))
        ))
      )

    # bias corrected

    bias.diff <-
      sum((intvar1 - intvar2) > est.diff.var) / length(intvar1)

    diff.LQ <- stats::qnorm(.025) - 2 * stats::qnorm(bias.diff)
    diff.UQ <- stats::qnorm(.975) - 2 * stats::qnorm(bias.diff)

    diff.var.bias <-
      c(
        est = est.diff.var,
        LL = unname(stats::quantile(
          intvar1 - intvar2,
          stats::pnorm(diff.LQ)
        )),
        UL = unname(stats::quantile(
          intvar1 - intvar2,
          stats::pnorm(diff.UQ)
        ))
      )

    bias.ratio <-
      sum((intvar1 / intvar2) > est.ratio.var) / length(intvar1)

    ratio.LQ <- stats::qnorm(.025) - 2 * stats::qnorm(bias.ratio)
    ratio.UQ <- stats::qnorm(.975) - 2 * stats::qnorm(bias.ratio)

    ratio.var.bias <-
      c(
        est = est.ratio.var,
        LL = unname(stats::quantile(
          intvar1 / intvar2,
          stats::pnorm(ratio.LQ)
        )),
        UL = unname(stats::quantile(
          intvar1 / intvar2,
          stats::pnorm(ratio.UQ)
        ))
      )

    boot_var_diffs <-
      list(
        norm_boot_var_diff = diff.var.norm,
        perc_boot_var_diff = diff.var.perc,
        bias_boot_var_diff = diff.var.bias,
        norm_boot_var_ratio = ratio.var.norm,
        perc_boot_var_ratio = ratio.var.perc,
        bias_boot_var_ratio = ratio.var.bias
      )
  }

  # Compile output based on arguments

  if (scaled_estimates | re_cov_test | var_boot_test) {
    output <- c(output,
      vpc_at_reduced = list(vpc_at_reduced)
    )
  }

  if (scaled_estimates) {
    output <- c(output,
      scaled_estimates = list(scaled_estimates_df)
    )
  }

  if (var_boot_test) {
    output <- c(output,
      boot_var_diffs = list(boot_var_diffs)
    )
  }

  if (re_cov_test) {
    output <- c(output,
      re_cov_test = list(re_cov_test_df)
    )
  }

  return(output)
}
