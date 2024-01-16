#' Deconstructing difference score correlation with multi-level modeling
#'
#' Can be used for either pre-fitted lmer-models or to long format data.
#'
#' @param model Multilevel model fitted with lmerTest.
#' @param data Data frame.
#' @param predictor Character string. Variable name of independent variable predicting difference score (i.e., x).
#' @param moderator Character string. Variable name indicative of difference score components (w).
#' @param moderator_values Vector. Values of the component score groups in moderator (i.e., y1 and y2).
#' @param DV Character string. Name of the dependent variable (if model is not supplied as input).
#' @param lvl2_unit Character string. Name of the level-2 clustering variable (if model is not supplied as input).
#' @param covariates Character string or vector. Variable names of covariates (Default NULL).
#' @param scaling_sd Character string (either default "observed" or "model"). Are the simple slopes scaled with observed or model-based SDs?
#' @param re_cov_test Logical. Significance test for random effect covariation? (Default FALSE)
#' @param var_boot_test Logical. Compare variance by lower-level groups at the upper-level in a reduced model with bootstrap? (Default FALSE)
#' @param nsim Numeric. Number of bootstrap simulations.
#' @param seed Numeric. Seed number for bootstrap simulations.
#' @param level Numeric. The confidence level required for the var_boot_test output (Default .95)

#'
#' @return
#' \item{results}{Summary of key results.}
#' \item{descriptives}{Means, standard deviations, and intercorrelations at level 2.}
#' \item{vpc_at_moderator_values}{Variance partition coefficients for moderator values in the model without the predictor and interactions.}
#' \item{model}{Fitted lmer object.}
#' \item{reduced_model}{Fitted lmer object without the predictor.}
#' \item{lvl2_data}{Data summarized at level 2.}
#' \item{ddsc_sem_fit}{ddsc_sem object fitted to level 2 data.}
#' \item{re_cov_test}{Likelihood ratio significance test for random effect covariation.}
#' \item{boot_var_diffs}{List of different variance bootstrap tests.}
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(95332)
#' n1 <- 10 # groups
#' n2 <- 10 # observations per group
#' dat <- data.frame(
#'   group = rep(c(LETTERS[1:n1]), each = n2),
#'   w = sample(c(-0.5, 0.5), n1 * n2, replace = TRUE),
#'   x = rep(sample(1:5, n1, replace = TRUE), each = n2),
#'   y = sample(1:5, n1 * n2, replace = TRUE)
#' )
#' library(lmerTest)
#' fit <- lmerTest::lmer(y ~ x * w + (w | group),
#'                       data = dat
#' )
#' round(ddsc_ml(model=fit,
#'               predictor="x",
#'               moderator="w",
#'               moderator_values=c(0.5,-0.5))$results,3)
#'
#' round(ddsc_ml(data=dat,
#'               DV="y",
#'               lvl2_unit="group",
#'               predictor="x",
#'               moderator="w",
#'               moderator_values=c(0.5,-0.5))$results,3)
#'
#' }
ddsc_ml <- function(model = NULL,
                    data = NULL,
                    predictor,
                    moderator,
                    moderator_values,
                    DV = NULL,
                    lvl2_unit = NULL,
                    re_cov_test = FALSE,
                    var_boot_test = FALSE,
                    nsim = NULL,
                    level = .95,
                    seed = NULL,
                    covariates = NULL,
                    scaling_sd = "observed") {

  # reorder moderator_values
  # moderator_values <-
  #  moderator_values[order(moderator_values)]

  # if data not supplied, obtain it from the lmer -object
  # obtain also lvl2_unit and model dependent variable
  if (!is.null(model)) {
    data <- data.frame(model@frame)

    lvl2_unit <- attr(stats::terms(stats::formula(model)), "term.labels")[
      grepl("\\|", attr(stats::terms(stats::formula(model)), "term.labels"))
    ]

    lvl2_unit <- sub(".*\\|\\s*", "", lvl2_unit)
    lvl2_unit <- gsub("\\s+", "", lvl2_unit)

    DV <- as.character(stats::formula(model)[[2]])
  }
  # reconstruct the data so that y1 is 0.5 and y2 -0.5
  # this helps in getting correct signs for all estimates
  # as well as effect coding properly
  data[, moderator] <-
    ifelse(data[, moderator] == moderator_values[1], 0.5,
           ifelse(data[, moderator] == moderator_values[2], -0.5, NA)
    )

  # calculate lvl2 means and save to data frame

  lvl2_means_y1 <- tapply(data[data[, moderator] == moderator_values[1], DV],
                          data[data[, moderator] == moderator_values[1], lvl2_unit],
                          mean,
                          na.rm = TRUE
  )

  lvl2_means_y2 <- tapply(data[data[, moderator] == moderator_values[2], DV],
                          data[data[, moderator] == moderator_values[2], lvl2_unit],
                          mean,
                          na.rm = TRUE
  )

  lvl2_means_x <- tapply(data[data[, moderator] == moderator_values[2], predictor],
                         data[data[, moderator] == moderator_values[2], lvl2_unit],
                         mean,
                         na.rm = TRUE
  )

  lvl2_data <- data.frame(
    y1 = lvl2_means_y1,
    y2 = lvl2_means_y2,
    x = lvl2_means_x
  )

  names(lvl2_data) <- c(
    "means_y1",
    "means_y2",
    predictor
  )

  # obtain descriptives
  ddsc_sem_fit <-
    multid::ddsc_sem(
      data = lvl2_data[stats::complete.cases(lvl2_data), ],
      x = predictor,
      y1 = "means_y1",
      y2 = "means_y2"
    )
  descriptives <- ddsc_sem_fit$descriptives
  sem_variance_test <- ddsc_sem_fit$variance_test

  # check if the predictor is properly scaled if model was the input
  if (abs(descriptives[predictor, "SD"] - 1) > .005 &
      !is.null(model)) {
    warning(paste0(
      "Predictor not properly standardized, SD = ",
      descriptives[predictor, "SD"]
    ))
  }

  # construct and run a model if not provided as input
  if (is.null(model) & !is.null(DV)) {
    # standardize the predictor with country-level mean and sd
    data[, predictor] <- (data[, predictor] - descriptives[predictor, "M"]) /
      descriptives[predictor, "SD"]

    model_formula <-
      stats::as.formula(paste0(
        DV, "~",
        moderator, "*", predictor, "+",
        paste0(covariates, collapse = "+"),
        ifelse(is.null(covariates), "(", "+("),
        moderator, "|", lvl2_unit, ")"
      ))

    model <- lmerTest::lmer(
      formula = model_formula, data = data,
      control = lme4::lmerControl(optimizer = "bobyqa")
    )
  }

  # fit a reduced model without predictor and cross-level interaction

  # obtain fixed effects as character vector
  FEs <- attributes(stats::terms(model))$term.labels

  # obtain random effects and DV as a character vector
  DV <- as.character(stats::formula(model,
                                    random.only = TRUE
  ))[2]
  RE <- as.character(stats::formula(model,
                                    random.only = TRUE
  ))[3]

  # define the interaction term correctly

  if (paste0(predictor, ":", moderator) %in% FEs) {
    int.term <- paste0(predictor, ":", moderator)
  } else {
    int.term <- paste0(moderator, ":", predictor)
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

  vpc_at_reduced <-
    vpc_at(
      model = reduced_model,
      lvl1.var = moderator,
      lvl1.values = moderator_values
    )

  # get scaling SDs from the reduced model if requested
  # if not, use observed (default)

  if (scaling_sd == "model") {
    # get sds from variance partition

    slope_sd_reduced <-
      vc[vc$var1 == moderator &
           !is.na(vc$var1) &
           is.na(vc$var2), "sdcor"]

    # compile scaling SDs

    scaling_SDs <- c(
      vpc_at_reduced$Intercept.sd,
      mean(vpc_at_reduced$Intercept.sd),
      slope_sd_reduced
    )
  } else {
    sd_y1 <- descriptives["means_y1", "SD"]
    sd_y2 <- descriptives["means_y2", "SD"]
    r_y1y2 <- descriptives["means_y1", "means_y2"]
    sd_diff <- descriptives["diff_score", "SD"]
    sd_pooled <- (sd_y1 + sd_y2) / 2

    scaling_SDs <- c(sd_y1, sd_y2, sd_pooled, sd_diff)
  }

  names(scaling_SDs) <- c("SD_y1", "SD_y2", "SD_pooled", "SD_diff_score")

  # format list for contrast values
  at.list <- list(moderator_values)
  names(at.list) <- moderator

  # obtain unstandardized slopes
  trends <- emmeans::emtrends(
    object = model,
    specs = moderator,
    var = predictor,
    lmerTest.limit = nrow(model@frame),
    disable.pbkrtest = TRUE,
    infer = c(FALSE, TRUE),
    at = at.list
  )

  # obtain slopes with pooled sd
  trends_bscale <- emmeans::emtrends(
    object = model,
    specs = moderator,
    var = paste0("scale(", predictor, ",1,1/", scaling_SDs[3], ")"),
    lmerTest.limit = nrow(model@frame),
    disable.pbkrtest = TRUE,
    infer = c(FALSE, TRUE),
    at = at.list
  )

  # obtain separately standardized slopes

  trend_1_rscale <- emmeans::emtrends(
    object = model,
    specs = moderator,
    var = paste0("scale(", predictor, ",1,1/", scaling_SDs[1], ")"),
    lmerTest.limit = nrow(model@frame),
    disable.pbkrtest = TRUE,
    infer = c(FALSE, TRUE),
    at = at.list
  )

  trend_2_rscale <- emmeans::emtrends(
    object = model,
    specs = moderator,
    var = paste0("scale(", predictor, ",1,1/", scaling_SDs[2], ")"),
    lmerTest.limit = nrow(model@frame),
    disable.pbkrtest = TRUE,
    infer = c(FALSE, TRUE),
    at = at.list
  )

  trends_rscale <-
    rbind(trend_1_rscale[1],
          trend_2_rscale[2],
          adjust = "none"
    )

  # obtain trends scaled with difference score SD

  trends_diffscale <- emmeans::emtrends(
    object = model,
    specs = moderator,
    var = paste0("scale(", predictor, ",1,1/", scaling_SDs[4], ")"),
    lmerTest.limit = nrow(model@frame),
    disable.pbkrtest = TRUE,
    infer = c(FALSE, TRUE),
    at = at.list
  )

  # obtain trend signs for correct DADAS contrasts

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

  # obtain contrasts for unstandardized slopes
  temp.cont <-
    emmeans::contrast(trends,
                      method = mlist,
                      side = ">"
    )

  temp.cont.df <- data.frame(temp.cont)
  rownames(temp.cont.df) <- c("abs_diff", "abs_sum")
  # obtain contrasts for pooled sd standardized slopes
  temp.cont_bscale <-
    emmeans::contrast(trends_bscale,
                      method = mlist,
                      side = ">"
    )
  temp.cont_bscale.df <- data.frame(temp.cont_bscale)
  rownames(temp.cont_bscale.df) <- c("abs_diff_bscale", "abs_sum_bscale")
  # obtain contrasts for separately standardized slopes
  temp.cont_rscale <-
    emmeans::contrast(trends_rscale,
                      method = mlist,
                      side = ">"
    )
  temp.cont_rscale.df <- data.frame(temp.cont_rscale)
  rownames(temp.cont_rscale.df) <- c("abs_diff_rscale", "abs_sum_rscale")

  # one-sided dadas test

  ml_abstest <-
    data.frame(emmeans::contrast(temp.cont,
                                 method = list(dadas = c(1, -1)),
                                 side = ">"
    ))
  rownames(ml_abstest) <- "dadas"
  ml_abstest_bscale <-
    data.frame(emmeans::contrast(temp.cont_bscale,
                                 method = list(dadas = c(1, -1)),
                                 side = ">"
    ))
  rownames(ml_abstest_bscale) <- "dadas_bscale"
  ml_abstest_rscale <-
    data.frame(emmeans::contrast(temp.cont_rscale,
                                 method = list(dadas = c(1, -1)),
                                 side = ">"
    ))
  rownames(ml_abstest_rscale) <- "dadas_rscale"
  # test of magnitude difference between interaction and main effect

  interaction_vs_main <-
    data.frame(emmeans::contrast(temp.cont,
                                 method = list(
                                   interaction_vs_main =
                                     c(1, -1 / 2)
                                 )
    ))
  rownames(interaction_vs_main) <- "interaction_vs_main"

  interaction_vs_main_bscale <-
    data.frame(emmeans::contrast(temp.cont_bscale,
                                 method = list(
                                   interaction_vs_main =
                                     c(1, -1 / 2)
                                 )
    ))
  rownames(interaction_vs_main_bscale) <- "interaction_vs_main_bscale"

  interaction_vs_main_rscale <-
    data.frame(emmeans::contrast(temp.cont_rscale,
                                 method = list(
                                   interaction_vs_main =
                                     c(1, -1 / 2)
                                 )
    ))
  rownames(interaction_vs_main_rscale) <- "interaction_vs_main_rscale"

  # obtain main effect and interaction for the output
  model.coefs <- summary(model)$coefficients
  main_effect <- model.coefs[predictor, ]
  moderator_effect <- model.coefs[moderator, ]
  interaction.term <-
    ifelse(paste0(predictor, ":", moderator) %in% rownames(model.coefs),
           paste0(predictor, ":", moderator),
           paste0(moderator, ":", predictor)
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
  rownames(trends.df) <- c("w_11", "w_21")

  trends_bscale.df <- data.frame(trends_bscale)
  colnames(trends_bscale.df) <- colnames(data.frame(temp.cont))
  rownames(trends_bscale.df) <- c("b_11", "b_21")

  trends_rscale.df <- data.frame(trends_rscale)
  colnames(trends_rscale.df) <- colnames(data.frame(temp.cont))
  rownames(trends_rscale.df) <- c("r_xy1", "r_xy2")

  # cross-over point

  # moderator effect exactly at predictor midpoint
  at.list.mod <- list()
  at.list.mod[[predictor]] <- 0
  at.list.mod[[moderator]] <- moderator_values

  mod_eff <-
    emmeans::emmeans(
      object = model,
      specs = moderator,
      # var = predictor,
      lmerTest.limit = nrow(model@frame),
      disable.pbkrtest = TRUE,
      infer = c(FALSE, TRUE),
      at = at.list.mod
    )

  # compile cross-over point effect
  cop_effs <-
    data.frame(rbind(
      emmeans::contrast(mod_eff,
                        adjust = "none",
                        method = list(mod_eff = c(1, -1))
      ),
      emmeans::contrast(trends,
                        adjust = "none",
                        method = list(int = c(1, -1))
      ),
      adjust = "none"
    ))

  cross_over_point <-
    (-1) * cop_effs[1, "estimate"] /
    cop_effs[2, "estimate"]

  # absolute difference test

  # adt <-
  #  emmeans::contrast(trends,
  #                    method = mlist,
  #                    side = ">",
  #                    null = abs_diff_test
  #  )
  # adt <- data.frame(adt)[1, ]
  # adt$contrast <-
  #  paste0(adt$contrast, "_test_null", adt$null)
  # adt <-
  #  adt[, -which(colnames(adt) == "null")]

  # Estimates for slope non-parallelism

  # obtain difference score correlation
  r_xy1y2 <- data.frame(emmeans::contrast(trends_diffscale,
                                          adjust = "none",
                                          method = list(r_xy1y2 = c(1, -1))
  ))
  rownames(r_xy1y2) <- c("r_xy1y2")

  # obtain Cohen's q
  q_rxy1_rxy2 <-
    atanh(trends_rscale.df[trends_rscale.df$contrast == moderator_values[1], "estimate"]) -
    atanh(trends_rscale.df[trends_rscale.df$contrast == moderator_values[2], "estimate"])

  # obtain Cohen's q for pooled sd standardized slopes
  q_b11_b21 <-
    atanh(trends_bscale.df[trends_bscale.df$contrast == moderator_values[1], "estimate"]) -
    atanh(trends_bscale.df[trends_bscale.df$contrast == moderator_values[2], "estimate"])

  # combine to same result output

  results <-
    rbind(
      r_xy1y2,
      trends.df,
      trends_rscale.df,
      trends_bscale.df
    )

  results <- results[, -which(names(results) == "contrast")]
  results <- rbind(
    results,
    mi.coefs
  )

  # combine singe number outputs
  sn_results <-
    data.frame(
      estimate =
        c(q_b11_b21, q_rxy1_rxy2, cross_over_point),
      SE = NA,
      df = NA,
      t.ratio = NA,
      p.value = NA
    )

  rownames(sn_results) <-
    c("q_b11_b21", "q_rxy1_rxy2", "cross_over_point")

  results <- rbind(
    results,
    sn_results,
    rbind(
      interaction_vs_main,
      interaction_vs_main_bscale,
      interaction_vs_main_rscale
    )[, 2:6],
    rbind(
      ml_abstest,
      ml_abstest_bscale,
      ml_abstest_rscale
    )[, 2:6],
    rbind(
      temp.cont.df,
      temp.cont_bscale.df,
      temp.cont_rscale.df
    )[, 2:6]
  )

  scaling_SDs["VR"] <- unname(scaling_SDs["SD_y1"]^2)/
    (scaling_SDs["SD_y2"]^2)

  output <- list(
    results = results,
    descriptives = descriptives,
    vpc_at_moderator_values = vpc_at_reduced,
    model = model,
    reduced_model = reduced_model,
    lvl2_data = lvl2_data,
    ddsc_sem_fit = ddsc_sem_fit,
    SDs = scaling_SDs
  )

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
    vc <-
      as.data.frame(lme4::VarCorr(reduced_model))

    re_cov_reduced <-
      vc[vc$var1 == "(Intercept)" &
           vc$var2 == moderator &
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
      2 * boot.fit$t[, 2] * moderator_values[1] +
      boot.fit$t[, 4] * moderator_values[1]^2

    intvar2 <- boot.fit$t[, 1] +
      2 * boot.fit$t[, 2] * moderator_values[2] +
      boot.fit$t[, 4] * moderator_values[2]^2

    # obtain the values in the estimated model

    intvar1.t0 <-
      unname(boot.fit$t0[1] +
               2 * boot.fit$t0[2] * moderator_values[1] +
               boot.fit$t0[4] * moderator_values[1]^2)

    intvar2.t0 <- unname(boot.fit$t0[1] +
                           2 * boot.fit$t0[2] * moderator_values[2] +
                           boot.fit$t0[4] * moderator_values[2]^2)

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
