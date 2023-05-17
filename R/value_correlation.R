#' Testing and quantifying how much ipsatization (profile centering) influence associations between value and a correlate
#'
#' @param data A data frame.
#' @param rv Character string or vector. Variable name(s) of the non-ipsatized value variable(s) (raw value score).
#' @param cf Character string. Variable name of the common factor that is used for ipsatizing raw value scores.
#' @param correlate Character string. Name of the variable to which associations with values are examined.
#' @param estimator Character string. Estimator used in SEM (Default "ML").
#' @param scale_by_rv Logical. Is standard deviation of the raw non-ipsatized value score used for scaling the common factor as well? (Default FALSE)
#' @param standardize_correlate Logical. Should the correlate be standardized? (Default FALSE)
#' @param level Numeric. The confidence level required for the result output (Default .95)
#' @param sampling.weights Character string. Name of sampling weights variable.
#' @param sesoi Numeric. Smallest effect size of interest. Used for equivalence testing differences in ipsatized and non-ipsatized value associations (Default 0).
#'
#' @return
#' \item{parameter_estimates}{Parameter estimates from the structural equation model.}
#' \item{transformed_data}{Data frame with variables used in SEM (after scaling is applied).}
#' \item{results}{Summary of key results.}
#
#'
#' @export
#'
#' @examples
#' \dontrun{
#'set.seed(342356)
#'d <- data.frame(
#'  rv1 = rnorm(50),
#'  rv2 = rnorm(50),
#'  rv3 = rnorm(50),
#'  rv4 = rnorm(50),
#'  x = rnorm(50)
#')
#'d$cf<-rowMeans(d[,c("rv1","rv2","rv3","rv4")])
#'fit<-value_correlation(
#'  data = d, rv = c("rv1","rv2","rv3","rv4"), cf = "cf",
#'  correlate = "x",scale_by_rv = TRUE,
#'  standardize_correlate = TRUE,
#'  sesoi = 0.10
#')
#'round(fit$variability_summary,3)
#'round(fit$association_summary,3)
#' }

value_correlation <- function(data,
                              rv,
                              cf,
                              correlate,
                              scale_by_rv = FALSE,
                              standardize_correlate = FALSE,
                              estimator = "ML",
                              level = .95,
                              sampling.weights = NULL,
                              sesoi = 0) {
  if (scale_by_rv) {
    if (length(rv) == 1) {
      sd_rv <- stats::sd(data[, rv])
    } else {
      sd_rv <- mean(sapply(data[, rv], stats::sd))
    }

    data[, rv] <- data[, rv] / sd_rv
    data[, cf] <- data[, cf] / sd_rv
  }

  if (standardize_correlate) {
    data[, correlate] <-
      (data[, correlate] - mean(data[, correlate])) /
        stats::sd(data[, correlate])
  }

  # list where all results are saved
  result.list <- list()

  # repeat the same analysis for each value at a time
  for (i in 1:length(rv)) {
    # value name as temporary
    rv.temp <- rv[i]

    # descriptive observed statistics
    data[, "diff"] <- data[, rv.temp] - data[, cf]

    descriptives <-
      cbind(
        rbind(
          c(mean(data[, rv.temp]), stats::sd(data[, rv.temp])),
          c(mean(data[, cf]), stats::sd(data[, cf])),
          c(mean(data[, "diff"]), stats::sd(data[, "diff"])),
          c(mean(data[, correlate]), stats::sd(data[, correlate]))
        ),
        stats::cor(data[, c(rv.temp, cf, "diff", correlate)])
      )

    colnames(descriptives) <-
      c("M", "SD", rv.temp, cf, "diff", correlate)

    model <-
      paste0(
        paste0(rv.temp, "~~v_rv*", rv.temp), "\n",
        paste0(cf, "~~v_cf*", cf), "\n",
        paste0(correlate, "~~v_pr*", correlate), "\n",
        paste0(rv.temp, "~~cov_rv_pr*", correlate), "\n",
        paste0(cf, "~~cov_cf_pr*", correlate), "\n",
        paste0(rv.temp, "~~cov_rv_cf*", cf), "\n",
        paste0(rv.temp, "~mean_rv*1"), "\n",
        paste0(cf, "~mean_cf*1"), "\n",
        paste0(correlate, "~mean_pr*1"), "\n",
        paste0("mean_ip:=mean_rv-mean_cf"), "\n",
        paste0("sd_rv:=sqrt(v_rv)"), "\n",
        paste0("sd_cf:=sqrt(v_cf)"), "\n",
        paste0("sd_pr:=sqrt(v_pr)"), "\n",
        paste0("v_ip:=v_rv+v_cf-2*cov_rv_cf"), "\n",
        paste0("sd_ip:=sqrt(v_ip)"), "\n",
        paste0("var_diff_rv_cf:=v_rv-v_cf"), "\n",
        paste0("var_diff_rv_ip:=v_rv-v_ip"), "\n",
        paste0("var_diff_cf_ip:=v_cf-v_ip"), "\n",
        paste0("cor_rv_pr:=cov_rv_pr/(sd_rv*sd_pr)"), "\n",
        paste0("cor_cf_pr:=cov_cf_pr/(sd_cf*sd_pr)"), "\n",
        paste0("cor_rv_cf:=cov_rv_cf/(sd_rv*sd_cf)"), "\n",
        paste0("cor_ip_pr:=(cov_rv_pr-cov_cf_pr)/sd_ip"), "\n",
        paste0("cov_ip_pr:=(cov_rv_pr-cov_cf_pr)"), "\n",
        paste0("cor_diff_rv_pr_vs_ip_pr:=cor_rv_pr-cor_ip_pr"), "\n",
        paste0("cor_equi_low_rv_pr_vs_ip_pr:=cor_diff_rv_pr_vs_ip_pr+", sesoi), "\n",
        paste0("cor_equi_hig_rv_pr_vs_ip_pr:=cor_diff_rv_pr_vs_ip_pr-", sesoi), "\n",
        paste0("cov_diff_rv_pr_vs_ip_pr:=cov_rv_pr-cov_ip_pr"), "\n",
        paste0("cov_equi_low_rv_pr_vs_ip_pr:=cov_diff_rv_pr_vs_ip_pr+", sesoi), "\n",
        paste0("cov_equi_hig_rv_pr_vs_ip_pr:=cov_diff_rv_pr_vs_ip_pr-", sesoi), "\n"
      )

    # fit model

    fit <-
      lavaan::sem(
        model = model,
        data = data,
        estimator = estimator,
        sampling.weights = sampling.weights
      )

    # output.data <-
    #  data[, c(rv.temp, cf, "diff", correlate)]

    pars <- data.frame(lavaan::parameterestimates(fit,
      level = level
    ))

    # vector of parameters to be included in the results

    res.pars <- c(
      "mean_rv", "mean_ip", "mean_cf", "mean_pr",
      "sd_rv", "sd_ip", "sd_cf", "sd_pr",
      "v_rv", "v_ip", "v_cf", "v_pr",
      "cov_rv_pr", "cov_ip_pr", "cov_cf_pr", "cov_rv_cf",
      "cor_rv_pr", "cor_ip_pr", "cor_cf_pr", "cor_rv_cf",
      "var_diff_rv_cf",
      "var_diff_rv_ip",
      "cor_diff_rv_pr_vs_ip_pr",
      "cor_equi_low_rv_pr_vs_ip_pr",
      "cor_equi_hig_rv_pr_vs_ip_pr",
      "cov_diff_rv_pr_vs_ip_pr",
      "cov_equi_low_rv_pr_vs_ip_pr",
      "cov_equi_hig_rv_pr_vs_ip_pr"
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

    # add equivalence test to results with two one-sided tests
    results["cor_equi_low_rv_pr_vs_ip_pr", "pvalue"] <-
      stats::pnorm(results["cor_equi_low_rv_pr_vs_ip_pr", "z"], lower.tail = F)

    results["cor_equi_low_rv_pr_vs_ip_pr", "ci.lower"] <- NA
    results["cor_equi_low_rv_pr_vs_ip_pr", "ci.upper"] <- NA

    results["cor_equi_hig_rv_pr_vs_ip_pr", "pvalue"] <-
      stats::pnorm(results["cor_equi_hig_rv_pr_vs_ip_pr", "z"], lower.tail = T)

    results["cor_equi_hig_rv_pr_vs_ip_pr", "ci.lower"] <- NA
    results["cor_equi_hig_rv_pr_vs_ip_pr", "ci.upper"] <- NA

    results["cov_equi_low_rv_pr_vs_ip_pr", "pvalue"] <-
      stats::pnorm(results["cov_equi_low_rv_pr_vs_ip_pr", "z"], lower.tail = F)

    results["cov_equi_low_rv_pr_vs_ip_pr", "ci.lower"] <- NA
    results["cov_equi_low_rv_pr_vs_ip_pr", "ci.upper"] <- NA

    results["cov_equi_hig_rv_pr_vs_ip_pr", "pvalue"] <-
      stats::pnorm(results["cov_equi_hig_rv_pr_vs_ip_pr", "z"], lower.tail = T)

    results["cov_equi_hig_rv_pr_vs_ip_pr", "ci.lower"] <- NA
    results["cov_equi_hig_rv_pr_vs_ip_pr", "ci.upper"] <- NA

    result.list[[rv.temp]] <- list(
      descriptive_statistics = descriptives,
      parameter_estimates = pars,
      # transformed_data = output.data,
      results = results
    )
  }

  variability_summary <- list()
  association_summary <- list()

  for (j in 1:length(rv)) {
    variability_summary[[rv[j]]] <-
      result.list[[rv[j]]]$results[c(
        "sd_rv", "sd_ip", "sd_cf",
        "v_rv", "v_ip", "v_cf",
        "var_diff_rv_ip",
        "var_diff_rv_cf",
        "cor_rv_cf"
      ), ]

    association_summary[[rv[j]]] <-
      result.list[[rv[j]]]$results[c(
        "cor_rv_pr", "cor_ip_pr", "cor_cf_pr",
        "cov_rv_pr", "cov_ip_pr", "cov_cf_pr",
        "cor_diff_rv_pr_vs_ip_pr",
        "cor_equi_low_rv_pr_vs_ip_pr",
        "cor_equi_hig_rv_pr_vs_ip_pr",
        "cov_diff_rv_pr_vs_ip_pr",
        "cov_equi_low_rv_pr_vs_ip_pr",
        "cov_equi_hig_rv_pr_vs_ip_pr"
      ), ]
  }

  result.list[["variability_summary"]] <- do.call(rbind, variability_summary)
  result.list[["association_summary"]] <- do.call(rbind, association_summary)

  return(result.list)
}
