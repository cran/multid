#' Difference between two dependent Pearson's correlations (with common index)
#'
#' Calculates Cohen's q effect size statistic for difference between two correlations, r_yx1 and r_yx2.
#' Tests if Cohen's q is different from zero while accounting for dependency between the two correlations.
#'
#' @param data Data frame.
#' @param y Character. Variable name of the common index variable.
#' @param x1 Character. Variable name.
#' @param x2 Character. Variable name.
#' @param level Numeric. The confidence level required for the result output (Default .95)
#' @param missing Character. Treatment of missing values (e.g., "ML", default = listwise deletion)
#'
#' @return Parameter estimates from the fitted structural path model.
#' @export
#'
#' @examples
#' set.seed(3864)
#' d<-data.frame(y=rnorm(100),x=rnorm(100))
#' d$x1<-d$x+rnorm(100)
#' d$x2<-d$x+rnorm(100)
#' diff_two_dep_cors(data=d,y="y",x1="x1",x2="x2")
diff_two_dep_cors <-
  function(data, y, x1, x2, level = .95, missing = "default") {
    mod <- paste0(
      paste0(y, "~~var_y*", y), "\n",
      paste0(x1, "~~var_x1*", x1), "\n",
      paste0(x2, "~~var_x2*", x2), "\n",
      paste0(y, "~~cov_yx1*", x1), "\n",
      paste0(y, "~~cov_yx2*", x2), "\n",
      paste0(x1, "~~cov_x1x2*", x2), "\n",
      paste0(y, "~m_y*1"), "\n",
      paste0(x1, "~m_x1*1"), "\n",
      paste0(x2, "~m_x2*1"), "\n",
      paste0("r_yx1:=cov_yx1/(sqrt(var_y)*sqrt(var_x1))"), "\n",
      paste0("r_yx2:=cov_yx2/(sqrt(var_y)*sqrt(var_x2))"), "\n",
      paste0("r_x1x2:=cov_x1x2/(sqrt(var_x1)*sqrt(var_x2))"), "\n",
      paste0("diff_r_yx1_vs_r_yx2:=r_yx1-r_yx2"), "\n",
      paste0("z_yx1:=(1/2)*log((1+r_yx1)/(1-r_yx1))"), "\n",
      paste0("z_yx2:=(1/2)*log((1+r_yx2)/(1-r_yx2))"), "\n",
      paste0("Cohens_q:=z_yx1-z_yx2")
    )

    fit <- lavaan::sem(mod,
      data = data,
      estimator = "ML", missing = missing
    )

    n <- c(
      nobs = lavaan::lavInspect(fit, "nobs"),
      norig = lavaan::lavInspect(fit, "norig"),
      npar = lavaan::lavInspect(fit, "npar")
    )
    coverage <- lavaan::lavInspect(fit, "coverage")

    par.list <- unclass(lavaan::parameterestimates(fit,
      ci = TRUE,
      level = level,
      output = "data.frame"
    ))
    pars <- cbind(
      est = par.list$est,
      se = par.list$se,
      z = par.list$z,
      pvalue = par.list$pvalue,
      ci.lower = par.list$ci.lower,
      ci.upper = par.list$ci.upper
    )
    rownames(pars) <- par.list$label

    output <- list(estimates = pars, n = n, coverage = coverage)

    return(output)
  }
