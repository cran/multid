#' Standardized mean difference with pooled standard deviation
#'
#' @param data A data frame.
#' @param var A continuous variable for which difference is estimated.
#' @param group.var The name of the group variable.
#' @param group.values Vector of length 2, group values (e.g. c("male", "female) or c(0,1)).
#' @param rename.output Logical. Should the output values be renamed according to the group.values? Default TRUE.
#' @param infer Logical. Statistical inference with Welch test? (default FALSE)
#'
#' @return Descriptive statistics and mean differences
#' @export
#'
#' @examples
#' d_pooled_sd(iris[iris$Species == "setosa" | iris$Species == "versicolor", ],
#'   var = "Petal.Length", group.var = "Species",
#'   group.values = c("setosa", "versicolor"), infer = TRUE
#' )
d_pooled_sd <-
  function(data,
           var,
           group.var,
           group.values,
           rename.output = TRUE,
           infer = FALSE) {
    dat1 <- data[data[, group.var] == group.values[1], ]
    dat2 <- data[data[, group.var] == group.values[2], ]

    sd.1 <- stats::sd(dat1[, var])
    n.1 <- length(dat1[, var])

    sd.2 <- stats::sd(dat2[, var])
    n.2 <- length(dat2[, var])

    pooled.sd <- sqrt(((n.1 - 1) * sd.1^2 + (n.2 - 1) * sd.2^2) / (n.1 + n.2 - 2))
    diff <- mean(dat1[, var]) - mean(dat2[, var])
    m.1 <- mean(dat1[, var])
    m.2 <- mean(dat2[, var])
    D <- diff / pooled.sd
    output <- cbind(n.1, n.2, m.1, m.2, sd.1, sd.2, pooled.sd, diff, D)

    if (rename.output) {
      colnames(output) <-
        c(
          paste0("n.", group.values[1]),
          paste0("n.", group.values[2]),
          paste0("m.", group.values[1]),
          paste0("m.", group.values[2]),
          paste0("sd.", group.values[1]),
          paste0("sd.", group.values[2]),
          "pooled.sd",
          "diff",
          "D"
        )
    }

    if (infer) {
      t_test <-
        stats::t.test(dat1[, var], dat2[, var], var.equal = FALSE)

      t_stat <- unname(t_test$statistic)
      df <- unname(t_test$parameter)
      p <- unname(t_test$p.value)

      output <- cbind(output, t_stat, df, p)
    }


    return(output)
  }
