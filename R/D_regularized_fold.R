#' Use manually defined data folds for regularization and obtain estimates for each separately.
#'
#' @param data A data frame.
#' @param mv.vars Character vector. Variable names in the multivariate variable set.
#' @param group.var The name of the group variable.
#' @param group.values Vector of length 2, group values (e.g. c("male", "female) or c(0,1)).
#' @param alpha Alpha-value for penalizing function ranging from 0 to 1: 0 = ridge regression, 1 = lasso, 0.5 = elastic net (default).
#' @param s Which lambda value is used for predicted values? Either "lambda.min" (default) or "lambda.1se".
#' @param type.measure Which measure is used during cross-validation. Default "deviance".
#' @param rename.output Logical. Should the output values be renamed according to the group.values? Default TRUE.
#' @param fold.var Character string. Name of the fold variable.
#' @param append.data Logical. If TRUE, the original data is appended to the predicted variables.
#'
#' @return
#' \item{D}{Multivariate descriptive statistics and differences.}
#' \item{pred.dat}{A data.frame with predicted values.}
#' \item{cv.mod}{Regularized regression model from cv.glmnet.}
#' @seealso \code{\link[glmnet]{cv.glmnet}}
#' @export
#'
#' @examples set.seed(34246)
#' n1 <- 100
#' n2 <- 10
#' d <-
#'   data.frame(
#'     sex = sample(c("male", "female"), n1 * n2, replace = TRUE),
#'     fold = sample(x = LETTERS[1:n2], size = n1 * n2, replace = TRUE),
#'     x1 = rnorm(n1 * n2),
#'     x2 = rnorm(n1 * n2),
#'     x3 = rnorm(n1 * n2)
#'   )
#'
#' D_regularized_fold(
#'   data = d,
#'   mv.vars = c("x1", "x2", "x3"),
#'   group.var = "sex",
#'   group.values = c("female", "male"),
#'   fold.var = "fold"
#' )$D
D_regularized_fold <-
  function(data,
           mv.vars,
           group.var,
           group.values,
           alpha = 0.5,
           s = "lambda.min",
           type.measure = "deviance",
           rename.output = TRUE,
           fold.var,
           append.data = FALSE) {
    data$group.var.num <-
      ifelse(data[, group.var] == group.values[1], 1,
        ifelse(data[, group.var] == group.values[2], 0,
          NA
        )
      )

    fold.num.data <-
      cbind.data.frame(
        fold.num = c(1:length(unique(data[, fold.var]))),
        fold = unique(data[, fold.var])
      )

    # data frame joining information
    join_vars <- colnames(fold.num.data)[2]
    names(join_vars) <- fold.var

    data <- dplyr::left_join(
      x = data,
      y = fold.num.data,
      by = join_vars
    )

    foldid <- data[, "fold.num"] # $fold.num

    cv.mod <-
      glmnet::cv.glmnet(
        x = as.matrix(data[, c(mv.vars)]),
        y = data$group.var.num,
        family = c("binomial"),
        foldid = foldid,
        type.measure = type.measure,
        alpha = alpha
      )

    preds <- data.frame(
      group = data[, group.var],
      fold = data[, fold.var],
      pred = as.numeric(
        stats::predict(cv.mod,
          newx = as.matrix(data[, c(mv.vars)]),
          s = s
        )
      )
    )

    if (append.data) {
      preds <- cbind(preds, data)
    }

    D.folded <- list()
    fold.names <- unique(preds[, "fold"])

    for (i in fold.names) {
      D.folded[[i]] <- multid::d_pooled_sd(
        data = preds[preds$fold == i, ],
        var = "pred",
        group.var = "group",
        group.values = group.values,
        rename.output = rename.output
      )
    }

    D.folded.df <- do.call(rbind.data.frame, D.folded)

    D.folded.df <- cbind(
      D.folded.df,
      colwise_pool(
        data = D.folded.df,
        n1 = names(D.folded.df)[1],
        n2 = names(D.folded.df)[2],
        m1 = names(D.folded.df)[3],
        m2 = names(D.folded.df)[4],
        sd1 = names(D.folded.df)[5],
        sd2 = names(D.folded.df)[6]
      )
    )

    D.folded.df$d.sd.total <- D.folded.df$diff /
      D.folded.df$pooled.sd.total

    # rename pooled.sd columns if requested
    if (rename.output) {
      names(D.folded.df)[names(D.folded.df) == "pooled.sd.1"] <- paste0("pooled.sd.", group.values[1])
      names(D.folded.df)[names(D.folded.df) == "pooled.sd.2"] <- paste0("pooled.sd.", group.values[2])
    }

    D.folded.df <- D.folded.df[order(row.names(D.folded.df)), ]

    comb.output <- list(D = D.folded.df, preds = preds, cv.mod = cv.mod)

    return(comb.output)
  }
