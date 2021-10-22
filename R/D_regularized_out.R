#' Use separate data partition for regularization and estimation.
#'
#' @param data A data frame.
#' @param mv.vars Character vector. Variable names in the multivariate variable set.
#' @param group.var The name of the group variable.
#' @param group.values Vector of length 2, group values (e.g. c("male", "female) or c(0,1)).
#' @param alpha Alpha-value for penalizing function ranging from 0 to 1: 0 = ridge regression, 1 = lasso, 0.5 = elastic net (default).
#' @param nfolds Number of folds used for obtaining lambda (range from 3 to n-1, default 10).
#' @param s Which lambda value is used for predicted values? Either "lambda.min" (default) or "lambda.1se".
#' @param type.measure Which measure is used during cross-validation. Default "deviance".
#' @param rename.output Logical. Should the output values be renamed according to the group.values? Default TRUE.
#' @param size Integer. Size of regularization data per each group. Default 1/4 of cases.
#' @param pcc Logical. Include probabilities of correct classification? Default FALSE.
#' @param auc Logical. Include area under the receiver operating characteristics? Default FALSE.
#' @param pred.prob Logical. Include table of predicted probabilities? Default FALSE.
#' @param prob.cutoffs Vector. Cutoffs for table of predicted probabilities. Default seq(0,1,0.20).
#'
#' @return
#' \item{D}{Multivariate descriptive statistics and differences.}
#' \item{pred.dat}{A data.frame with predicted values.}
#' \item{cv.mod}{Regularized regression model from cv.glmnet.}
#' \item{P.table}{Table of predicted probabilities by cutoffs.}
#' @export
#'
#' @examples D_regularized_out(
#'   data = iris[iris$Species == "setosa" |
#'     iris$Species == "versicolor", ],
#'   mv.vars = c(
#'     "Sepal.Length", "Sepal.Width",
#'     "Petal.Length", "Petal.Width"
#'   ),
#'   group.var = "Species",
#'   group.values = c("setosa", "versicolor"),
#'   size = 40,
#'   pcc = TRUE
#' )$D
D_regularized_out <-
  function(data,
           mv.vars,
           group.var,
           group.values,
           alpha = 0.5,
           nfolds = 10,
           s = "lambda.min",
           type.measure = "deviance",
           rename.output = TRUE,
           size = NULL,
           pcc = FALSE,
           auc = FALSE,
           pred.prob = FALSE,
           prob.cutoffs = seq(from = 0, to = 1, by = 0.20)) {
    data$group.var.num <-
      ifelse(data[, group.var] == group.values[1], 1,
        ifelse(data[, group.var] == group.values[2], 0,
          NA
        )
      )

    if (is.null(size)) {
      size <- round(nrow(data) / 4, 0)
    } else {
      size <- size
    }

    data$row.nmbr <- rownames(data)

    data.grouped <- dplyr::group_by(data, group.var.num)

    train.data <- dplyr::sample_n(data.grouped,
      size = size,
      replace = F
    )

    test.data <- data[!(data$row.nmbr %in% train.data$row.nmbr), ]
    train.data <- dplyr::ungroup(train.data)

    cv.mod <-
      glmnet::cv.glmnet(
        x = as.matrix(train.data[, c(mv.vars)]),
        y = train.data$group.var.num,
        family = c("binomial"),
        nfolds = nfolds,
        type.measure = type.measure,
        alpha = alpha
      )

    preds <- data.frame(
      group = test.data[, group.var],
      pred = as.numeric(
        stats::predict(cv.mod,
          newx = as.matrix(test.data[, c(mv.vars)]),
          s = s
        )
      )
    )

    D <- multid::d_pooled_sd(
      data = preds,
      var = "pred",
      group.var = "group",
      group.values = group.values,
      rename.output = rename.output
    )

    # Add pcc

    if (pcc) {
      pcc.out <-
        pcc(
          data = preds,
          pred.var = "pred",
          group.var = "group",
          group.values = group.values
        )

      # inherit naming from D frame
      colnames(pcc.out) <-
        c(
          paste0("pcc.", substr(colnames(D)[1:2], 3, stop = 999)),
          "pcc.total"
        )

      D <- cbind(D, pcc.out)
    }

    if (auc) {
      auc <- pROC::roc(
        response = preds[, "group"],
        predictor = preds[, "pred"],
        direction = ">",
        levels = group.values,
        quiet = TRUE
      )$auc[1]
      D <- cbind(D, auc)
    }

    if (pred.prob) {
      # calculate probability
      preds$P <- exp(preds$pred) / (1 + exp(preds$pred))
      # cutoffs frequencies
      preds$cut.groups <-
        cut(preds$P,
          breaks = prob.cutoffs,
          include.lowest = TRUE, right = FALSE
        )
      # probability table
      P.table <-
        prop.table(
          table(
            as.character(preds$group),
            preds$cut.groups
          ),
          margin = 1
        )
    } else {
      P.table <- NULL
    }

    comb.output <- list(
      D = D,
      pred.dat = preds,
      cv.mod = cv.mod,
      P.table = P.table
    )
    return(comb.output)
  }
