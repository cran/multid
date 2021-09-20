#' Multivariate group difference estimation with regularized binomial regression
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
#' @param out Logical. Should results and predictions be calculated on out-of-bad data set? (Default FALSE)
#' @param size Integer. Number of cases in regularization data per each group. Default 1/4 of cases.
#' @param fold Logical. Is regularization applied across sample folds with separate predictions for each fold? (Default FALSE)
#' @param fold.var Character string. Name of the fold variable. (default NULL)
#'
#' @return
#' \item{D}{Multivariate descriptive statistics and differences.}
#' \item{pred.dat}{A data.frame with predicted values.}
#' \item{cv.mod}{Regularized regression model from cv.glmnet.}
#' @seealso \code{\link[glmnet]{cv.glmnet}}
#' @export
#'
#' @examples D_regularized(
#'   data = iris[iris$Species == "setosa" | iris$Species == "versicolor", ],
#'   mv.vars = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   group.var = "Species", group.values = c("setosa", "versicolor")
#' )$D
#'
#' # out-of-bag predictions
#' D_regularized(
#'   data = iris[iris$Species == "setosa" | iris$Species == "versicolor", ],
#'   mv.vars = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   group.var = "Species", group.values = c("setosa", "versicolor"),
#'   out = TRUE, size = 15
#' )$D
#'
#' # separate sample folds
#' # generate data for 10 groups
#' set.seed(34246)
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
#' # Fit and predict with same data
#' D_regularized(
#'   data = d,
#'   mv.vars = c("x1", "x2", "x3"),
#'   group.var = "sex",
#'   group.values = c("female", "male"),
#'   fold.var = "fold",
#'   fold = TRUE,
#'   rename.output = TRUE
#' )$D
#'
#' # Out-of-bag data for each fold
#' D_regularized(
#'   data = d,
#'   mv.vars = c("x1", "x2", "x3"),
#'   group.var = "sex",
#'   group.values = c("female", "male"),
#'   fold.var = "fold",
#'   size = 17,
#'   out = TRUE,
#'   fold = TRUE,
#'   rename.output = TRUE
#' )$D
D_regularized <-
  function(data,
           mv.vars,
           group.var,
           group.values,
           alpha = 0.5,
           nfolds = 10,
           s = "lambda.min",
           type.measure = "deviance",
           rename.output = TRUE,
           out = FALSE,
           size = NULL,
           fold = FALSE,
           fold.var = NULL) {

    # out-of-bag and folds (fold_out)?
    if (out & fold) {
      D_regularized_fold_out(
        data = data,
        mv.vars = mv.vars,
        group.var = group.var,
        group.values = group.values,
        alpha = alpha,
        s = s,
        type.measure = type.measure,
        rename.output = rename.output,
        size = size,
        fold.var = fold.var
      )
    } # out-of-bag and no folds (out)?
    else if (out & !fold) {
      D_regularized_out(
        data = data,
        mv.vars = mv.vars,
        group.var = group.var,
        group.values = group.values,
        alpha = alpha,
        size = size,
        s = s,
        type.measure = type.measure,
        rename.output = rename.output
      )
    } # not out-of-bag and folds (fold)?
    else if (!out & fold) {
      D_regularized_fold(
        data = data,
        mv.vars = mv.vars,
        group.var = group.var,
        group.values = group.values,
        alpha = alpha,
        s = s,
        type.measure = type.measure,
        rename.output = rename.output,
        fold.var = fold.var
      )
    } # vanilla version
    else {
      D_regularized_vanilla(
        data = data,
        mv.vars = mv.vars,
        group.var = group.var,
        group.values = group.values,
        alpha = alpha,
        s = s,
        type.measure = type.measure,
        rename.output = rename.output
      )
    }
  }
