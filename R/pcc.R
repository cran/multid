#' Returns probabilities of correct classification for both groups in independent data partition.
#'
#' @param data Data frame including predicted values (e.g., pred.dat from D_regularized_out).
#' @param pred.var Character string. Variable name for predicted values.
#' @param group.values Vector of length 2, group values (e.g. c("male", "female) or c(0,1)).
#' @param group.var The name of the group variable.
#'
#' @return Vector of length 2. Probabilities of correct classification.
#' @export
#'
#' @examples
#' D_out <- D_regularized_out(
#'   data = iris[iris$Species == "versicolor" | iris$Species == "virginica", ],
#'   mv.vars = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   group.var = "Species", group.values = c("versicolor", "virginica"),
#'   size = 15
#' )
#'
#' pcc(
#'   data = D_out$pred.dat,
#'   pred.var = "pred",
#'   group.var = "group",
#'   group.values = c("versicolor", "virginica")
#' )
pcc <- function(data,
                pred.var,
                group.var,
                group.values) {
  pcc.1 <-
    sum(data[data[, group.var] == group.values[1], pred.var] > 0) /
      length(data[data[, group.var] == group.values[1], pred.var])

  pcc.2 <- sum(data[data[, group.var] == group.values[2], pred.var] < 0) /
    length(data[data[, group.var] == group.values[2], pred.var])

  pcc.total <- (sum(data[data[, group.var] == group.values[1], pred.var] > 0) +
    sum(data[data[, group.var] == group.values[2], pred.var] < 0)) /
    nrow(data)

  pcc.out <- cbind(
    pcc.1,
    pcc.2,
    pcc.total
  )

  return(pcc.out)
}
