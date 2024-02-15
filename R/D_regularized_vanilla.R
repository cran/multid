D_regularized_vanilla <-
  function(data,
           mv.vars,
           group.var,
           group.values,
           alpha = 0.5,
           nfolds = 10,
           s = "lambda.min",
           type.measure = "deviance",
           rename.output = TRUE,
           append.data = FALSE) {
    data$group.var.num <-
      ifelse(data[, group.var] == group.values[1], 1,
        ifelse(data[, group.var] == group.values[2], 0,
          NA
        )
      )

    cv.mod <-
      glmnet::cv.glmnet(
        x = as.matrix(data[, c(mv.vars)]),
        y = data$group.var.num,
        family = c("binomial"),
        nfolds = nfolds,
        type.measure = type.measure,
        alpha = alpha
      )

    preds <- data.frame(
      group = data[, group.var],
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

    D <- multid::d_pooled_sd(
      data = preds,
      var = "pred",
      group.var = "group",
      group.values = group.values,
      rename.output = rename.output
    )

    comb.output <- list(
      D = D,
      pred.dat = preds,
      cv.mod = cv.mod
    )
    return(comb.output)
  }
