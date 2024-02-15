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
           prob.cutoffs = seq(from = 0, to = 1, by = 0.20),
           append.data = FALSE) {

    if (is.null(size)) {
      size <- round(nrow(data) / 4, 0)
    } else {
      size <- size
    }

    if (is.data.frame(data)){
      data$group.var.num <-
        ifelse(data[, group.var] == group.values[1], 1,
               ifelse(data[, group.var] == group.values[2], 0,
                      NA
               )
        )

      data$row.nmbr <- rownames(data)

      data.grouped <- dplyr::group_by(data, group.var.num)

      train.data <- dplyr::sample_n(data.grouped,
                                    size = size,
                                    replace = FALSE
      )

      test.data <- data[!(data$row.nmbr %in% train.data$row.nmbr), ]
      train.data <- dplyr::ungroup(train.data)

    } else {
      train.data <- data[[1]]
      test.data <- data[[2]]

      train.data$group.var.num <-
        ifelse(train.data[, group.var] == group.values[1], 1,
               ifelse(train.data[, group.var] == group.values[2], 0,
                      NA
               )
        )

      test.data$group.var.num <-
        ifelse(test.data[, group.var] == group.values[1], 1,
               ifelse(test.data[, group.var] == group.values[2], 0,
                      NA
               )
        )
    }

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

    if (append.data) {
      preds <- cbind(preds, test.data)
    }

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
