colwise_pool <- function(data, n1, n2, m1, m2, sd1, sd2) {
  n1 <- data[, n1]
  n2 <- data[, n2]
  m1 <- data[, m1]
  m2 <- data[, m2]
  sd1 <- data[, sd1]
  sd2 <- data[, sd2]

  numerator1 <- (n1 - 1) * sd1^2
  numerator2 <- (n2 - 1) * sd2^2
  denominator1 <- (sum(n1) - length(n1))
  denominator2 <- (sum(n2) - length(n2))

  pooled.sd.1 <- sqrt(sum(numerator1) / denominator1)
  pooled.sd.2 <- sqrt(sum(numerator2) / denominator2)

  total.n1 <- sum(n1)
  total.n2 <- sum(n2)

  pooled.sd.total <- sqrt(((total.n1 - 1) * pooled.sd.1^2 +
    (total.n2 - 1) * pooled.sd.2^2) /
    (total.n1 + total.n2 - 2))
  output <- cbind(pooled.sd.1, pooled.sd.2, pooled.sd.total)
  return(output)
}
