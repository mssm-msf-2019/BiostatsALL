#' Functions to get p-value from Pooled t-test (equal variance)
#' @description Functions to get p-value from Pooled t-test (equal variance), used by createNiceTable()

printTtestPoolPVal <- function(fmla, data) {
  fit <- t.test(fmla, data=data, var.equal = TRUE)
  pv <- fit$p.value
  printPVal(pv)
}

