#' Functions to get p-value from Satterthwaite t-test (un-equal variance)
#' @description Functions to get p-value from Satterthwaite t-test (un-equal variance), used by createNiceTable()

printTtestSattPVal <- function(fmla, data) {
  fit <- t.test(fmla, data=data, var.equal = FALSE)
  pv <- fit$p.value
  printPVal(pv)
}
