#' Functions to get p-value from  Wilcoxon test
#' @description Functions to get p-value from Wilcoxon test, used by createNiceTable()

printWilcoxPVal <- function(fmla, data) {
  k <- wilcox.test(fmla, data, paired=FALSE)
  printPVal(k$p.value)
}

