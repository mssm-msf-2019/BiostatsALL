#' Functions to get p-value from Chi-sq
#' @description Functions to get p-value from Chi-sq, used by createNiceTable()

printChiSqPVal <- function(fmla, data) {
  tbl <- ftable(fmla, data)
  Xsq <- chisq.test(tbl)
  printPVal(Xsq$p.value)
}
