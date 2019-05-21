#' Functions to build pretty-prints p values
#' @description Functions to build pretty-prints p values, used by createNiceTable()

printPVal <- function(pval) {
  if (is.na(pval)) {
    "NA"
  } else {
    log10pval <- log10(pval)
    if (log10pval < -5) {
      sprintf("< 0.00001")
    } else if (log10pval < -4) {
      sprintf("< 0.0001")
    } else if (log10pval < -3) {
      sprintf("< 0.001")
    } else {
      sprintf("%0.3f", round(pval, digits=3))
    }
  }
}
