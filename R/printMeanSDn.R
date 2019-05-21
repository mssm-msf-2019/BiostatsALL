#' Functions to get mean and SD with different Ns
#' @description Functions to get mean and SD used by createNiceTable()

printMeanSDn <- function (vec) {
  sprintf("%.1f (%.1f), %.0f", mean(vec, na.rm = TRUE), sd(vec, na.rm = TRUE), sum(!is.na(vec)))
}
