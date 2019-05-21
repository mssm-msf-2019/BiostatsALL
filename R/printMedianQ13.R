#' Functions to get median and IQR
#' @description Functions to get median and IQR used by createNiceTable()

printMedianQ13 <- function(vec) {
  q <- quantile(vec, na.rm=TRUE)
  sprintf("%.1f [%.1f-%.1f]", median(vec,na.rm=TRUE), q["25%"], q["75%"])
}

