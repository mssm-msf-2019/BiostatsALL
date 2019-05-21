#' Functions to get mean and SD
#' @description Functions to get mean and SD used by createNiceTable()

printMeanSD <- function(vec) {
  sprintf("%.1f (%.1f)", mean(vec,na.rm=TRUE), sd(vec,na.rm=TRUE))
}
