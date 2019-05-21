#' Functions to get mean and range (Min and Max)
#' @description Functions to get mean and Range (Min Max) used by createNiceTable()

printMeanMinMax <- function(vec) {
  sprintf("%.1f [%.1f-%.1f]", mean(vec,na.rm=TRUE), range(vec,na.rm=TRUE)[1], range(vec,na.rm=TRUE)[2])
}












# Lewis made this function, he is awesome
