#' A function that replace values in a vector 
#' @description This function replaces in a vector the values that satisfice a given pattern by other values also given
#' @param pattern vector of values
#' @param replacement vector of values
#' @param x vector of values
#' @examples 
#' mgsub(c(0,1),c("No","Yes"),x=c(1,1,1,0,1,0,0,1,0,1))
 

mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}