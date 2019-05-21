#' Returns the number of missing values in a vector
#' @description Given a vector the function calculates the number of missing values
#' @param x vector

nmissing <- function(x) sum(is.na(x))
