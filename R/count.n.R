
#' counting the non-missing value
#' @description Given a vector with missing value, return the number of non-missing
#' @param x: a vector

count.n<-function(x){sum(!is.na(x))}
