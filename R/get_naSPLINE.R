
#' Functions to replace each NA with interpolated values
#' @description Functions to replace each NA with interpolated values
#' @param x Variables to be used for interpolation 
#' @examples 
#' z <- c(11, NA, 13, NA, 15, NA)
#' get_naSPLINE(z)

get_naSPLINE<-function(x){SPLINE= na.spline(x)
                               return(SPLINE)}
