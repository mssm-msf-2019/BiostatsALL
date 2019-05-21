#' Function for Last Observation Carried Forward
#' @description Function for Last Observation Carried Forward: Replacingg each NA with the mostrecent non-NA prior to it. 
#' @param x Variables to be used for replacing each NA  
#' @examples 
#' z <- c(11, NA, 13, NA, 15, NA)
#' get_naLOCF(z)
 
get_naLOCF<-function(x){LOCF=zoo::na.locf(x)
                        return(LOCF)}
