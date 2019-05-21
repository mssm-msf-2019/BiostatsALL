
#' Functions for Best Observation Carried Forward
#' @description for Best Observation Carried Forward
#' @param x Variables to be used for populates missing values with the subject's best-case nonmissing value.
#' @examples 
#' z <- c(11, NA, 13, NA, 15, NA)
#' get_naBestOCF(z)
#' 
get_naBestOCF<-function(x){BestOCF<-x;
       BestOCF[is.na(x)]<-min(x[!is.na(x)])
       return(BestOCF)}
