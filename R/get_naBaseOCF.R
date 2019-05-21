
#' Functions for baseline observation carried forward 
#' @description for baseline observation carried forward 
#' @param x Variables to be used for replacing NAs with the first observation 
#' @examples 
#' z <- c(11, NA, 13, NA, 15, NA)
#' get_naBaseOCF(z)
#' 
    get_naBaseOCF<-function(x){BaseOCF<-x;
                              BaseOCF[is.na(x)]<-first(x[!is.na(x)])
                              return(BaseOCF)}
    