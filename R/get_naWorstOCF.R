
#' Functions for Worst Observation Carried Forward
#' @description for Worst Observation Carried Forward
#' @param x Variables to be used for populates missing values with the subject's Worst-case nonmissing value.
#' @examples 
#' z <- c(11, NA, 13, NA, 15, NA)
#' get_naWorstOCF(z)
 
get_naWorstOCF<-function(x){WorstOCF<-x;
                            WorstOCF[is.na(x)]<-max(x[!is.na(x)])
                            return(WorstOCF)}
