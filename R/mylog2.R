#' A function that returns the base 2 logarithm of the components of a numeric vector
#' @description This function given a numeric vector calculates the base 2 logarithm of the components
#' @param x numeric vector

mylog2<-function(x){y<-x; y[(x==0)]<-min(x[(x!=0)],na.rm=TRUE)/5; return(log2(y))}
