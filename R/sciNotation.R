#' A function that returns a scientific notation of a p value
#' @description This function creates an expression of a p value using scientific notation 
#' @param x Numeric object. Number to be expressed using scientific notation.
#' @param digits Numeric object of type integer. Number of digits. 
#' @param maxexplogical Numeric object of type integer. Maximun exponential admissible.
#' @examples 
#' sciNotation(x=0.00000000000068,digits=2,maxexp=16)

sciNotation<-function(x,digits=2,maxexp=16){
  
  if (length(x) > 1) {return(append(sciNotation(x[1]), sciNotation(x[-1])))} 
  if (is.na(x)) return('')
  if (x>10^(-maxexp)){
    exponent<-floor(log10(x))
    base <-round(x/(10^exponent), digits)
    if (base==1) out<-as.expression(substitute(10^exponent,list(base=base, exponent=exponent)))
    else out<-as.expression(substitute(base %*% 10^exponent,list(base=base, exponent=exponent)))
  }
  else {	  out<-expression(paste("p<",10^{-16},sep=''))}
  return(out)
}

