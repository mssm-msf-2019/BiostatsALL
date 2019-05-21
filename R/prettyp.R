#' A function that express a number using scientific notation or an specific number of decimals
#' @description This function express a number using scientific notation if it is smaller than a given power of 10 or rounded to a given number of decimals
#' @param p Numeric object

prettyp<-function(p,dig=2){
  if (p<10^(-dig)) {out<-sciNotation(p,digits=1)} else {out<-as.character(round(p,dig))}
  return(as.character(out))
}