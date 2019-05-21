#' A function that given a matriz of log2FCH calculates the (signed) FCH 
#' @description A function that reads GMT files and returns a list where each element is a list of genes
#' @param coefs.lg: matrix of log fold changes, usually the estimated coeficients of the limma model (ebfit$coef)
#' @examples 
#' A=rnorm(cbind(runif(10,-2,2),runif(10,-1,1)))
#' calcFCH(A)

calcFCH<-function(coefs.lg){
  coefs.fch<-apply(coefs.lg,2,function(x){sign(x)*2^(abs(x))})
  colnames(coefs.fch)<-paste('FCH',colnames(coefs.lg),sep='_')
  coefs.fch<-round(coefs.fch,2)
  return(coefs.fch)
}