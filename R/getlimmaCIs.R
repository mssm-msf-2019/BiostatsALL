#' Functions to get CI from the limma model
#' @description Functions to get  CI from the limma model

getlimmaCIs<-function(fit){
  width.ci<-qt(0.975, fit$df.residual + fit$df.prior) * fit$stdev.unscaled * sqrt(fit$s2.post)
  out<-list(Lower.CI=(fit$coef-width.ci),Upper.CI=(fit$coef+width.ci), SE=fit$stdev.unscaled * sqrt(fit$s2.post))
  return(out)
}
