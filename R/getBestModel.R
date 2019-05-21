
#' A function to compare two model fits 
#' @description Compare 2 models based on AIC and model complexity, application can be to compare two mixed model with same parameters but different correlartion structure
#' @param fit1: height(cm) and weight(kg)
#' @param fit2: height(cm) and weight(kg)
#' @examples 
#' 
getBestModel<-function(fit1,fit2){
  aic<- AIC(fit1,fit2)
  anova.out<-try(anova(fit1,fit2))
  if (any(colnames(anova.out)=='p-value')&(class(anova.out)!="try-error")){
    if (anova.out[2,'p-value']<=0.05) {fit<-list(fit1,fit2)[[which.min(aic$AIC)]] }
    else {fit<-list(fit1,fit2)[[which.min(aic$df)]]}
  } else  {fit<-list(fit1,fit2)[[which.min(aic$AIC)]]}
  return(fit)
}

