
#' Function for logistic regression with predictors that has more than 2 levels. 
#' @description predictor is independent descrete variable that has more than 2 levels.
#' 


## for those discrete predictor with more than 2 levels
RunLogitXwith3levels<- function(predictor, data, outcome, covariates)
{
  library(mvtnorm)
  library(emmeans)
  library(multcompView)
  
  db <- subset(data, select = c(outcome,covariates,predictor))
  frml <- paste0(outcome,'~',paste0(c(covariates,predictor),collapse = '+'))
  lsmean.frml<-paste0('~',paste0(c(predictor),collapse = '+'))
  
  f <- glm(frml,family = binomial, data = db)
  lsm<-emmeans(f,as.formula(lsmean.frml), adjust='none')
  plot(lsm,comparison=T)
  
  Cst<-as.data.frame(contrast(lsm,method = "revpairwise",adjust='none'))[,c('contrast','estimate','SE','p.value')]
  Cst.CI<-as.data.frame(confint(contrast(lsm,method = "revpairwise",adjust='none')))[,c('contrast','asymp.LCL','asymp.UCL')]
  
  af<-as.data.frame(anova(f,test='Chisq'));
  af<-af[grep(predictor,rownames(af),value=T),grep('pr',colnames(af),ignore.case = T)]
  Res<-merge(Cst,Cst.CI,by='contrast')
  Res<-cbind(predictor=predictor,Res,Anova.p=af)
  return(Res)
} 