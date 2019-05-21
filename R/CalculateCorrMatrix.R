#' Internal function of DoCorrelationPlotwithStarts. 
#' @description It calculates the correlation matrix and pvalues and returns a list with components. 
#' @param db: data base of continuos variables; columns are variables to be clustered.
#' @param method: correlation method (spearman, perason etc, default to spearman)
#' @param abs.cor: should the distance be 1-cor (FALSE,default) or 1-|cor| (TRUE). In teh later case, variables negative correlated will be as near as positive correlated
#' @examples
#' # in this examples, I want to cluster samples
#' a<-CalculateCorrMatrix(exprs(reset)[1:10,])
#' print(round(a$cor,2))
#' print(1*(a$p<=0.05))

CalculateCorrMatrix<-function(db,method,abs.cor=FALSE){
  db.cor<-cor(db,use='p',m= method)
  if (method=='spearman') db.pcor<-cor.prob.spearman(db) else db.pcor<-cor.prob.pearson(db)
  return(list(cor=db.cor,p=db.pcor))
}