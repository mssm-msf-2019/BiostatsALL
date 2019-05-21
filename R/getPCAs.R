#' Extraction of principal components
#' @description a function that extracts principal components from an expression set and output the scores of most important components and variance explained
#' the original function would output all results in a data frame. This function outputs to a list
#' @param eset is an expression set
#' @param maxpc is the maximum number of principal components to be extracted
#' @examples 
#' getPCAs(eset,maxpc=3)

getPCAs<-function(eset,maxpc=3){
  if (!is.matrix(eset)) {ex<-exprs(eset); out<-pData(eset)} else {out<-NULL; ex<-eset}
  pca.res <- prcomp(t(ex))
  x<-pca.res$x[,1:maxpc]; colnames(x)<-paste('PC',1:maxpc,sep='.')
  varex=round(100*summary(pca.res)$importance[2,])
  db<-cbind.data.frame(x)
  if (!is.null(out)) db<-cbind(db,out)
  return(list(db=as.data.frame(db),varex=varex))
}
