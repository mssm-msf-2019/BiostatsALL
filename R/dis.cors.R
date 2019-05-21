#' Internal function of DoCorrelationPlotwithStarts. 
#' @description Distance matrix based on correlations. 
#' @param x: symetrix matrix indicating the  correlation coefficient for each column/row pair of variables
#' @param method: correlation method (spearman, perason etc, default to spearman)
#' @param abs: should the distance be 1-cor (FALSE,default) or 1-|cor| (TRUE). In teh later case, variables negative correlated will be as near as positive correlated
#' @examples
#' # in this examples, I want to cluster samples
#' a<-dis.cors(cor(exprs(reset)[1:10,]), abs=TRUE,method='pearson')

dis.cors <- function(x, abs=FALSE, method = "spearman",...) {
  cm<-cor(as.matrix(t(x)), method = method, use = "p")
  if (abs){ cm<-abs(cm)}
  res <- as.dist((1-cm)/2)
  # attr(res, "method") <- method
  return(res)
}
