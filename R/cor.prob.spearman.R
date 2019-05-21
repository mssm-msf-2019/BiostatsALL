#' calculates the p-values for the pair-wise spearman correlation coefficient (called by CorrelationHeatmap)
#' @description Given a matrix of observations by variables, calculates the spearman correlation between the variable(columns). columns must be named. It returns a simmetric matrix of number of variables, with p values for the correlation between row and column.
#' @param x: matrix of observations by variables
#' @examples
#' cor.prob.pearson(db)

cor.prob.spearman<- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs",method="spearman")
  r2 <- R^2
  Fstat <- r2 * dfr/(1 - r2)
  R<- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  return(R)
}
