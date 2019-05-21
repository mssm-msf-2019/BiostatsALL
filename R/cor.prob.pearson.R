#' calculates the p-values for the pair-wise pearson correlation coefficient (called by CorrelationHeatmap)
#' @description Given a matrix of observations by variables, calculates the pearson correlation between the variable(columns). columns must be named. It returns a simmetric matrix of number of variables, with p values for the correlation between row and column.
#' @param x: matrix of observations by variables
#' @examples
#' cor.prob.pearson(db)

cor.prob.pearson<-function (x)
{
  FUN <- function(x, y) cor.test(x, y, use='p')[[3]]
  R <- outer(colnames(x), colnames(x), Vectorize(function(i, j) FUN(x[, i], x[, j])))
  dimnames(R) <- list(colnames(x), colnames(x))
  R[row(R) == col(R)] <- NA
  return(R)
}
