#' Functions to get p-value from Kruskal-Wallis test
#' @description Functions to get p-value from Kruskal-Wallis test, used by createNiceTable()

printKruskalPVal <- function(fmla, data) {
  k <- kruskal.test(fmla, data)
  printPVal(k$p.value)
}
