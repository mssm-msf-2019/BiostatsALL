#' A function to draw Venn Diagrams from Decision Matrices
#' @description a function that is supposed to draw a heatmap based on a matrix of gene expressions
#' @param theD is the decision matrix that indicates if genes are differentially expressed
#' @param ps is the vector with identification for the probesets
#' @param fname is the prefix for the file to save the output
#' @param shape is the parameter that defines the shape for the Venn diagram figure
#' @param circle.names is a vector with the identifiers for the circles
#' @param show is a list if some parameters for the Vennerable Venn Diagram
#' @param sign if =TRUE will rpduce Venn for Up/Down regulated geens separately
#' @examples
#' doVenn(theD = D, ps = probeset, fname = 'filename', shape='circles', circle.names=NULL, show = list(Faces=FALSE,DarkMatter=FALSE))

doVenn<-function (theD, ps, fname, shape = "circles", circle.names = NULL,
                  show = list(Faces = FALSE, DarkMatter = FALSE), sign=FALSE)
{
  library(Vennerable)
  L1 = L2 = v1 = v2 = NULL
  pdf(file = paste(fname, ".pdf", sep = ""), onefile = TRUE)
  L <- apply(abs(theD) > 0, 2, function(x, ps) {
    ps[x]
  }, ps)
  if (class(L) == "matrix")
    L <- as.list(as.data.frame(L))
  if (!is.null(circle.names))
    names(L) <- circle.names
  else names(L) <- paste(colnames(theD))
  v <- VennFromSets(L)
  plot(v, doWeights = TRUE, doEuler = TRUE, type = shape, show = show)
  if (sign){
    L1 <- apply((theD) > 0, 2, function(x, ps) {
      (ps[x])
    }, ps)
    if (class(L1) == "matrix")
      L1 <- as.list(as.data.frame(L1))
    names(L1) <- paste(colnames(theD), "(Up)")
    L1 <- L1[which(sapply(L1, length) > 0)]
    if (length(L1) > 1) {
      if (!is.null(circle.names))
        names(L1) <- paste(circle.names, "(Up)")
      v1 <- VennFromSets(L1)
      plot(v1, doWeights = TRUE, doEuler = TRUE, type = shape,
           show = show)
    }
    L2 <- apply(theD < 0, 2, function(x, ps) {
      ps[x]
    }, ps)
    if (class(L2) == "matrix")
      L2 <- as.list(as.data.frame(L2))
    names(L2) <- paste(colnames(theD), "(Down)")
    L2 <- L2[(sapply(L2, length) > 0)]
    if (length(L2) > 1) {
      if (!is.null(circle.names))
        names(L2) <- paste(circle.names, "(Down)")
      v2 <- VennFromSets(L2)
      plot(v2, doWeights = TRUE, doEuler = TRUE, type = shape,
           show = show)
    }
  }
  dev.off()
  return(invisible(list(L1, L2, v1, v2)))
}
