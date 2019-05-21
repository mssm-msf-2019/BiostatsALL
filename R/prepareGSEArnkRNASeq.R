#' A function that prepares GSEA ranked data
#' @description a function that prepares preRanked GSEA files for different contrasts for RNASeq Data
#' @param coef a matrix with coefficients of a linear model
#' @param Annpkg.db is the package with annotation for the platform that is being used
#' @param fname is the prefix for the files to be output
#' @param ATab is the annotation table
#' @examples
#' ps<-rownames(D2)[which( (D2[,'A.LSvsNormal'] != 0) )]
#' prepareGSEArnkRNASeq(coefs = cbind(A.LSvsN=ebfit2$coef[ps,c("A.LSvsNormal")]),
#'                      fname='../RNASeq Analysis/1. Excluding Pools  9 10/GSEA/GSEA.coef',
#'                      ATab=TableAnn[ps,])

prepareGSEArnkRNASeq<-function (coefs, fname, ATab = NULL) {
  coefs <- as.matrix(coefs)
  if (is.null(ATab)) {
    print("Please provide annotation")
  } else {
    ann <- as.data.frame(subset(ATab, ensembl_gene_id %in% rownames(coefs)))[,c("ensembl_gene_id", "external_gene_id")]
  }
  coef.names <- colnames(coefs)
  coefs <- as.matrix(coefs[as.character(ann$ensembl_gene_id), ])

  for (i in 1:ncol(coefs)) {
    print(i)
    o <- order(coefs[, i], decreasing = TRUE)
    mat <- cbind(ann[o, 1:2], d = coefs[o, i])
    write.table(file = paste(fname, "_PS_", coef.names[i],
                             ".rnk", sep = ""), x = as.matrix(mat[, -2]), sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    mat <- mat[which(mat[, "external_gene_id"] != ""), ]
    gene.vector<-as.character(mat[, "external_gene_id"])
    names(gene.vector)<-as.character(rownames(mat))
    l <- unlist(myfancy_findLargest(probe.ids=mat[, "ensembl_gene_id"],
                             testStat=abs(as.numeric(mat[,"d"])),
                             gene.names=gene.vector))
    mat2 <- mat[as.vector(unlist(l)), ]
    o <- order(as.numeric(mat2[, "d"]), decreasing = TRUE)
    write.table(file = paste(fname, "_Gene_", coef.names[i],
                             ".rnk", sep = ""), x = as.matrix(mat2[o, -1]), sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
