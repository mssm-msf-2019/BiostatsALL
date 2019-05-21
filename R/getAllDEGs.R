#' Given a Decision matrix (output of decideTests), it returns the vestoor of probes that appears significant in at least one comparison
#' @description When several contrats are tested in limma, a call over if the gene is differentially Expressed is made after defining the FCH and FDR/p cutoffs
#' using the decideTestFunctions. This function takes that output and returns of vector that are DEG in at least 1 comparison. 
#' @param coefs.lg: matrix of log fold changes, usually the estimated coeficients of the limma model (ebfit$coef)
#' @examples 
#' D<-decideTests(ebfit,method="separate",adjust.method="BH",p.value=0.05,lfc=log2(2))
#' DEGs<-getAllDEGs(D)
getAllDEGs<-function(D){ rownames(D)[rowSums(D!=0)>0]}