#' A function that creates genotypes for GSEA
#' @description a function that creates genotypes for GSEA
#' @param G a vector with the identification of the genotypes
#' @examples 
#'makeGSEApheno(G=D)

 
makeGSEApheno<-function(G){
  G<-as.character(G);  nG<-length(G)
  nf<-length(unique(G))
  phmat<-matrix('',3,nG)
  phmat[3,]<-G
  phmat[2,1]<-'#'
  phmat[2,2:(nf+1)]<-unique(G)
  phmat[1,1:3]<-c(nG,nf,1)
  return(phmat)
}