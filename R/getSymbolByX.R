#' A function that 
#' @description This function 
#' @param Pmat numeric matrix
#' @param X expression set
#' @param descend logical object. The default is TRUE. If TRUE the output will be sort on decreassing order.  

getSymbolByX<-function(Pmat,X,descend=TRUE){
  genebyX<-rep(NA,ncol(Pmat));names(X)<-colnames(Pmat)
  for (i in c(1:ncol(Pmat))){
    g<-which(Pmat[,i])
    og<-order(as.numeric(X[g]),decreasing=descend)
    genebyX[i]<-stackchar(X=paste(names(g)[og],'(',as.character(round(X[g][og],2)),')',sep=''),sep=', ')
  }
  return(genebyX)
}
