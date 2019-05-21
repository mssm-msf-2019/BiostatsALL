#' A function that creates ranked files for GSEA 
#' @description a function that eliminates ambiguity in gene symbols and prepares preRanked Gsea files for different contrasts
#' depends on getMat and myfindLargest
#' @param coefs a matrix with coefficients of a linear model
#' @param ann is the package with annotation for the platform that is being used
#' @param fname is the prefix for the files to be output
#' @examples 
#' prepareGSEArnk(coefs = coefs, Annpkg = hgu133plus2.db, fname='mycontrasts',symbs='NULL')


Create_GSEAInput<-function(coefs,ann,fname){
  for (i in c(1:ncol(coefs))){
    imatr<-getMat(coefs[,i],ann)
    write.table(file=paste(fname,colnames(coefs)[i],'_2GSEA.rnk',sep=''), x=as.matrix(imatr), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
    a<-myfindLargest(ann[rownames(coefs),'Symbol'],coefs[,i],ann[rownames(coefs),'Symbol'])
    imatr<-getMat(coefs[a,i],ann)
    write.table(file=paste(fname,colnames(coefs)[i],'_2GSEA_Unique.rnk',sep=''), x=as.matrix(imatr), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
}