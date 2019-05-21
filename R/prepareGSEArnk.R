#' A function that prepares GSEA ranked data
#' @description a function that prepares preRanked Gsea files for different contrasts
#' @param coef a matrix with coefficients of a linear model
#' @param Annpkg is the package with annotation for the platform that is being used
#' @param fname is the prefix for the files to be output
#' @param symbs is a vector with genes symbols
#' @examples
#' prepareGSEArnk(coefs = coefs, Annpkg = hgu133plus2.db, fname='mycontrasts',symbs=NULL)


prepareGSEArnk<-function(coefs, Annpkg.db=NULL, fname, ATab=NULL){
  coefs<-as.matrix(coefs)

  if (is.null(ATab)) {
    library(Annpkg.db$packageName, character.only=TRUE)
    ann<-select(Annpkg.db, rownames(coefs), c('SYMBOL'))
  } else {
    ann<-as.data.frame(subset(ATab, PROBEID%in% rownames(coefs)))[,c('PROBEID','SYMBOL')]
  }


  coefs<-as.matrix(coefs[as.character(ann$PROBEID),])
  for (i in 1:ncol(coefs)){
    print(i)
    o<-order(coefs[,i],decreasing=TRUE)
    mat<-cbind(ann[o,1:2],d=coefs[o,i]);

    write.table(file=paste(fname,'_PS_',colnames(coefs)[i],'.rnk',sep=''),x=as.matrix(mat[,-2]),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

    mat<-mat[which(mat[,'SYMBOL']!=""),]
    #l<-myfindLargest(mat[,'PROBEID'],-abs(as.numeric(mat[,'d'])),mat[,'SYMBOL'])
    #mat2<-mat[as.vector(l),]

    l<-myfancy_findLargest(mat[,'PROBEID'],abs(as.numeric(mat[,'d'])),mat[,'SYMBOL'])
    mat2<-mat[as.vector(unlist(l)),]
    o<-order(as.numeric(mat2[,'d']),decreasing=TRUE)
    write.table(file=paste(fname,'_SYMB_',colnames(coefs)[i],'.rnk',sep=''),
                x=as.matrix(mat2[o,-1]),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  }

}
