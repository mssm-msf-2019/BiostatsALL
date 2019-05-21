#' A function that creates an Annotation table for microorrays experiments. 
#' @description A function that creates an Annotation table for microorrays experiments. It must be associated with an existing Array annotation data (ie hgu133plus2.db) where annotations has been  
#' assembled using data from public repositories
#' @param IDS: type of methods. Options are Regression, Classification adn Dual use
#' @param AnnPkg: named of the Array annotation data (ie, hgu133plus2.db)
#' @param w.annot: name of the information to be included in the table
#' @param w.annot.unique: which columns should be used to solved duplicated annotation; ie when a probe have two gene names.
#' @examples 
#' librarycolors()
#' 
#' 
#' AnnPkg <- hgu133plus2.db
#' getAnnotationTable(IDS=featureNames(eset), AnnPkg, w.annot=c("ENTREZID","SYMBOL","GENENAME"),w.annot.unique=c("ENTREZID"))

getAnnotationTable<-function(IDS, AnnPkg, w.annot=c("ENTREZID","SYMBOL","GENENAME","OMIM","CHRLOC","PATH","GO"), 
                             w.annot.unique=c("ENTREZID","SYMBOL","OMIM","CHRLOC")){
  require(AnnotationDbi)
  #   annotationTable <-AnnotationDbi::select(x=AnnPkg.db, keys=IDS, columns=unique(c(w.annot, w.annot.unique)))  ##This is the old version
  annotationTable <-AnnotationDbi::select(x=eval(as.name(AnnPkg)), keys=IDS, columns=unique(c(w.annot, w.annot.unique)))  ##This is the new version of the line, NOTE: the addition of "as.name" and "eval"
  annotationTable<-annotationTable[,-which(colnames(annotationTable)=='CHRLOC')]
  auxsub<-function(x){x[(x=="CHRLOC")]<-"CHRLOCCHR"; return(x)}
  w.annot.unique<-auxsub(w.annot.unique)
  w.annot<- auxsub(w.annot)
  annotationTabcollapsed<-collapseduplicatesids(annotationTable[,c('PROBEID', setdiff(w.annot,w.annot.unique))])
  annotationTabDupls<-treatduplicatesids(annotationTable[,c('PROBEID', w.annot.unique)])
  
  annotationTabCombine<-cbind(annotationTabDupls$First, 
                              annotationTabcollapsed,
                              renamecols(annotationTabDupls$Synon[, -1],suffix='Other'))
  annotationTabCombine<-annotationTabCombine[,!duplicated(colnames(annotationTabCombine))]
  return(as.data.frame(annotationTabCombine))
  
}

