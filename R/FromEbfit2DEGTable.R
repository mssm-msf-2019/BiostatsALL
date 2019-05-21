#' creates an output table including genes that are declared differentially expresed (given the user-speceified cut=offs) in at least one of the comparisons embeded in ebfit
#' @description For each contratst in ebfit, a matrix with lgFCH, FCH,pvalue , fdr and Status (1 up,-1 doen,0 non sign) is created using pre-specified cut-oofs
#' @param ebfit output of eBayes
#' @param adj adjustment method for multiple hypothesis
#' @param mcut cutoff for FCH (to be used for Status columns)
#' @param pcut cutoff for pvalue/FDR  (to be used for Status columns)
#' @param annot should the dat be annotated? If TRUE it will check if tehre is annotTab or otherwise whatever annotation is stored in ebfit$genes
#' @param annotTab annotation database
#' @param first.ann.cols Indicates which elements of the annotation db, should appear before all the comparisons.
#' @examples
#' TabDEG<-FromEbfit2DEGTable(ebfit,pcut=0.05,mcut=2)
#'
#'

FromEbfit2DEGTable<-function(ebfit,adj='BH',mcut=2,pcut=0.05, annot=TRUE, annotTab=NULL,first.ann.cols=c('SYMBOL','GENENAME')){
  Tab<-FromEbfit2Table(ebfit,adj,mcut,pcut,annot=FALSE)
  StatusTab<-as.matrix(Tab[,grep('Status',colnames(Tab))])
  DEGs<-rownames(Tab)[which(rowSums(abs(StatusTab))>0)]
  TabDEGs<-Tab[DEGs,]

  if (annot){
    if ((is.null(annotTab))&(!is.null(ebfit$genes))) {
      if (ncol(ebfit$genes)>1) {annotTab<-ebfit$genes[rownames(TabDEGs),] }
    }

    if(!is.null(annotTab)) {
      print(head(annotTab))
      imp.ann.col<-which(colnames(annotTab)%in%first.ann.cols)
      TabDEGs<-cbind(annotTab[rownames(TabDEGs), imp.ann.col],
                     TabDEGs, annotTab[rownames(TabDEGs),-imp.ann.col])
    }
  }

  return(TabDEGs)
}
