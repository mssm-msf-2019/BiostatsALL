#' creates an output table including pathways that are declared differentially expresed (given the user-speceified cut-offs) in at least one of the comparisons embeded in ebfit. Like FromEbfit2DEGTable, but does not include FCH columns as they dont make sense in the GSVA context.
#' @description For each contratst in ebfit, a matrix with lgFCH,pvalue , fdr and Status (1 up,-1 doen,0 non sign) is created using pre-specified cut-oofs.
#' @param ebfit output of eBayes
#' @param adj adjustment method for multiple hypothesis
#' @param mcut cutoff for FCH (to be used for Status columns)
#' @param pcut cutoff for pvalue/FDR  (to be used for Status columns)
#' @param annot should the dat be annotated? If TRUE it will check if tehre is annotTab or otherwise whatever annotation is stored in ebfit$genes
#' @param annotTab annotation database
#' @param first.ann.cols Indicates which elements of the annotation db, should appear before all the comparisons.
#' @examples
#' TabDEG<-FromEbfit2DEGTable_GSVA(ebfit,pcut=0.05,mcut=2)
#'
#'

FromEbfit2DEGTable_GSVA<-function(ebfit,adj='BH',mcut=2,pcut=0.05, annot=TRUE, annotTab=NULL,first.ann.cols=c('SYMBOL','GENENAME')){
  Tab<-FromEbfit2Table_GSVA(ebfit,adj,mcut,pcut,annot=FALSE)
  StatusTab<-as.matrix(Tab[,grep('Status',colnames(Tab))])
  DEGs<-rownames(Tab)[which(rowSums(abs(StatusTab))>0)]
  TabDEGs<-Tab[DEGs,]
  return(TabDEGs)
}
