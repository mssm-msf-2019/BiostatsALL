#' Function given a collection of pathways, returns the summary ( score for)mean, CI, p value etc) for multiple comparisons at once
#' @description Given a collection of pathwyas in SETs and the value of an estimate (ie lg FCH) in each column of STATs, it returns the summary
#' @param SETs matrix of number_of_genes by number_of_pathways. TRUE of FALSE if the gene belong to the pathway
#' @param STATS matrix of number_of_genes by number_of_comparisons with the estimates
#' @param fname name of the file to save results
#' @examples
#' z=getMultiplMeanCI.by.Path(ebfit$coef,SETs,'MeanFCH4Comparisons')

getMultipleMeanCI.by.Path<-function(STATs,SETs,fname){
  L<-lapply(STATs,function(y,set){getMeanCI.by.Path(set,y)}, SETs)
  if (!is.null(fname)){
    for (i in c(1:ncol(STATs))) {
      write.csv(file=paste(fname,colnames(STATs)[i], '.csv',sep=''),L[[i]])
    }
  }
}