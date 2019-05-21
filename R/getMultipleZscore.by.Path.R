#' Function given a collection of pathways, returns the Z score for multiple comparisons at once
#' @description Given a collection of pathwyas in SETs and the value of a T statistic in each column of T, it returns the Z score in a matrix of number_of_pathwyas by number_of_comparisons
#' @param SETs matrix of number_of_genes by number_of_pathways. TRUE of FALSE if the gene belong to the pathway
#' @param STATS matrix of number_of_genes by number_of_comparisons with the T statistics
#' @param fname name of the file to save results
#' @examples
#' z=getMultipleZscore.by.Path(ebfit$t,SETs,'Zscore4Comparisons')
#' 
#' 

getMultipleZscore.by.Path<-function(STATs,SETs,fname){
  L<-lapply(STATs,function(y,set){getZscore.by.Path(set,y)}, SETs)
  if (!is.null(fname)){
    for (i in c(1:ncol(STATs))) {
         write.csv(file=paste(fname,colnames(STATs)[i], '.csv',sep=''),L[[i]])
    }
  }
}