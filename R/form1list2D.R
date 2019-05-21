#' A function that converts a list of genelist into a Decision matrix
#' @description Converts a list of L vectors containing features into a binary indicator matrix of nxL where L is the length of the universe of genes in the List
#' @param list: list of L vectors containing gene lists 
#' @examples 
#' L<-readGmt('KruegerSYMBOL2015.gmt')
#' D<-from1list2D(L)
from1list2D<-function(list){
  univ<-unlist(list)
  D<-matrix(0,nrow=length(univ),ncol=length(list),dimnames=list(univ,names(list)))
  for (i in c(1:length(list))){D[univ,i]<-1*(univ%in%list[[i]])}
  return(D)
}	