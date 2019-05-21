#' A function that converts a two class lists of up and down genes into a vector with value -1 if down,1 if up and zero otherwise
#' @description Input are the up and down-regulated list for a set of comparisons C1,...CL,the output will be a matrix of L columns, where 
#' each element if #' 1(up),0 (not),-1 (down). Useful to construct venn-diagram and status columns for some applications when you are reading lists 
#' from collections and you want to intersect with your analysis.
#' @param up: is a list of L vectors containing up regulated genes in each comparison C1,...CL
#' @param down: is a list of L vectors containing down regulated genes in each comparison C1,...CL
#' @examples 
#' up<-list(Set1=c('g1','g3'),Set2=c('g1','g4'))
#' down<-list(Set1=c('g5','g7'),Set2=c('g2','g6'))
#' D<-from2list2D(up,down)
from2list2D<-function(up,down){
  univ<-unique(c(unlist(up),unlist(down)))
  D<-matrix(0,nrow=length(univ),ncol=length(up),dimnames=list(univ,names(up)))
  for (i in c(1:length(up))){
    D[which(univ%in%down[[i]]),i]<-rep(-1,sum(univ%in%down[[i]])); 
    D[which(univ%in%up[[i]]),i]<-rep(1,sum(univ%in%up[[i]])); 
  }
  return(D)
}		 