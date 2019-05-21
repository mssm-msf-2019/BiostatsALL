#' A function that changes the column names of a matrix by adding a suffix function unique does keep NA in the output
#' @description A function that changes the column names of a matrix by adding a suffix function unique does keep NA in the output
#' @param mat: matrix
#' @param suffix: character
#' @examples 
#' givecolnames(ebfit$coef,'lgFCH')
#' givecolnames(ebfit$p,'p')
givecolnames<-function(mat,suffix){
  if (length(grep('FCH',suffix))>0) {mat<-round(mat,2)}
  colnames(mat)<-paste(suffix,colnames(mat),sep='_'); 
  return(mat)}