#' Functions to convert list object to a matrix
#' @description Functions to convert list object to a matrix

List2mat<-function(list){
  out<-t(do.call(rbind,list))
  colnames(out)<-names(list)
  return(out)
}
