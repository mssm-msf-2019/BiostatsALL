
#' convert specified values into missing 
#' @description Given a vector(could be factor), map specified value(s) into NA and return a vector
#' @param x: a vector(factor) with values that you want to convert into NA
#' @param values:a element/vector that will be converted as NA
#' @examples 
#' a<-c(999,0,1,99);mapvalues_toNA(a,999)
#' b<-c('salt', 'sugar','','pepper','unknown');mapvalues_toNA(b,c('','unknown'))

mapvalues_toNA<-function(x,values){
  x[(x%in%values)]<-NA;
  if (class(x)=='factor') x<-drop.levels(x)
  return(x)
}
