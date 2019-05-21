#' A function that given a vector returns only the unique elements that are not NA or epty characters ("").
#'  function unique does keep NA in the output
#' @description A function that given a vector returns only the unique elements that are not NA or epty characters (""). Function unique does keep NA in the output
#' @param x: vector of characters
#' @examples 
#' uniquenona(c(LETTERS[1:10],"",NA))

uniquenona<-function(x){
  x<-x[!is.na(x)]; 
  if (is.character(x)) x<-x[!(x=='')]; 
  if (length(x)==0) y="" else y<-unique(x)
  return(y)
}