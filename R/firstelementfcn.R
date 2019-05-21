#' It returns the first element non-NA from a vector. 
#' @description Internal function to use in annotation. 
#' @param y vector of characters
#' @examples firstelementfcn(LETTERS[1:10])
firstelementfcn<-function(y){ y<-uniquenona(y); y<-ifelse(length(y)==0,'',y[1])}
