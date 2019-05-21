#' A function that assigns levels of gray to a vector
#' @description This function creates a vector with differents levels of gray in concordance with the levels of a given vector 
#' @param g vector
#' @examples 
#' getgray(c("mild","moderate","severe"))

getgray<-function(g){
  n<-length(levels(g)); cc<-rep('',n)
  inc<-1/(n-1)
  for (i in c(1:n)){cc[(g==levels(g)[i])]<-gray((i-1)*inc)}
  return(cc)
}