#' A function that returns an object of type character equal to 1, 2, 3 or 4 stars according to the value of a given object of type numeric
#' @description This function creates an object of type character equal to 1, 2, 3 or 4 stars (*) for diferents ranges of p values. To smaller p values correspond more number of stars. 
#' @param p numeric object. Usually a p value of a Test.
#' @param pmax numeric object. Usually the maximun permisible Type I Error of a Test.
#' @param parenthesis logical TRUE if parenthesis are required to wrap the stars

getstar<-function(p,pmax=0.05,parenthesis=FALSE){
  star<-' ';
  if (is.nan(p)) p<-1
  if (p<=0.1) star <-'*'
  if (p<=0.05) star <-'**'
  if (p<=0.01) star <-'***'
  if (p<=0.001) star <-'****'
  if (pmax==0.05) star=substr(star,2,nchar(star))
  if ((parenthesis)&(p<pmax)) star<-paste('(',star,')',sep='')
  return(star)
}