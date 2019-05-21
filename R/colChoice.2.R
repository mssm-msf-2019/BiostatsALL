#' Defining colors from a pattern
#' @description a function that returns a vector of color from a given pattern
#' @param colpat is the color pattern that will guide choice of colors
#' @param colfacs is the factor that should be colored according to the color pattern
#' @examples 
#' colChoice.2.R(colplat="blue",colfacs = myfactor)


colChoice.2 <- function(colpat="blue",colfacs){
  require(plyr)
  n=length(levels(colfacs))
  nrs <-topo.colors(n)
  colfacs_out<-mapvalues(colfacs, from = levels(colfacs), to = nrs)
  colfacs_out <- sapply(colfacs_out, as.character)
  return(colfacs_out)
}
