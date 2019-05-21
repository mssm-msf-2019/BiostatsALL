#' A function that returns one element per symbol, returning the largest of the probes that match each symbol
#' @description This function  Outputs a list of gene names with the probe set id with the maximum logFCH for that gene.
#' where if more than one probe set matches a gene and one ore more
#' of these probe sets is a cross-matchning one it is removed
#' prior to detecting the max logFCH
#' @param testStat : logFCH data for the specified coefficient
#' @param gene.names : The respective gene names of the testStat dataframe
#' @param probe.ids : The probeset ids for the testStat data frame
#' @examples 
#' a<-myfindLargest(myfindLargest(YOUR.RESULT.DF[,COEF.OF.INTEREST],YOUR.RESULT.DF[,"Symbol"],YOUR.RESULT.DF[,"Probe"]))


myfancy_findLargest<-function (probe.ids, testStat, gene.names){
 
  if (length(testStat) != length(probe.ids))
    stop("testStat and probe.ids must be the same length")
  if (is.null(names(testStat))) names(testStat) = probe.ids
  
  probe.ids<-probe.ids[(!is.na(gene.names[probe.ids]))&(gene.names[probe.ids]!="")]
  testStat<-testStat[probe.ids]
 # probe.ids <- as.factor(probe.ids)
  
  tSsp = split.default(testStat, gene.names[probe.ids])
  gr_prob_x<-function(x){
    if (length(x) >1){
      x <- x[grep(pattern="_x_at",names(x),invert=TRUE)]
      return(names(which.max(x)))
    } else{
      return(names(which.max(x)))
    }
  }
  out<-sapply(tSsp,gr_prob_x,USE.NAMES=T,simplify=T)

  #  return(out[c("YBX3","YIF1B","ZZZ3")]) #Testing subset containg "_x_at" probesets 
  return(out)
}