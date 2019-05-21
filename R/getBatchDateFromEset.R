#' getBatchDateFromEset
#' @description Returns batch dates from an eset
#' @param eset The eset from which you wish to aquire batch dates from.
#' @export

getBatchDateFromEset<-function(eset){
  Batch.Date<-sapply(pData(protocolData(eset)[sampleNames(eset),])$ScanDate,
                     function(x){substr(x,1,10)}, simplify=TRUE,USE.NAMES=FALSE)
  return(Batch.Date)
}
