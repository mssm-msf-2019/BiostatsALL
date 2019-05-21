#' Function that runs the original Combat script where no model matrix of covariates is needed. It calls a version created by MSF to return several parts of the internal modeling
#' @description See ComBat_MSF. It accepts an ExpressionSet and returns also an ExpressionSet. 
#' @param eset: object of class ExpressionSet
#' @param g: factor defining the batches
#' @param covariates: names of the covariates to be found in pData(eset)
#' @examples 
#' abatch=ReadAffy()
#' abatch$Batch<-sapply(protocolData(abatch)$ScanDate,function(x){substr(x,1,10)}
#' eset<-rma(abatch)
#' RunCombat_MSF(eset,eset$Batch)
#' 



RunCombat_MSF<-function(eset,g,covariates=NULL){
  Z<-cbind(Array=sampleNames(eset),sampleNames(eset),Batch=factor(g))
  if (!is.null(covariates)) {Z<-cbind(Z,pData(eset)[,covariates])}
  write.table(file="BatchInfo.txt",Z,sep='\t',row.names=FALSE,col.names=TRUE)
  write.csv(file="ExprsBeforeBtachadj.csv",exprs(eset))
  output<-ComBat_MSF("ExprsBeforeBtachadj.csv","BatchInfo.csv",type='csv',skip=1,write=F,prior.plots=T,filter=F)
  aux<-as.matrix(output$exprs); dimnames(aux)<-dimnames(exprs(eset)); 
  ad.eset<-eset; exprs(ad.eset)<-aux
  return(ad.eset)
}