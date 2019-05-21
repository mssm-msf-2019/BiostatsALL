#' Function that runs the original Combat script where no model matrix of covariates is needed
#' @description It accpets an ExpressionSet and returns also an ExpressionSet. If you have a matrix, use RunCombatMat
#' @param eset: object of class ExpressionSet
#' @param g: factor defining the batches
#' @param covariates: names of the covariates to be found in pData(eset)
#' @examples 
#' abatch=ReadAffy()
#' abatch$Batch<-sapply(protocolData(abatch)$ScanDate,function(x){substr(x,1,10)}
#' eset<-rma(abatch)
#' RunCombat(eset,eset$Batch)
#' 
RunCombat<-function(eset,g,covariates=NULL){
  Z<-cbind(Array=sampleNames(eset),sampleNames(eset),Batch=factor(g))
  if (!is.null(covariates)) {Z<-cbind(Z,pData(eset)[,covariates])}
  write.table(file="BatchInfo.txt",Z,sep='\t',row.names=FALSE,col.names=TRUE)
  write.csv(file="ExprsBeforeBtachadj.csv",exprs(eset))
  Aeset<-ComBat("ExprsBeforeBtachadj.csv","BatchInfo.txt",type='csv',skip=1,write=F,prior.plots=T,filter=F)
  aux<-as.matrix(Aeset[,-1]);rownames(aux)<-as.character(Aeset[,1])
  esetA<-eset
  exprs(esetA)<-aux
  rm(aux)
  return(esetA)
}