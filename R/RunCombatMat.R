#' Function that runs the original Combat script where no model matrix of covariates is needed.
#' @description It accepts an a matrix and gives back a matrix
#' @param exprs: matrix of expression
#' @param g: factor defining the batches
#' #' @examples 
#' abatch=ReadAffy()
#' abatch$Batch<-sapply(protocolData(abatch)$ScanDate,function(x){substr(x,1,10)}
#' eset<-rma(abatch)
#' RunCombatMat(exprs(eset),eset$Batch)
#' 

RunCombatMat<-function(eset,g){
  Z<-cbind(Array=colnames(eset), colnames(eset),Batch=factor(g))
  write.table(file="BatchInfo.txt",Z,sep='\t',row.names=FALSE,col.names=TRUE)
  write.csv(file="ExprsBeforeBtachadj.csv",(eset))
  Aeset<-ComBat("ExprsBeforeBtachadj.csv","BatchInfo.txt",type='csv',skip=1,write=F,prior.plots=T,filter=F)
  aux<-as.matrix(Aeset[,-1]);rownames(aux)<-as.character(Aeset[,1])
  return(aux)
}