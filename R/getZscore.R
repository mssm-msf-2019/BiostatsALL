#' Internal function to calculate Z scores for gene expression data 
#' @description It calculates the Z scores after normilizing data
#' @param mat: expression matrix
#' @examples
#' m1up<-exprs(reset)[intersect(featureNames(reset),JKList[["MADAD (LSvsNL) Up "]]),]
#' m1down<-exprs(reset)[intersect(featureNames(reset),JKList[["MADAD (LSvsNL) Down "]]),]
#' MADAD_LSvsNL<-getzscore(rbind(m1up,-1*m1down))


getZscore<-function(mat){
  zscores_mat<-t(apply(mat,1,function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)}))
  zscore<-apply(zscores_mat,2,function(x){x<-x[!is.na(x)]; sum(x/sqrt(length(x)))});
  return(zscore)
}
