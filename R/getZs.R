#' A function that calculates Z-scores and p-values
#' @description This function given a set of values, returns Z-scores, p-values and p-values adjusted 
#' @param SETs numeric matrix

getZs<-function(SETs,T){
  Zscore<-function(set,x){y<-x[set==TRUE]; t=sum(y)/sqrt(length(y)); p=2*pnorm(-abs(t)); return(c(t,p))}
  Res<-t(apply(SETs,2,Zscore,T))
  colnames(Res)<-c('Z','p')
  Res<-cbind(Res,FDR=p.adjust(Res[,2]))
  return(Res)
}
