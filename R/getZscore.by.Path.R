#' Function given a T statistic and a collection of pathways, returns the Z score
#' @description Given a collection of pathwyas in SETs and the value of a T statistic in T, it returns the Z score for each patway
#' @param SETs matrix of number_of_genes by number_of_pathways. TRUE of FALSE if the gene belong to the pathway
#' @param y named vector of size number_of_genes. 
#' @examples
#' z=getZscore.by.Path(SETs,ebfit$t[,1])

getZscore.by.Path<-function(SETs,T){
  Zscore<-function(set,x){y<-x[set==TRUE]; t=sum(y)/sqrt(length(y)); p=2*pnorm(-abs(t)); return(c(t,p))}
  Res<-t(apply(SETs,2,Zscore,T))
  colnames(Res)<-c('Z','p')
  Res<-cbind(Res,FDR=p.adjust(Res[,2]))
  return(Res)
}
