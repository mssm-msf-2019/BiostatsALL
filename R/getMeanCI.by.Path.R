#' Function given a value in y and a collection of pathways
#' @description Given a collection of pathwyas in SETs and a value of an attribute in y, it summarizes the mean, CI and p values foe y for each patway
#' @param SETs matrix of number_of_genes by number_of_pathways. TRUE of FALSE if the gene belong to the pathway
#' @param y named vector of size number_of_genes. 

getMeanCI.by.Path<-function(SETs,y){
  n<-colSums(SETs)
  res<-apply(SETs[names(y),n>2],2,function(x,y){summary(lm(y~0+x,data=data.frame(x=x,y=y)))$coef[2,]},as.vector(y))
  res<-cbind(n=n[colnames(res)],t(res))
  res<-cbind(res,uiw=qt(0.975, res[,3]))
  res<-cbind(res,FDRs=p.adjust(res[,5],'BH'))
  return(res)
}
