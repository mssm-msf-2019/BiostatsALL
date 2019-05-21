#' Function to scale marix of evalues (used for pHeatmap)
#' @description Function to scale marix of evalues (used for pHeatmap)


scale2breaks<-function(mat,breaks,scale=T){
  mat<-as.matrix(mat)
  dnn<-dimnames(mat)
  if (scale) mat_scale<-t(apply(as.matrix(mat),1,scale)) else mat_scale=mat
  mat_scale[which(mat_scale>max(breaks))]<-max(breaks);
  mat_scale[which(mat_scale<min(breaks))]<-min(breaks);
  dimnames(mat_scale)<-dnn
  return(mat_scale)
}
