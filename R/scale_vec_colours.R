#' Function to replace values outside the breaks with max/min color (pheatmap)
#' @description Function to replace values outside the breaks with max/min color (pheatmap)

scale_vec_colours = function (x, col = rainbow(10), breaks = NA, na_col){
  x[which(x>max(breaks))]<-max(breaks);
  x[which(x<min(breaks))]<-min(breaks);
  x<-as.numeric(cut(x, breaks = breaks, include.lowest = T))
  return(col[x])

}
