#' Create output from a list of pathways
#' @description create an output from a list of pathways
#' @param STATs is the stats used
#' @param SETs are the pathways we are interested in
#' @param fname is the file name
#' @examples 
#' dogetPathByCoef(STATs,SETs,fname)

dogetPathByCoef<-function(STATs,SETs,fname){
  for (i in c(1:ncol(STATs)))    
  {write.csv(file=paste(fname,colnames(STATs)[i], '.csv',sep=''),getPath(SETs,STATs[,i]))}
}


