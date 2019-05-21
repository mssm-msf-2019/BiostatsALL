#' Function that plots the results of PVCA batch analysis
#' @description Function plots the results of PVCA batch analysis
#' @param pvcaObj: PVCA object from pvcaBatchAssess.MSF() function
#' @param fname: name of the file to be saved
#' @param ht: height of the plot
#' @param wd: width of the plot
#' @examples
#' pvcaBatch<-pvcaBatchAssess.MSF(eset,batch.factors= c('Batch','Tissue','WEEK','Treatment'), threshold=09,include.inter=c('WEEK:Treatment'))
#' plot.pvca(pvcaBatch,fname='./PVCA/PVCA_Batch')

plot.pvca<-function(pvcaObj, cex.percentage=1, fname=NULL,ht=4,wd=5,title=fname){

  require(stringr)
  require(ggplot2)

  title <- mgsub(c('/','.','PVCA'),rep('',3),fixed=TRUE,fname)
  labels <- gsub(":"," x ",pvcaObj$label)

  db <- data.frame(perc=t(pvcaObj$dat),labels= factor(labels,levels=labels))

  p <- ggplot(db, aes(x=labels,y=perc)) +
    geom_bar(stat='identity',fill='blue') +
    theme_bw() +
    scale_y_continuous(limits=c(0,1.1)) +
    geom_text(aes(x=labels,y=perc+0.05,label=paste0(as.character(round(100*perc,1)),'%'))) +
    labs(x='Effects',y="Weighted average proportion variance",title=paste("PVCA ",title)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 8))

  if (!is.null(fname)) {ggsave(file=paste0(fname,'.pdf'), height=ht,width=wd,plot=p)}

  return(p)

}
