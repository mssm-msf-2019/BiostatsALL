#' Function is called by auxBiplot().
#' @description this function is called by auxBiplot().
#' @param dir=getwd()  directory to save plot and table with resutls
#' @examples
#' auxBiplot(X,qVar='IL13',BiVar='Tissue',ylab=expression(log[2]('IL-13'),model=TRUE,dir=getwd(),logbase=NULL)

doBoxplot<-function(db,qLab='Quantivative Variable',BiLab='Categorical Variable',
                    db.avg=NULL,logbase=NULL){
  # db.avg<-mutate(db.avg,qVar=M)
  rg<-range(db$qVar)
  rg<-c(rg[1]-0.1*abs(diff(rg)),rg[2]+0.1*abs(diff(rg)))
  p<-ggplot(db,aes(y=qVar,x=BiVar))+geom_boxplot(outlier.size=0,outlier.color='gray') +
    geom_jitter(position = position_jitter(width = .15, height=0.01)) +
    geom_point(data=db.avg,aes(x=BiVar,y=qVar),color='darkgray',size=5) +
    geom_errorbar(data=db.avg,aes(ymin = qVar-SEM,ymax = qVar+SEM),
                  width=0.25, show_guide=FALSE, color='darkgray') +
    labs(y=qLab,x=BiLab) +
    coord_cartesian(ylim=rg)+ theme_bw()

  if (!is.null(logbase)){
    rg<-print(p)$panel$ranges[[1]]$y.minor_source
    lrg<-as.character(round(logbase^rg,2))
    p<-p+scale_y_continuous(breaks=rg, labels=lrg)
  }
  return(p)
}







