#' A function that will plot roc curves
#' @description a function that will plot roc curves
#' @param data is a data frame with the columns named pred and obs
#' @param cols is the color to plot the roc curve
#' @param output is the name to save the graph
#' @param wd is the width of the figure
#' @param ht is the height of the figure
#' @param graftitle is the title for the graphics
#' @examples
#' plotROC(x,y)

plotROCbyGroups<-function(data,
                          cols=c('lightblue','green','cyan','darkblue','purple'),
                          output ='roc.pdf',
                          group=group,wd=8,ht=6,graftitle=''){

  # loading ROCR
  require(ROCR)
  require(pROC)


  fit<-by(data[,c('obs','pred')],data[,group],function(d)pROC::roc(d$obs,d$pred))

  getSpec<-function(L){L$specificities}
  getSens<-function(L){L$sensitivities}
  getGroupLength<-function(L){length(L$specificities)}

   spec<-unlist(lapply(fit,getSpec))
   sens<-unlist(lapply(fit,getSens))
   n<-unlist(lapply(fit,getGroupLength))


  rocdata<-data.frame(FPR = 1-spec, TPR = sens,Group=rep(names(n),n))




  diag = data.frame(x = seq(0, 1, by = 0.01 ), y = seq(0, 1, by = 0.01))

 p <- ggplot(data = rocdata, aes(x = FPR, y = TPR,colour=Group)) +
             scale_colour_manual(values=cols,name=group)+
             geom_point(size=0.3) +
             geom_line(aes(colour = Group),width=0.5) +
             geom_line(data =diag, aes(x = x, y = y), color = "red",width=0.5)



  f <- p  +
  theme_bw()+
  theme(axis.text = element_text(size = 16),
        title = element_text(size = 18)) +
  labs(x = "False Positive Rate", y = "True Positive Rate", title = graftitle)

ggsave(f, filename = output, width = wd, height = ht)

return(f)
}
