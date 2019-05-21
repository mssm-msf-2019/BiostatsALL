#' Creates a wordcloud representin Fold changes and color coded to significance
#' @description the features in vars will have size proportional to the column Estimate in Dmat and color coded to column pvalues. Positive estimates (up-regulated) will be in shades of red depending on significance, and negative estimates (down-regulated) will be in shades of blue)
#' @param vars variables to add to wordcloud
#' @param Dmat matrix where rows are features to appear as words and columns are elements of the comparison to be used. column Estimate and pvalue will be used for size/colors of words
#' @param leg_wcloud legend for colors
#' @param col_wcloud colors  
#' @param addFCH it contains numbers that will be ploted as to give an idea of which fold changes represents teh size.
#' @param logEst should be False unless to want to size it by log2 scale
#' @keywords wordcloud ; plotting
#' @examples 
#' dowordcloud() 

dowordcloud <- function(vars, Dmat	, title, addFCH=NULL, logEst=FALSE, leg_wcloud=c('Up-(p>0.1)','Up-(p<0.1)','Up-(p<0.05)','Down-(p>0.1)','Down-(p<0.1)','Down-(p<0.05)'), col_wcloud=c('lightpink2','red','darkred','lightblue','deepskyblue2','blue3')){
  
  freq<-abs(Dmat[vars,'Estimate'])
  
  if (logEst==TRUE) {freq<-log2(freq)}
  col.pval<-fromPtoCol(abs(Dmat[vars,intersect(colnames(Dmat),c("Pr(>|t|)","p","P"))]), Dmat[vars,'Estimate'], col_wcloud)
  words<-sapply(gsub(".Total","", vars),makenames_MSF)
  
  if (!is.null(addFCH)){
    words<-c(as.character(addFCH),words); 
    freq<-c(log2(addFCH),freq); 
    col.pval<-c(rep('black',length(addFCH)),col.pval)
  }
  
  wordcloud(words, freq=freq, colors=col.pval, ordered.colors=TRUE, min.freq=0, max.freq=10000, random.order=FALSE,use.r.layout=TRUE)
  legend("bottom",c(1:6),leg_wcloud[c(1,4,2,5,3,6)], text.col=col_wcloud[c(1,4,2,5,3,6)], border='white', box.col='white',ncol=3)
  mtext(title,3,line=1,adj=.5)
}