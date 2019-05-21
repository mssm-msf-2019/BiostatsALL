#' Internal function of dowordcloud. 
#' @description Transforms pvalues to colors
#' @param p p values determining colors
#' @param fch estimates (fold changes) red/blue colors will be used for positive/negatives values
#' @param leg_wcloud legend for colors
#' @param col_wcloud colors  
#' 
fromPtoCol<-function(p,fch=NULL,col_wcloud=c('mistyrose3','red3','darkred','greenyellow','green','darkgreen')){
  lgp<-(-log10(p))
  lgp.col<-rep('gray',length(lgp)); '#FFFFFF'
  # lgp.col[(fch>0)]<-'pink';lgp.col[(lgp>1)]<-'red3';lgp.col[(lgp>1.3)]<-'darkred'
  lgp.col[(fch>0)]<-col_wcloud[1]; lgp.col[(lgp>1)]<-col_wcloud[2]; lgp.col[(lgp>1.3)]<-col_wcloud[3]
  if (!is.null(fch)){
    lgp.col[(fch<0)]<-col_wcloud[4]; lgp.col[(lgp>1)&(fch<0)]<-col_wcloud[5]; lgp.col[(lgp>1.3)&(fch<0)]<-col_wcloud[6]
  }
  return( lgp.col)
}