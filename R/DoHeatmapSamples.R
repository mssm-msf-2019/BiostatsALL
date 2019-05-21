#' A function that draws a heatmap based on a matrix
#' @description a function that is supposed to draw a heatmap based on a matrix of gene expressions
#' @param exprs is a matrix of expressions
#' @param corm is the metric for correlation evaluation
#' @param colfac a factor in the phenotype data that will add colors to the top according to different levels
#' @param dmeth the method to measure distance between matrices rows
#' @param hmeth the method for the hierarchical clustering
#' @param cex.genes is the proportion to which genes symbols should be expanded in the graphic
#' @param margins is the definition of vertical and horizontal margins to plot the heatmap
#' @param main is the title for the heatmap
#' @param ColD is a logical parameter that indicates if the columns should be clustered
#' @param scl is an indicator if data should be scaled by row, column or both
#' @param thebreaks enable color transition at specified limits
#' @param samplelabels keeps a vector of characters to be used in the heatmap columns
#' @param colcol is the color schema for the rows 
#' @param colscheme is the schema for colors in the body of the heatmap
#' @examples 
#' DoHeatmap(mat = exprs(eset),facsepr = mygroup1, colfac = mygroup2, symb = genesymbols, dmeth=cor.dist, hmeth='average', 
#' cex.genes=0.5, margins = c(5,15), main='', scl = row, breaks = NULL, cexCol = 1, colscheme = 'gray')

DoHeatmapSamples<-function(exprs, corm='pearson',colfac=NULL, dmeth=euc, hmeth='average', 
                           cex.genes=0.5, margins=c(5,15), main='',scl='none',thebreaks=NULL,
                           samplelabels=NULL,colcol=NULL,...){
  require(geneplotter)
  require(bioDist)
  mat<-cor(exprs,method=corm,use='p')
  
  #cc<-rainbow(length(levels(factor(colfac))),start=3/6,end=6/6)[as.numeric(factor(colfac))]
  if (is.null(colfac)) colcol=NULL
  else {
    if (is.null(colcol)) colcol <-getgray(factor(colfac)) 
    else colcol<-colcol[as.numeric(factor(colfac))]
  }
  print(colcol)
  if (is.null(thebreaks)) cols<-dChip.colors(50) else cols<-dChip.colors(length(thebreaks)-2);
  
  if (!is.null(samplelabels)&(length(samplelabels)==ncol(mat))) {colnames(mat)<-samplelabels; rownames(mat)<-samplelabels}
  
  if (!is.null(colfac)) 
    hl<-heatmap.2(mat, distfun=dmeth,hclustfun=function(c){hclust(c,method=hmeth)},
                  col=cols, trace='none', density.info='none', keysize=1, scale=scl, 
                  Colv=TRUE, RowSideColors= colcol,  ColSideColors= colcol, cexRow=cex.genes, 
                  cexCol=cex.genes, margins=margins, main=main, breaks=thebreaks,...)
  else   	
    hl<-heatmap.2(mat, distfun=dmeth,hclustfun=function(c){hclust(c,method=hmeth)},col=cols, trace='none',
                  density.info='none', keysize=1, scale=scl, Colv=TRUE, cexRow=cex.genes, cexCol=cex.genes,
                  margins=margins, main=main, breaks=thebreaks,...)
  return(invisible(hl))
}
