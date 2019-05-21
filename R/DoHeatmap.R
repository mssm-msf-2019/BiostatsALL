#' A function that draws a heatmap based on a matrix
#' @description a function that is supposed to draw a heatmap based on a matrix of gene expressions
#' @param mat a matrix with genes in the rows and samples in the columns
#' @param facsepr a factor in the phenotype data that separates columns in the heatmap
#' @param colfac a factor in the phenotype data that will add colors to the top according to different levels
#' @param symb a character vector that indicates gene symbols for each row
#' @param dmeth the method to measure distance between matrices rows
#' @param hmeth the method for the hierarchical clustering
#' @param cex.genes is the proportion to which genes symbols should be expanded in the graphic
#' @param margins is the definition of vertical and horizontal margins to plot the heatmap
#' @param ColD is a logical parameter that indicates if the columns should be clustered
#' @param scl is an indicator if data should be scaled by row, column or both
#' @param breaks enable color transition at specified limits
#' @param cexCol is the proportion to which columns labels should be expanded in the graphic
#' @param main is the title for the heatmap
#' @param colscheme is the schema for colors in the body of the heatmap
#' @examples 
#' DoHeatmap(mat = exprs(eset),facsepr = mygroup1, colfac = mygroup2, symb = genesymbols, dmeth=cor.dist, hmeth='average', 
#' cex.genes=0.5, margins = c(5,15), main='', scl = row, breaks = NULL, cexCol = 1, colscheme = 'gray')

DoHeatmap<-function(mat, facsepr=NULL, colfac, symb=NULL, dmeth=cor.dist, hmeth='average', cex.genes=0.5,margins=c(5,15), 
                    ColD=FALSE,main='',scl='row', breaks=NULL,cexCol=1,colscheme='gray',...){
  library(geneplotter)
  library(bioDist)
  library(gplots)
  AveExp<-apply(mat,1,mean)
  o<-order(AveExp, decreasing=TRUE)
  if (class(colfac)!='factor') colfac<-factor(colfac)
  if (length(colscheme)==length(levels(colfac)))	
    cc<-colscheme[as.numeric(colfac)]
  else cc<-rainbow(length(levels(colfac)),start=3/6,end=6/6)[as.numeric(colfac)]
  if (colscheme=='gray') cc<-getgray(colfac)
  rc<-heat.colors(nrow(mat))[o][o]
  cl<-colnames(mat)
  
  if (is.null(breaks)) cols<-dChip.colors(50) else cols<-dChip.colors(length(breaks)-2);	if (!ColD) colsep<-which(facsepr[-1]!= facsepr[-length(facsepr)]) else colsep<-NULL
  print(colsep)
  if (!is.null(symb)) symb<-symb[rownames(mat)] 
  hl<-heatmap.2(mat, distfun=dmeth,hclustfun=function(c){hclust(c,method=hmeth)},col=cols, trace='none', density.info='none', keysize=1, scale=scl, Colv=ColD, RowSideColors=rc,  ColSideColors=cc, cexRow=cex.genes, cexCol= cexCol, labRow=symb, labCol=cl, margins=margins, main=main, breaks=breaks, colsep=colsep,...)
  return(invisible(hl))
}


