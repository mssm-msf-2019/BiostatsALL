#' A function that draws a heatmap based on a matrix
#' @description a function that is supposed to draw a heatmap based on a matrix of gene expressions
#' @param mat a matrix with genes in the rows and samples in the columns
#' @param facsepr a factor in the phenotype data that separates columns in the heatmap
#' @param colfac a factor in the phenotype data that will add colors to the top according to different levels
#' @param dendogram indicates if dendogram should  be added to rows, columns or both
#' @param symb a character vector that indicates gene symbols for each row
#' @param dmeth the method to measure distance between matrices rows
#' @param hmeth the method for the hierarchical clustering
#' @param cex.genes is the proportion to which genes symbols should be expanded in the graphic
#' @param cexCol is the proportion to which columns labels should be expanded in the graphic
#' @param margins is the definition of vertical and horizontal margins to plot the heatmap
#' @param ColD is a logical parameter that indicates if the columns should be clustered
#' @param scl is an indicator if data should be scaled by row, column or both
#' @param breaks enable color transition at specified limits
#' @param main is the title for the heatmap
#' @param colscheme is the schema for colors in the body of the heatmap
#' @examples 
#' DoHeatmap.3_MSF(mat = exprs(eset),facsepr = mygroup1, colfac = mygroup2, symb = genesymbols, dmeth=cor.dist, hmeth='average', 
#' cex.genes=0.5, margins = c(5,15), main='', scl = row, breaks = NULL, cexCol = 1, colscheme = 'gray')

DoHeatmap.3_MSF<-function(mat, facsepr=NULL, colfac, dendogram='row',symb=NULL, dmeth=function(x){cor.dist(x,abs=F)}, hmeth='average',
                          cex.genes=0.5,cexCol=1, margins=c(5,15), ColD=FALSE, main='',scl='row', breaks=c(seq(-2,2,0.1)), colscheme='gray', labRow="", ...){
  require(geneplotter)
  require(bioDist)
  require(gplots)
  require(GMD)
  colChoice.2_MSF <- function(colpat="blue",colfacs){
    require(plyr)
    n=length(levels(colfacs))
    nrs <-topo.colors(n)
    print(levels(colfacs))
    colfacs_out<-mapvalues(colfacs, from = levels(colfacs), to = nrs)
    #    colfacs_out <- sapply(colfacs_out, as.character)
    return(colfacs_out)
  }
  AveExp<-apply(mat,1,mean)
  o<-order(AveExp, decreasing=TRUE)
  
  nfac<-ncol(colfac)
  rc<-heat.colors(nrow(mat))[o][o]
  cl<-colnames(mat)
  if (is.null(breaks)) cols<-dChip.colors(50) else cols<-dChip.colors(length(breaks)-2); 
  if (!ColD) colsep<-which(facsepr[-1]!= facsepr[-length(facsepr)]) else colsep<-NULL
  print(colsep)
  if (!is.null(symb)) symb<-symb[rownames(mat)]
  
  if (nfac==1){
    if (class(colfac)!='factor') colfac<-factor(colfac)
    if (length(colscheme)==length(levels(colfac)))
      cc<-colscheme[as.numeric(colfac)]
    else cc<-rainbow(length(levels(colfac)),start=3/6,end=6/6)[as.numeric(colfac)]
    if (colscheme=='gray') cc<-getgray(colfac)
    hl<-heatmap.2(mat, distfun=dmeth,hclustfun=function(c){hclust(c,method=hmeth)},col=cols, trace='none', density.info='none', keysize=1, scale=scl, Colv=ColD, RowSideColors=rc, ColSideColors=cc,
                  cexRow=cex.genes, cexCol= cexCol, labRow=symb, labCol=cl, margins=margins, main=main, breaks=breaks, colsep=colsep,...)
    
  } else if (nfac==2) {
    
    index=1; 
    rownames(colfac)<-colnames(mat)
    colfac.asfac<-colfac
    for (i in colnames(colfac)) {
      if (class(colfac[,i])=='numeric') {
        qs<-round(quantile(colfac[,i],na.rm=TRUE),1);
        colfac[,i]<-cut(colfac[,i],qs)                                        
      }
      if (class(colfac.asfac[,i])!='factor') colfac.asfac[,i]<-factor(colfac.asfac[,i])
      if (index==1) {cc1<-getgray(colfac.asfac[,i])}
      else if (index==2) cc2<-colChoice.2_MSF(colfacs=colfac.asfac[,i])
      index=index+1
    }
    clab<-cbind(sapply(cc1, as.character),sapply(cc2, as.character))
    colnames(clab)<-colnames(colfac);  rownames(clab) <- rownames(colfac)
    clab<-clab[colnames(mat),]
    
    legend_c<-c(levels(cc1),levels(cc2))
    legend_i<-c(levels(colfac[,1]),levels(colfac[,2]))
    print('bbb')
    print(legend_c)
    print(legend_i)
    
    hl<-BiostatsALL::heatmap.3(x=mat, density.info="none", Colv=ColD, dendogram=dendogram, scale=scl, col=cols, trace="none", 
                               ColSideColors=clab, margins=margins, main=main, labCol=cl, labRow=symb, distfun=dmeth, 
                               hclustfun=function(c){hclust(c,method=hmeth)},cexRow=cex.genes, cexCol=cexCol, breaks=breaks,
                               symbreaks=TRUE, colsep=colsep, ...)
    
    #,RowSideColors=as.matrix(t(rc)))
    legend("topright",legend=legend_i,fill=legend_c, border=T, bty="n", y.intersp = 0.7, cex=0.7, x.intersp = 0.25,ncol=3)
  } else {
    stop("Error in DoHeatmap: nfac MUST be either 1 or 2!")
  }#' @param main is the title for the heatmap
  
  return(invisible(hl))
}
