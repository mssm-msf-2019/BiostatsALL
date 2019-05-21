#' Save output of a heatmap built from a decision matrix and coefficients from limma model
#' @description a function that preprocess output from limma model and draws a heatmap according to a set of parameters
#' @param Decision is a matrix indicating up(1), down(-1) regulation for each contrast is a gene list.
#' @param ebfit is a matrix of coefficients for each contrast
#' @param reset is the expression set used as input for the limma model
#' @param dmeth is the method for distance
#' @param hmeth is the agglomerative method used for the hierarchical clustering
#' @param cols is a vector defining how the columns in the hierarchical heatmap should be ordered
#' @param Group is a vector indicating the levels associated with columns in the expression sets
#' @param Patients is a vector indicating the patients IDs associated with columns in the expression sets
#' @examples 
#' DoCluster2(Decision = D, ebfit = fit, reset = eset, dmeth = "euclidean" , hmeth = "complete" , cols =NULL,Group=Group, Patients = IDs)


DoCluster2<-function(Decision,ebfit,reset,dmeth="euclidean",hmeth="complete",cols=NULL,Group,Patients){
  somediff<-(rowSums(Decision!=0)>0)
  
  symb<-unlist(mget(featureNames(reset),hgu133a2SYMBOL))
  somediff[is.na(symb)]<-FALSE
  
  mat<-exprs(reset)[somediff,]
  symb<-unlist(mget(featureNames(reset)[somediff],hgu133a2SYMBOL))
  name<-unlist(mget(featureNames(reset)[somediff],hgu133a2GENENAME))
  
  g<-(pData(reset)$Group)
  Patients<-(pData(reset)$Patient)
  mat2<-mat
  for (pat in levels(Patients)){
    mat2[,(Patients==pat)]<-mat[,(Patients==pat)] - rowMeans(mat[,(Patients==pat)])
  }
  
  mat<-mat2[,order(g)]
  g<-g[order(g)]
  cn<-Patients[order(g)]
  
  
  AveExp<-ebfit$Amean[somediff]
  o<-order(AveExp, decreasing=TRUE)
  rc<-heat.colors(nrow(mat))[o][o]
  
  desc<-paste(symb,name,sep=" - ")
  cc<-getgray(g)
  
  ho<-heatmap.2_MSF(mat, distfun=function(c){Dist(c,method=dmeth)},
                    hclustfun=function(c){hclust(c,method=hmeth)},
                    Colv=cols, col=greenred.colors(50),labRow=symb,
                    labCol=cn, RowSideColors=rc,ColSideColors=cc,
                    cexRow=0.4,cexCol=1,trace='none', density.info='none',
                    keysize=0.8, scale='row')
  
}