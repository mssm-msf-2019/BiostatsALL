#' Function that given a data base of continuos variables, clusters it by variable and represents the associations as a heatmap-colored correlation matrix.
#' @description Function that given a data base of continuos variables, clusters it by variable and represents the associations as a heatmap-colored correlation matrix in the order established by undupervised clustering analysis. The correlation coeficient appears in the lower diagonal and the p-values  based * are in the upper diagonal.
#' @param db: data base of continuos variables; columns are variables to be clustered.
#' @param method: correlation method (spearman, perason etc, default to spearman)
#' @param abs.cor: should the distance be 1-cor (FALSE,default) or 1-|cor| (TRUE). In teh later case, variables negative correlated will be as near as positive correlated
#' @param hmeth: aglomeration method to do the k-means clustering. default to average
#' @param max.stars: pvalue that will be represented by one asterisc (*) . If max.stars=0.05 (default), ***=0.001,**=0.01, *=0.05
#' @param breaks: color breaks for heatmap, defualt seq(-1,1,0.25);
#' @param print.corrcoef: Parameters to control the heatmap aestetics. Logicalin dicating if correlation coeficiets should be printed. Default to TRUE.
#' @param title: Parameters to control the heatmap aestetics. title for the plot
#' @param stars.size: Parameters to control the heatmap aestetics. size of the text overlay in the figure (correlations and p values)
#' @examples
#' # in this examples, I want to cluster samples
#' CorrelationHeatmap_Samples(exprs(reset)[1:10,],print.corrcoef=F) #supress text
#' # in this examples, I want to cluster genes
#' CorrelationHeatmap_Samples(t(exprs(reset)[1:10,]), max.stars=0.1,print.corrcoef=F) #***=0.01,**=0.05, *=0.1

CorrelationHeatmap_Samples<-function (db, method = "spearman", abs.cor = FALSE, hmeth = "average",
                                      max.stars = 0.05, breaks = seq(-1, 1, 0.25), stars.size = 2,
                                      cors.size = 2, print.corrcoef = TRUE, pcuts = c(0.001, 0.01,  0.05),
                                      psymbs = c("***", "**", "*"),...) {
  library(geneplotter)
  library(gplots)
  print('CorrelationHeatmap_Samples')
  print(pcuts)
  print(psymbs)
  db.cor <- cor(db, use = "p", m = method)
  if (method == "spearman")
    db.pcor <- cor.prob.spearman(db)
  else db.pcor <- cor.prob.pearson(db)
  dd <- dis.cors(t(db), method = method, abs = abs.cor)
  h <- hclust(dd, hmeth)
  title <- paste("Distance: ", ifelse(abs.cor, "|", ""), method,
                 ifelse(abs.cor, "|", ""), " / ", hmeth, sep = "")
  nv <- length(h$order)
  db.pcor <- db.pcor[h$order[nv:1], h$order[nv:1]]
  db.cor <- db.cor[h$order[nv:1], h$order[nv:1]]
  cormat <- db.cor
  cormat[lower.tri(db.cor)] <- db.pcor[lower.tri(db.cor)]
  p <- DoCorrelationPlotwithStarts(db.cor = db.cor, db.pcor = db.pcor,
                                   title = title, breaks = breaks, stars.size = stars.size,
                                   cors.size = cors.size, print.corrcoef = print.corrcoef,
                                   pcuts=pcuts, psymbs=psymbs,...)
  return(list(plot = p, cormat = cormat, dendogram = h))
}


