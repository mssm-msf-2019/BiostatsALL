#' Function that given a correlation matrix and the pvalues, makes a heatmap representing the pairwise correlation.
#' @description Function that given a correlation matrix and the pvalues, makes a heatmap representing the pairwise correlation. It's used by CorrelationHeatmap_Samples functions. The correlation coeficient appears in the lower diagonal and the p-values  based * are in the upper diagonal.
#' @param db.cor: symetrix matrix indicating the  correlation coefficient for each column/row pair of variables
#' @param db.pcor: symetrix matrix indicating the p-value associated with the correlation coefficient for each column/row pair of variables
#' @param max.stars: pvalue that will be represented by one asterisc (*) . If max.stars=0.05 (default), ***=0.001,**=0.01, *=0.05
#' @param print.corrcoef: Logicalindicating if correlation coeficiets should be printed. Default to TRUE.
#' @param title: title for the plot
#' @param breaks: color breaks for heatmap, defualt seq(-1,1,0.25);
#' @param stars.size: size of the text overlay in the figure (correlations and p values)
#' @param pcuts: cuts for p-values
#' @param psymbs: stars for p-value
#' @param cors.round.places decimal places to round the correlations
#' @examples
#' DoCorrelationPlotwithStarts(db.cor, db.pcor, max.stars=0.05, print.corrcoef=TRUE)

DoCorrelationPlotwithStarts<-function (db.cor, db.pcor, max.stars = 0.05, print.corrcoef = TRUE,
                                       title = "Correlations", breaks = seq(-1, 1, 0.25), stars.size = 2,
                                       cors.size = 2, pcuts = c(0.001, 0.01, 0.05), psymbs = c("***",  "**", "*"),
                                       cors.round.places=2) {
  library(reshape)
  library(weights)
  print('Docorr')
  print(pcuts)
  print(psymbs)
  db.stars <- apply(db.pcor, 2, weights::starmaker, pcuts, psymbs)
  rownames(db.stars) <- colnames(db.stars)
  order.vars <- make.names(((rownames(db.cor))))
  db.stars <- db.stars[rownames(db.cor), rownames(db.cor)]
  db.stars[lower.tri(db.stars)] <- NA
  nbar <- data.frame(row = order.vars, as.data.frame(db.cor))
  nbap <- data.frame(row = order.vars, as.data.frame(db.stars))
  rownames(nbar) <- NULL
  rownames(nbap) <- NULL
  nbar.m <- rename(melt(nbar, id.vars = "row"), c(value = "cor"))
  nbap.m <- rename(melt(nbap, id.vars = "row"), c(value = "p.stars"))
  nbar.m$cor.breaks <- cut(nbar.m$cor, breaks = breaks, include.lowest = TRUE,
                           label = paste("(", breaks[-length(breaks)], ",", breaks[-1],
                                         ")", sep = ""))
  db.plot <- merge(nbar.m, nbap.m, by = c("row", "variable"),
                   sorted = T, all.x = T, all.y = T)
  db.plot <- mutate(db.plot, signif.stars = as.character(p.stars),
                    print.value = cor)
  db.plot[!is.na(db.plot$p.stars), "print.value"] <- NA
  db.plot[is.na(db.plot$p.stars), "signif.stars"] <- ""
  db.plot <- mutate(db.plot, signif.stars = as.character(signif.stars),
                    print.value = as.character(round(print.value, cors.round.places)))
  db.plot[with(db.plot, which(row == variable)), "print.value"] <- ""
  db.plot <- mutate(db.plot, row = factor(as.character(row),
                                          levels = rev(order.vars)), variable = factor(as.character(variable),
                                                                                       levels = order.vars))
  db.plot <- arrange(db.plot, row, variable)
  po.nopanel <- theme(panel.background = element_blank(), panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90,
                                                                                     hjust = 1, size = 12, colour = "black"), axis.text.y = element_text(size = 12,
                                                                                                                                                         colour = "black"))
  pa <- ggplot(db.plot, aes(row, variable)) +
    geom_tile(aes(fill = cor),  colour = "white") +
    scale_fill_gradient2(low = "blue",  high = "red", mid = "white", midpoint = median(breaks), limits = range(breaks)) +
    scale_x_discrete(breaks = levels(nbar.m$row),  labels = as.vector(sapply(levels(nbar.m$row), makenames_MSF))) +
    scale_y_discrete(breaks = levels(nbar.m$variable),
                     labels = as.vector(sapply(levels(nbar.m$variable), makenames_MSF))) +
    geom_text(aes(label = signif.stars),  colour = "black", size = stars.size, na.rm = TRUE) +
    labs(x =stackchar(paste0(psymbs,'(p<',pcuts,')'),sep=' '), y = "", title = title) +
    po.nopanel
  if ((print.corrcoef == TRUE) & (cors.size > 0))
    pa <- pa + geom_text(aes(label = print.value), colour = "black", size = cors.size, na.rm = TRUE)
  return(pa)
}
