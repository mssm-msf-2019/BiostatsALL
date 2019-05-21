#' Function that prints a dendogram using ggplot. Was modified from ggdendrogram 
#' @description Function that prints a dendogram using ggplot. Was modified from ggdendrogram . Refer to ggdendro package http://www.icesi.edu.co/CRAN/web/packages/ggdendro/vignettes/ggdendro.pdf . Full description of how to print dendograms here: http://rstudio-pubs-static.s3.amazonaws.com/1876_df0bf890dd54461f98719b461d987c3d.htmlThe correlation coeficient appears in the lower diagonal and the p-values  based * are in the upper diagonal. 
#' @param Either a dendro object or an object that can be coerced to class dendro using the dendro_data function, i.e. objects of class dendrogram, hclust or tree
#' @param segments If TRUE, show line segments
#' @param labels if TRUE, shows segment labels
#' @param leaf_labels if TRUE, shows leaf labels
#' @param rotate if TRUE, rotates plot by 90 degrees
#' @param theme_dendro if TRUE, applies a blank theme to plot (see theme_dendro)
#' @param ... other parameters passed to geom_text
#' @examples
#' ggDendrogram_MSF(DATA, rotate=FALSE)



ggDendrogram_MSF<-function (data, segments = TRUE, labels = TRUE, leaf_labels = TRUE, 
                            rotate = FALSE, theme_dendro = TRUE, size=12, ...) 
{
	require(ggplot2)
	require(ggdendro)
	
	
	# library(ggplot2)
	# hc <- hclust(dist(USArrests), "ave")
	# ### demonstrate plotting directly from object class hclust
	# p <- ggDendrogram_MSF(hc, rotate=FALSE)
	# print(p)
  dataClass <- if (inherits(data, "dendro")) 
    data$class
  else class(data)
  angle <- if (dataClass %in% c("dendrogram", "hclust")) {
    ifelse(rotate, 0, 90)
  }
  else {
    ifelse(rotate, 90, 0)
  }
  hjust <- if (dataClass %in% c("dendrogram", "hclust")) {
    ifelse(rotate, 0, 1)
  }
  else {
    0.5
  }
  #if (!is.dendro(data)) 
  data <- dendro_data(data)
  p <- ggplot()
  if (all(segments, !is.null(data$segments))) {
    p <- p + geom_segment(data = segment(data), aes_string(x = "x", 
                                                           y = "y", xend = "xend", yend = "yend"))
  }
  if (all(leaf_labels, !is.null(data$leaf_labels))) {
    p <- p + geom_text(data = leaf_label(data), aes_string(x = "x", 
                                                           y = "y", label = "label"), hjust = hjust, angle = angle, ...)
  }
  if (rotate) {
    p <- p + scale_x_discrete(labels = sapply(data$labels$label,makenames_MSF))
  }
  else {
    p <- p + scale_x_discrete(labels = sapply(data$labels$label,makenames_MSF))
  }
  if (rotate) {
    p <- p + coord_flip()
    p <- p + scale_y_continuous()
  }
  else {
    p <- p + scale_y_continuous()
  }
  if (theme_dendro) 
    p <- p + theme_dendro()
  p <- p + theme(axis.text.x = element_text(angle = angle, hjust = 1, size=size, colour='black'))
  p <- p + theme(axis.text.y = element_text(angle = angle, hjust = 1, size=size, colour='black'))
  p
}