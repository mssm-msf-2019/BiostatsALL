#' Function that produces the correlation plots and dendograms
#' @description Function that given a correlation matrix makes a correlation plot, dendograms, and phylo trees
#' @param method: pearson or spearman correlation
#' @param hmeth: method to calculate the distance
#' @param mat: matrix of values
savecorrplot<-function (mat, fname = "", dir.path = getwd(), method = "pearson",
          max.stars = 0.05, ht = 10, wd = 12, ht.d = 4.5, wd.d = 0.9 *
            wd, abs.cor=F,hmeth='average',...)
{
  l <- CorrelationHeatmap_Samples(mat, method = method, max.stars = 0.05,
                                  ...)
  fname <- gsub("__", "_", paste(fname, ifelse(abs.cor, "AbsCor",
                                               ""), method, hmeth, sep = "_"))
  ggsave(paste0(dir.path, "Correlations_", fname, ".pdf"),
         l$plot, height = ht, width = wd)
  pdf(paste0(dir.path, "Dendograms", fname, ".pdf"), height = ht.d,
      width = wd.d)
  plot(l$dendogram, hang = -1)
  dev.off()
}
