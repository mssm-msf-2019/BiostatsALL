#' Plots a boxplot of qVar over different categories of BiVar. If model=TRUE, it overlays the estimated mean and SEM.
#' @description It selects the N most dissimilar models among those in caret suite by maximizing the Jaccard dissimilarity between sets of models.
#' @param dat: database
#' @param qVar: quantitative variable
#' @param BiVar: categorical variable
#' @param qLab: label for quantitative variable
#' @param BiLab: label for categorical variable
#' @param model=TRUE,  overlays mean and SEM
#' @param logbase=NULL: if the data was transfromed using logb, this willl allows to plot results in the naturale scale, logbase is the base used in the log transformation
#' @param fdata=''  Suffix Name for the plot; It would be composed as paste(dir,'/Boxplot',fdata,'_',qVar,'By',BiVar,sep='')
#' @param dir=getwd()  directory to save plot and table with resutls
#' @examples
#' auxBiplot(X,qVar='IL13',BiVar='Tissue',ylab=expression(log[2]('IL-13'),model=TRUE,dir=getwd(),logbase=NULL)

auxBiplot<-function (dat, qVar, BiVar, qLab = qVar, BiLab = BiVar, model = TRUE,
                     logbase = NULL, fdata = "", dir = getwd()) {
  db <- data.frame(qVar = dat[, qVar], BiVar = dat[, BiVar])
  db <- db[(!is.na(db$qVar)) & (!is.na(db$BiVar)), ]
  db.avg <- ddply(db, "BiVar", function(df) return(c(qVar = mean(df$qVar),
                                                     SD = sd(df$qVar), n = length(df$qVar), SEM = sd(df$qVar)/sqrt(length(df$qVar)),
                                                     MIN = min(df$qVar))))
  p <- doBoxplot(db=db, qLab=qLab, BiLab=BiLab, db.avg=db.avg, logbase = logbase)
  if (model) {
    f <- lm(qVar ~ BiVar, db)
    s <- summary(f)
    write.csv(file = paste(dir, "Summary", fdata, "_", qVar, "By", BiVar, ".csv", sep = ""),
              x = rbind(summary(lm(qVar ~ 0 + BiVar, db))$coef, summary(f)$coef[2, ]))

    p <- p + annotate("text", x = 0.75, y = 1.1 * max(db$qVar, na.rm = T),
                      label = prettyp(anova(f)['BiVar', 5]), parse = TRUE)
  }
  ggsave(file = paste(dir, "/Boxplot", fdata, "_", qVar, "By",
                      BiVar, ".pdf", sep = ""), plot = p, height = 5, width = 5)
  return(p)
}
