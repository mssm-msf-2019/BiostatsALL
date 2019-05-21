#' Functions to get p-value from  ANOVA
#' @description Functions to get p-value from ANOVA, used by createNiceTable()

printAOVPVal <- function(fmla, data) {
  fit <- aov(fmla, data=data)
  pv <- summary(fit)[[1]][["Pr(>F)"]][1]
  printPVal(pv)
}
