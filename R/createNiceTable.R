#' A function to create descriptive table (Table 1)
#' @description a function that creates nicely looking descriptive table (Table 1)
#' @param desc - list describing variable and test. Available tests: printTtestPoolPVal, printTtestSattPVal, printWilcoxPVal, printAOVPVal,
#' printChiSqPVal, printFisherPVal, printKruskalPVal.
#' Available summaries: printMeanSD, printMedianQ13, printLogicalPct.
#' @param byvar - group variable (has to be a factor)
#' @param data - dataframe
#' @param dmeth the method to measure distance between matrices rows
#' @examples
#' table1 <- list(age = c("Age", grpFn=printMeanSD, summaryFn=printTtestSattPVal),
#' sex = c("Female",grpFn=printLogicalPct, summaryFn=printChiSqPVal))  # categorical variable have to be factors!
#' m <- createNiceTable(table1, "Treatment", df)


createNiceTable<-function (desc, byvar, data) {
  require(xtable)
  require(plyr)

  # Available tests: printTtestPoolPVal, printTtestSattPVal, printWilcoxPVal,
  # printAOVPVal, printChiSqPVal, printFisherPVal, printKruskalPVal.
  # Available summaries: printMeanSD, printMedianQ13, printLogicalPct.

  tbl <- table(data[, byvar])
  rv <- rbind(c("", t(names(tbl)), ""), c("", t(paste("n=", pInt(tbl))), ""),
              do.call(rbind, lapply(names(desc), createNiceRow, desc, byvar, data)))
  row.names(rv) <- NULL
  out<-rv[-1,]; colnames(out)<-c('variable',rv[1,c(2:(ncol(rv)-1))],'p-value')
  return(out)
}



# createNiceTable <- function(desc, byvar, data) {
#
#   require(xtable)
#   require(plyr)
#
#   # Available tests: printTtestPoolPVal, printTtestSattPVal, printWilcoxPVal,
#   # printAOVPVal, printChiSqPVal, printFisherPVal, printKruskalPVal.
#   # Available summaries: printMeanSD, printMedianQ13, printLogicalPct.
#
#   tbl <- table(data[,byvar])
#   # begin building the result matrix
#   # format for multigroup will be:
#   # desc | grpA  | grpB | grpC | summary
#   rv <- rbind(
#     c("", t(names(tbl)), ""),
#     c("", t(paste("n=",pInt(tbl))), ""),
#     do.call(rbind, lapply(names(desc), createNiceRow, desc, byvar, data)))
#   row.names(rv) <- NULL
#   rv
#
# }










