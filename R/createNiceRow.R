#' Function to pretty rows for createNiceTable()
#' @description Function to pretty rows for createNiceTable()

createNiceRow <- function(var, desc, byvar, data) {

  require(xtable)
  require(plyr)

  params <- as.list(desc[[var]])
  label <- params[[1]]
  if (is.null(params$grpFn)) {
    if (is.logical(data[,var]) || is.factor(data[,var])) {
      params$grpFn <- printLogicalPct
    } else {
      params$grpFn <- printMedianQ13
    }
  }

  if (is.null(params$summaryFn)) {
    if (is.logical(data[,var]) || is.factor(data[,var])) {
      params$summaryFn <- printChiSqPVal
    } else {
      params$summaryFn <- printKruskalPVal
    }
  }

  # Create a vector for the table row:
  # first element is the text
  # next elements are a summary statistic for each group
  # final element is a p-value for between-group significance
  rv <- c(label,
          laply(split(data[,var], data[,byvar]), params$grpFn),
          params$summaryFn(formula(paste(var, "~", byvar)), data))
  if (is.factor(data[,var])) {
    # how this works:
    # apply each factor level to an inline function
    # that takes data as an argument and adds a truth vector "$.tmpvec"
    # then do the same row-building as above
    # using the factor name, the group summarization function,
    # and then summarize using a contingency table based on $.tmpvec
    rv <- t(cbind(rv,
                  sapply(levels(data[,var]),
                         function(lvl, df) {
                           df$.tmpvec <- df[,var]==lvl
                           c(paste(" ",lvl),
                             laply(split(df$.tmpvec, df[,byvar]), params$grpFn),
                             params$summaryFn(formula(paste(".tmpvec ~", byvar)), df))
                         }, data)))
  }
  rbind(rv)
}
