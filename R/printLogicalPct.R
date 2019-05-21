#' Functions to get frequency and column percent
#' @description Functions to get frequency and column percent, used by createNiceTable()

printLogicalPct <- function(vec) {
  if (is.factor(vec)) {
    ""
  } else {
    sprintf("%s (%.1f%%)",pInt(sum(vec, na.rm=TRUE)),sum(vec,na.rm=TRUE)/length(vec)*100)
  }
}
