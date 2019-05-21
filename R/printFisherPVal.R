#' Functions to get p-value from Fisher test
#' @description Functions to get p-value from Fisher test, used by createNiceTable()

printFisherPVal<-function (fmla, data) {
  tbl <- ftable(fmla, data)
  print(tbl)
  if (any(dim(tbl)!=c(2,2))) {printPVal(NA)}
  else{
    Xsq <- fisher.test(tbl)
    printPVal(Xsq$p.value)}
}
