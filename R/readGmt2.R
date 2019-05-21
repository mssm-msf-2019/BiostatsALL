#' Function for reading a GMT file of gene set collection
#' @description This function reads a gene set collection GMT file, and stores the output in a list of genesets, geneset.names, and  geneset.descriptions. This version is adapted from the GSA package by Robert Tibshirani, with some modifications with respect to the possible empty cells, making it more efficient.
#' @param fname File name of the .gmt file. The format of the file should be: 1.Column=Geneset Name ; 2.Column=Geneset Description ; 3.-N.Column=GENE SYMBOLS
#' @keywords gmt, gene sets
#' @examples 
#' jkr.list <- readGmt(fname = "./KruegerSYMBOL.gmt")

readGmt2<-function(fname){
{
  #Hi
  require(stringr)
  a = scan(fname, what = list("", ""), sep = "\t", quote = NULL, 
           fill = T, flush = T, multi.line = F)
  geneset.names = a[1][[1]]
  geneset.descriptions = a[2][[1]]
  dd = scan(fname, what = "", sep = "\t", quote = NULL)
  nn = length(geneset.names)
  n = length(dd)
  ox = rep(NA, nn)
  ii = 1
  for (i in 1:nn) {
    cat(i," ")
    while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] != geneset.descriptions[i])) {
      ii = ii + 1
    }
    ox[i] = ii
    ii = ii + 1
  }
  genesets = vector("list", nn)
  for (i in 1:(nn - 1)) {
    cat(i ," ", fill = T)
    i1 = ox[i] + 2
    i2 = ox[i + 1] - 1

    geneset.descriptions[i] = dd[ox[i] + 1]
    dd_tmp <- dd[i1:i2]
    dd_tmp <- dd_tmp[dd_tmp!=""]
    genesets[[i]] = dd_tmp
    
  }
  geneset.descriptions[nn] = dd[ox[nn] + 1]
  genesets[[nn]] = dd[(ox[nn] + 2):n]
  out = list(genesets = genesets, geneset.names = str_trim(string=geneset.names,side="right"), 
             geneset.descriptions = str_trim(string=geneset.descriptions,side="right"))
  class(out) = "smc"
  return(out)
}
}



