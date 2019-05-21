#' A function that reads GMT files 
#' @description A function that reads GMT files and returns a list where each element is a list of genes
#' @param fname: filename
#' @examples 
#' readGmt('KruegerSYMBOL2015.gmt')
readGmt<-function(fname){
  require(stringr)
  f <- readLines(fname)
  mc <- list()
  for (i in 1:length(f)) {
    dat <- unlist(strsplit(f[i], "\t", fixed = TRUE))
    dat.tmp<-dat[1]
    if(is.na(dat.tmp)){
    }else{
      m <- new("smc")
      m@reference <- str_trim(string=dat[1],side="right")
      if (is.na(dat[2])){
        m@desc <- ""
        
      } else{
        m@desc <- str_trim(string=dat[2],side="right")
      } 
      ids <- dat[3:length(dat)]
      m@ids <- ids[!(ids == "NA")]
      mc <- c(mc, list(m))
    }}
  names(mc) <- unlist(lapply(mc, function(x) paste(x@reference,x@desc)))
  return(mc)
}