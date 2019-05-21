#' A function to solve ambiguities in duplicated symbols for the same probeid
#' @description a function that eliminates duplicates ids and return two symbols, the first found and a synonim
#' @param gns is a dataframe with annotation about a geneset including probeset IDS, GO , PATH and CHRLOC
#' @examples 
#' treatduplicatesids(gns=mygeneset)

treatduplicatesids<-function(gns){
  require(parallel)
  gnlst <- tapply(1:nrow(gns), gns$PROBEID, function(x) gns[x,])
  gnlst <- mclapply(gnlst, function(x){apply(x, 2, function(y) collapsefcn(setdiff(y, firstelementfcn(y))))})
  gns1 <- do.call("rbind", gnlst)
  
  gnlst2 <- tapply(1:nrow(gns), gns$PROBEID, function(x) gns[x,setdiff(colnames(gns),c('PATH','GO','CHRLOC'))])
  gnlst2 <- mclapply(gnlst2, function(x) apply(x, 2, firstelementfcn))
  gns2 <- do.call("rbind", gnlst2)
  
  return(list(First=gns2, Synon=gns1))
}