#' A function that collapse annotation of probes that mapp to more than one gene
#' @description This function given a matrix of genes collapses annotation of probes that mapp to more than one gene
#' @param gns matrix of genes
collapseduplicatesids<-function(gns){
  require(parallel)
  gnlst <- tapply(1:nrow(gns), gns$PROBEID, function(x) gns[x,])
  gnlst <- mclapply(gnlst, function(x) apply(x, 2, function(y) paste(uniquenona(y), collapse = " | ")))
  gns <- do.call("rbind", gnlst)
}