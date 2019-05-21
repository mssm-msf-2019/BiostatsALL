#' Function for parsing a GMT gene set list
#' @description This function takes a GMT gene set list (from "readGmt") as input, and outputs a list of gene sets, with names being the geneset.names, and values the corresponding genesets. Note this function also works on a GSA .gmt file object.
#' @param gene_sets Name of geneset list that was created by "readGmt"
#' @param  trim_snam TRUE/FALSE if the geneset.names should be trimmed using "str_trim" from the "stringr" package. Defaults to FALSE.
#' @keywords gene sets, gmt
#' @examples 
#' jkr.symbol <- processGMT(jkr.list)
processGMT <- function(gene_sets,trim_snam=F){
  require(stringr)
  
  gene_sets$genesets <- gene_sets$genesets[!gene_sets$geneset.names==""]
  gene_sets$geneset.descriptions <- gene_sets$geneset.descriptions[!gene_sets$geneset.names==""]
  gene_sets$geneset.names <- gene_sets$geneset.names[!gene_sets$geneset.names==""]
  
  if(trim_snam==TRUE){
    gene_sets_out <- list()
    for(i in c(1:length(gene_sets$genesets))){
      gene_s1 <- as.character(unlist(gene_sets$genesets[i]))
      gene_s1 <- gene_s1[!gene_s1==""]
      gene_snam <- stringr::str_trim(as.character(unlist(gene_sets$geneset.name[i])))
      gene_snam_long <- as.character(unlist(gene_sets$geneset.descriptions[i]))
      gene_snam_long <- as.character(gene_snam_long)
      gene_s_len <- length(gene_s1)
      cat(i,gene_snam,"\n")
      gene_sets_out[gene_snam] <- list(gene_s1)
    } 
    
  }else{
    gene_sets_out <- list()
    for(i in c(1:length(gene_sets$genesets))){
      gene_s1 <- as.character(unlist(gene_sets$genesets[i]))
      gene_s1 <- gene_s1[!gene_s1==""]
      gene_snam <- as.character(unlist(gene_sets$geneset.name[i]))
      gene_snam_long <- as.character(unlist(gene_sets$geneset.descriptions[i]))
      gene_snam_long <- as.character(gene_snam_long)
      gene_s_len <- length(gene_s1)
      cat(i,gene_snam,"\n")
      gene_sets_out[gene_snam] <- list(gene_s1)
    } }
  
  return(gene_sets_out)
  
}

