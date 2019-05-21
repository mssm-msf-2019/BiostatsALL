#' A function that changes the names of a dataframe columns
#' @description This function changes the names of a dataframe columns adding a suffix to the current names (uses "_" to separate coulmn name and suffix)
#' @param db dataframe 
#' @param suffixb dataframe 

renamecols<-function(db,suffix){colnames(db)<-paste(colnames(db),suffix,sep='_'); return(db)}