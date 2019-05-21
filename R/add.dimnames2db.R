
#' a function to add colnames and rownames
#' @description Adding rownames and colnames; A essential step for pheatmap; Lots of plyr steps will lose dimnames
#' @param x: a matrix/dataframe
#' @param cn: colnames to be added
#' @param rn: rownames to be added

add.dimnames2db<-function(x,cn=colnames(x),rn=rownames(x)){
                     x<-as.data.frame(x)
                     colnames(x)<-cn
                     rownames(x)<-rn
                 return(x)}
