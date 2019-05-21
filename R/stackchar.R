#' A function that concatenates a vector of character using sperator sep
#' @description A function that concatenates a vector of character using sperator sep
#' @param x: vector of characters
#' @param sep: separator to be used
#' @examples
#' stackchar(LETTERS[1:10],",")

stackchar<-function(X,sep){
  if (length(X)==0) return("");
  XX<-NULL
  for (x in X){ XX<-paste(XX,x,sep=sep); 
                # Y<-substr(XX,2,nchar(XX))
  }
  return(XX)
}

