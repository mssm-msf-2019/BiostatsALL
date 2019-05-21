#' Given a vector of characters a b c it returns a single character of the form a |b |c 
#' @description Given a vector of characters a b c it returns a single character of the form a |b |c  . Internal function to use in annotation. 
#' @param y vector of characters
#' @examples
#' collapsefcn()
collapsefcn <- function(y=c("a","b","c")){
	y <- uniquenona(y) 
	y <- ifelse(length(y)==0,'',y[1])
	return(y)
	}
