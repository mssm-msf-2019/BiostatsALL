#' function for pheatmap
#' @description function for pheatmap
#' @export

is.na2 = function(x){
  if(is.list(x) | length(x) > 1){
    return(FALSE)
  }
  if(length(x) == 0){
    return(TRUE)
  }

  return(is.na(x))
}
