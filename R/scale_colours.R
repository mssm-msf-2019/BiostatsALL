#' function for pheatmap
#' @description function for pheatmap
#' @export

scale_colours = function(mat, col = rainbow(10), breaks = NA, na_col){
  mat = as.matrix(mat)
  return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks, na_col = na_col), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}
