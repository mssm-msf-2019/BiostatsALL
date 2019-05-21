#' function for pheatmap
#' @description function for pheatmap
#' @export

find_gaps = function(tree, cutree_n){
  v = cutree(tree, cutree_n)[tree$order]
  gaps = which((v[-1] - v[-length(v)]) != 0)

}
