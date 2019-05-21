#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_main = function(text, ...){
  res = textGrob(text, gp = gpar(fontface = "bold", ...))

  return(res)
}
