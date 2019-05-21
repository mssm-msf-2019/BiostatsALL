#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_connection = function(x1, x2, y1, y2, y){
  res = list(
    x = c(x1, x1, x2, x2),
    y = c(y1, y, y, y2)
  )

  return(res)
}
