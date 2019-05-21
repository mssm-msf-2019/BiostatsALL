#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_colnames = function(coln, gaps, ...){
  coord = find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size

  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))

  return(res)
}
