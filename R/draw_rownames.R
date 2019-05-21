#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_rownames = function(rown, gaps, ...){
  coord = find_coordinates(length(rown), gaps)
  y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)

  res = textGrob(rown, x = unit(3, "bigpts"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))

  return(res)
}
