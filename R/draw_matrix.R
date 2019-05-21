#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_matrix = function(matrix, border_color, gaps_rows, gaps_cols, fmat, fontsize_number, number_color){
  n = nrow(matrix)
  m = ncol(matrix)

  coord_x = find_coordinates(m, gaps_cols)
  coord_y = find_coordinates(n, gaps_rows)

  x = coord_x$coord - 0.5 * coord_x$size
  y = unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)

  coord = expand.grid(y = y, x = x)

  res = gList()

  res[["rect"]] = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = coord_y$size, gp = gpar(fill = matrix, col = border_color))

  if(attr(fmat, "draw")){
    res[["text"]] = textGrob(x = coord$x, y = coord$y, label = fmat, gp = gpar(col = number_color, fontsize = fontsize_number))
  }

  res = gTree(children = res)

  return(res)
}
