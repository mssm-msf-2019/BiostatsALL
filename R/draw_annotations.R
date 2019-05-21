#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_annotations = function(converted_annotations, border_color, gaps, fontsize, horizontal){
  n = ncol(converted_annotations)
  m = nrow(converted_annotations)

  coord_x = find_coordinates(m, gaps)

  x = coord_x$coord - 0.5 * coord_x$size

  # y = cumsum(rep(fontsize, n)) - 4 + cumsum(rep(2, n))
  y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1
  y = unit(y, "bigpts")

  if(horizontal){
    coord = expand.grid(x = x, y = y)
    res = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = unit(fontsize, "bigpts"), gp = gpar(fill = converted_annotations, col = border_color))
  }
  else{
    a = x
    x = unit(1, "npc") - y
    y = unit(1, "npc") - a

    coord = expand.grid(y = y, x = x)
    res = rectGrob(x = coord$x, y = coord$y, width = unit(fontsize, "bigpts"), height = coord_x$size, gp = gpar(fill = converted_annotations, col = border_color))
  }

  return(res)
}
