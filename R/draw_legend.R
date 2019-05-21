#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_legend = function(color, breaks, legend, ...){
  height = min(unit(1, "npc"), unit(150, "bigpts"))

  legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
  legend_pos = height * legend_pos + (unit(1, "npc") - height)

  breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
  breaks = height * breaks + (unit(1, "npc") - height)

  h = breaks[-1] - breaks[-length(breaks)]

  rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
  text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))

  res = grobTree(rect, text)

  return(res)
}
