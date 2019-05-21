#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_annotation_legend = function(annotation, annotation_colors, border_color, ...){
  y = unit(1, "npc")
  text_height = unit(1, "grobheight", textGrob("FGH", gp = gpar(...)))

  res = gList()
  for(i in names(annotation)){
    res[[i]] = textGrob(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = "bold", ...))

    y = y - 1.5 * text_height
    if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
      n = length(annotation_colors[[i]])
      yy = y - (1:n - 1) * 2 * text_height

      res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = 2 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = annotation_colors[[i]]))
      res[[paste(i, "t")]] = textGrob(names(annotation_colors[[i]]), x = text_height * 2.4, y = yy - text_height, hjust = 0, vjust = 0.5, gp = gpar(...))

      y = y - n * 2 * text_height

    }
    else{
      yy = y - 8 * text_height + seq(0, 1, 0.25)[-1] * 8 * text_height
      h = 8 * text_height * 0.25

      res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = 2 * text_height, gp = gpar(col = NA, fill = colorRampPalette(annotation_colors[[i]])(4)))
      res[[paste(i, "r2")]] = rectGrob(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = 8 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = NA))

      txt = rev(range(grid.pretty(range(annotation[[i]], na.rm = TRUE))))
      yy = y - c(1, 7) * text_height
      res[[paste(i, "t")]]  = textGrob(txt, x = text_height * 2.4, y = yy, hjust = 0, vjust = 0.5, gp = gpar(...))
      y = y - 8 * text_height
    }
    y = y - 1.5 * text_height
  }

  res = gTree(children = res)

  return(res)
}
