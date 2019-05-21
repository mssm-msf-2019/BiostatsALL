#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_annotation_names = function(annotations, fontsize, horizontal){
  n = ncol(annotations)

  x = unit(3, "bigpts")

  y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1
  y = unit(y, "bigpts")

  if(horizontal){
    res = textGrob(colnames(annotations), x = x, y = y, hjust = 0, gp = gpar(fontsize = fontsize, fontface = 2))
  }
  else{
    a = x
    x = unit(1, "npc") - y
    y = unit(1, "npc") - a

    res = textGrob(colnames(annotations), x = x, y = y, vjust = 0.5, hjust = 0, rot = 270, gp = gpar(fontsize = fontsize, fontface = 2))
  }

  return(res)
}
