#' function for pheatmap
#' @description function for pheatmap
#' @export

lo = function(rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, treeheight_col, treeheight_row, legend, annotation_row, annotation_col, annotation_colors, annotation_legend, annotation_names_row, annotation_names_col, main, fontsize, fontsize_row, fontsize_col, gaps_row, gaps_col, ...){
  # Get height of colnames and length of rownames
  if(!is.null(coln[1]) | (!is.na2(annotation_row) & annotation_names_row)){
    if(!is.null(coln[1])){
      t = coln
    } else {
      t = ""
    }
    tw = strwidth(t, units = 'in', cex = fontsize_col / fontsize)
    if(annotation_names_row){
      t = c(t, colnames(annotation_row))
      tw = c(tw, strwidth(colnames(annotation_row), units = 'in'))
    }
    longest_coln = which.max(tw)
    gp = list(fontsize = ifelse(longest_coln <= length(coln), fontsize_col, fontsize), ...)
    coln_height = unit(1, "grobheight", textGrob(t[longest_coln], rot = 90, gp = do.call(gpar, gp))) + unit(10, "bigpts")
  }
  else{
    coln_height = unit(5, "bigpts")
  }

  if(!is.null(rown[1])){
    t = rown
    tw = strwidth(t, units = 'in', cex = fontsize_row / fontsize)
    if(annotation_names_col){
      t = c(t, colnames(annotation_col))
      tw = c(tw, strwidth(colnames(annotation_col), units = 'in'))
    }
    longest_rown = which.max(tw)
    gp = list(fontsize = ifelse(longest_rown <= length(rown), fontsize_row, fontsize), ...)
    rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], gp = do.call(gpar, gp))) + unit(10, "bigpts")
  }
  else{
    rown_width = unit(5, "bigpts")
  }

  gp = list(fontsize = fontsize, ...)
  # Legend position
  if(!is.na2(legend)){
    longest_break = which.max(nchar(names(legend)))
    longest_break = unit(1.1, "grobwidth", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
    title_length = unit(1.1, "grobwidth", textGrob("Scale", gp = gpar(fontface = "bold", ...)))
    legend_width = unit(12, "bigpts") + longest_break * 1.2
    legend_width = max(title_length, legend_width)
  }
  else{
    legend_width = unit(0, "bigpts")
  }

  # Set main title height
  if(is.na(main)){
    main_height = unit(0, "npc")
  }
  else{
    main_height = unit(1.5, "grobheight", textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))
  }

  # Column annotations
  textheight = unit(fontsize, "bigpts")

  if(!is.na2(annotation_col)){
    # Column annotation height
    annot_col_height = ncol(annotation_col) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")

    # Width of the correponding legend
    t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
    annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
    if(!annotation_legend){
      annot_col_legend_width = unit(0, "npc")
    }
  }
  else{
    annot_col_height = unit(0, "bigpts")
    annot_col_legend_width = unit(0, "bigpts")
  }

  # Row annotations
  if(!is.na2(annotation_row)){
    # Row annotation width
    annot_row_width = ncol(annotation_row) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")

    # Width of the correponding legend
    t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
    annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
    if(!annotation_legend){
      annot_row_legend_width = unit(0, "npc")
    }
  }
  else{
    annot_row_width = unit(0, "bigpts")
    annot_row_legend_width = unit(0, "bigpts")
  }

  annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)

  # Tree height
  treeheight_col = unit(treeheight_col, "bigpts") + unit(5, "bigpts")
  treeheight_row = unit(treeheight_row, "bigpts") + unit(5, "bigpts")

  # Set cell sizes
  if(is.na(cellwidth)){
    mat_width = unit(1, "npc") - rown_width - legend_width - treeheight_row - annot_row_width - annot_legend_width
  }
  else{
    mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * unit(4, "bigpts")
  }

  if(is.na(cellheight)){
    mat_height = unit(1, "npc") - main_height - coln_height - treeheight_col - annot_col_height
  }
  else{
    mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * unit(4, "bigpts")
  }

  # Produce gtable
  gt = gtable(widths = unit.c(treeheight_row, annot_row_width, mat_width, rown_width, legend_width, annot_legend_width), heights = unit.c(main_height, treeheight_col, annot_col_height, mat_height, coln_height), vp = viewport(gp = do.call(gpar, gp)))

  cw = convertWidth(mat_width - (length(gaps_col) * unit(4, "bigpts")), "bigpts", valueOnly = T) / ncol
  ch = convertHeight(mat_height - (length(gaps_row) * unit(4, "bigpts")), "bigpts", valueOnly = T) / nrow

  # Return minimal cell dimension in bigpts to decide if borders are drawn
  mindim = min(cw, ch)

  res = list(gt = gt, mindim = mindim)

  return(res)
}
