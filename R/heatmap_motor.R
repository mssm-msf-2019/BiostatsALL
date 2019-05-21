#' function for pheatmap
#' @description function for pheatmap
#' @export

heatmap_motor = function(matrix, border_color, cellwidth, cellheight, tree_col, tree_row, treeheight_col, treeheight_row, filename, width, height, breaks, color, legend, annotation_row, annotation_col, annotation_colors, annotation_legend, annotation_names_row, annotation_names_col, main, fontsize, fontsize_row, fontsize_col, fmat, fontsize_number, number_color, gaps_col, gaps_row, labels_row, labels_col, ...){
  # Set layout
  lo = lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, gaps_row = gaps_row, gaps_col = gaps_col,  ...)

  res = lo$gt
  mindim = lo$mindim

  if(!is.na(filename)){
    if(is.na(height)){
      height = convertHeight(gtable_height(res), "inches", valueOnly = T)
    }
    if(is.na(width)){
      width = convertWidth(gtable_width(res), "inches", valueOnly = T)
    }

    # Get file type
    r = regexpr("\\.[a-zA-Z]*$", filename)
    if(r == -1) stop("Improper filename")
    ending = substr(filename, r + 1, r + attr(r, "match.length"))

    f = switch(ending,
               pdf = function(x, ...) pdf(x, ...),
               png = function(x, ...) png(x, units = "in", res = 300, ...),
               jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
               jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
               tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...),
               bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
               stop("File type should be: pdf, png, bmp, jpg, tiff")
    )

    # print(sprintf("height:%f width:%f", height, width))

    # gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)

    f(filename, height = height, width = width)
    gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)
    grid.draw(gt)
    dev.off()

    return(gt)
  }

  # Omit border color if cell size is too small
  if(mindim < 3) border_color = NA

  # Draw title
  if(!is.na(main)){
    elem = draw_main(main, fontsize = 1.3 * fontsize, ...)
    res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", clip = "off")
  }

  # Draw tree for the columns
  if(!is.na2(tree_col) & treeheight_col != 0){
    elem = draw_dendrogram(tree_col, gaps_col, horizontal = T)
    res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
  }

  # Draw tree for the rows
  if(!is.na2(tree_row) & treeheight_row != 0){
    elem = draw_dendrogram(tree_row, gaps_row, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
  }

  # Draw matrix
  elem = draw_matrix(matrix, border_color, gaps_row, gaps_col, fmat, fontsize_number, number_color)
  res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", name = "matrix")

  # Draw colnames
  if(length(labels_col) != 0){
    pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, ...)
    elem = do.call(draw_colnames, pars)
    res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", name = "col_names")
  }

  # Draw rownames
  if(length(labels_row) != 0){
    pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, ...)
    elem = do.call(draw_rownames, pars)
    res = gtable_add_grob(res, elem, t = 4, l = 4, clip = "off", name = "row_names")
  }

  # Draw annotation tracks on cols
  if(!is.na2(annotation_col)){
    # Draw tracks
    converted_annotation = convert_annotations(annotation_col, annotation_colors)
    elem = draw_annotations(converted_annotation, border_color, gaps_col, fontsize, horizontal = T)
    res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", name = "col_annotation")

    # Draw names
    if(annotation_names_col){
      elem = draw_annotation_names(annotation_col, fontsize, horizontal = T)
      res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", name = "col_annotation_names")
    }
  }

  # Draw annotation tracks on rows
  if(!is.na2(annotation_row)){
    # Draw tracks
    converted_annotation = convert_annotations(annotation_row, annotation_colors)
    elem = draw_annotations(converted_annotation, border_color, gaps_row, fontsize, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", name = "row_annotation")

    # Draw names
    if(annotation_names_row){
      elem = draw_annotation_names(annotation_row, fontsize, horizontal = F)
      res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", name = "row_annotation_names")
    }
  }

  # Draw annotation legend
  annotation = c(annotation_col[length(annotation_col):1], annotation_row[length(annotation_row):1])
  annotation = annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]

  if(length(annotation) > 0 & annotation_legend){
    elem = draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)

    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, clip = "off", name = "annotation_legend")
  }

  # Draw legend
  if(!is.na2(legend)){
    elem = draw_legend(color, breaks, legend, fontsize = fontsize, ...)

    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, clip = "off", name = "legend")
  }

  return(res)
}
