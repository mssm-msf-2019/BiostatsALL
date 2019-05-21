
#' A function to draw clustered heatmaps.
#'
#' A function to draw clustered heatmaps where one has better control over some graphical
#' parameters such as cell size, etc.
#'
#' The function also allows to aggregate the rows using kmeans clustering. This is
#' advisable if number of rows is so big that R cannot handle their hierarchical
#' clustering anymore, roughly more than 1000. Instead of showing all the rows
#' separately one can cluster the rows in advance and show only the cluster centers.
#' The number of clusters can be tuned with parameter kmeans_k.
#'
#' @param mat numeric matrix of the values to be plotted.
#' @param color vector of colors used in heatmap.
#' @param kmeans_k the number of kmeans clusters to make, if we want to agggregate the
#' rows before drawing heatmap. If NA then the rows are not aggregated.
#' @param breaks a sequence of numbers that covers the range of values in mat and is one
#' element longer than color vector. Used for mapping values to colors. Useful, if needed
#' to map certain values to certain colors, to certain values. If value is NA then the
#' breaks are calculated automatically.
#' @param border_color color of cell borders on heatmap, use NA if no border should be
#' drawn.
#' @param cellwidth individual cell width in points. If left as NA, then the values
#' depend on the size of plotting window.
#' @param cellheight individual cell height in points. If left as NA,
#' then the values depend on the size of plotting window.
#' @param scale character indicating if the values should be centered and scaled in
#' either the row direction or the column direction, or none. Corresponding values are
#' \code{"row"}, \code{"column"} and \code{"none"}
#' @param cluster_rows boolean values determining if rows should be clustered or \code{hclust} object,
#' @param cluster_cols boolean values determining if columns should be clustered or \code{hclust} object.
#' @param clustering_distance_rows distance measure used in clustering rows. Possible
#' values are \code{"correlation"} for Pearson correlation and all the distances
#' supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none
#' of the above it is assumed that a distance matrix is provided.
#' @param clustering_distance_cols distance measure used in clustering columns. Possible
#' values the same as for clustering_distance_rows.
#' @param clustering_method clustering method used. Accepts the same values as
#' \code{\link{hclust}}.
#' @param clustering_callback callback function to modify the clustering. Is
#' called with two parameters: original \code{hclust} object and the matrix
#' used for clustering. Must return a \code{hclust} object.
#' @param cutree_rows number of clusters the rows are divided into, based on the
#'  hierarchical clustering (using cutree), if rows are not clustered, the
#' argument is ignored
#' @param cutree_cols similar to \code{cutree_rows}, but for columns
#' @param treeheight_row the height of a tree for rows, if these are clustered.
#' Default value 50 points.
#' @param treeheight_col the height of a tree for columns, if these are clustered.
#' Default value 50 points.
#' @param legend logical to determine if legend should be drawn or not.
#' @param legend_breaks vector of breakpoints for the legend.
#' @param legend_labels vector of labels for the \code{legend_breaks}.
#' @param annotation_row data frame that specifies the annotations shown on left
#'  side of the heatmap. Each row defines the features for a specific row. The
#' rows in the data and in the annotation are matched using corresponding row
#'  names. Note that color schemes takes into account if variable is continuous
#'  or discrete.
#' @param annotation_col similar to annotation_row, but for columns.
#' @param annotation deprecated parameter that currently sets the annotation_col if it is missing
#' @param annotation_colors list for specifying annotation_row and
#' annotation_col track colors manually. It is  possible to define the colors
#' for only some of the features. Check examples for  details.
#' @param annotation_legend boolean value showing if the legend for annotation
#' tracks should be drawn.
#' @param annotation_names_row boolean value showing if the names for row annotation
#' tracks should be drawn.
#' @param annotation_names_col boolean value showing if the names for column annotation
#' tracks should be drawn.
#' @param drop_levels logical to determine if unused levels are also shown in
#' the legend
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param main the title of the plot
#' @param fontsize base fontsize for the plot
#' @param fontsize_row fontsize for rownames (Default: fontsize)
#' @param fontsize_col fontsize for colnames (Default: fontsize)
#' @param display_numbers logical determining if the numeric values are also printed to
#' the cells. If this is a matrix (with same dimensions as original matrix), the contents
#' of the matrix are shown instead of original values.
#' @param number_format format strings (C printf style) of the numbers shown in cells.
#' For example "\code{\%.2f}" shows 2 decimal places and "\code{\%.1e}" shows exponential
#' notation (see more in \code{\link{sprintf}}).
#' @param number_color color of the text
#' @param fontsize_number fontsize of the numbers displayed in cells
#' @param gaps_row vector of row indices that show shere to put gaps into
#'  heatmap. Used only if the rows are not clustered. See \code{cutree_row}
#'  to see how to introduce gaps to clustered rows.
#' @param gaps_col similar to gaps_row, but for columns.
#' @param labels_row custom labels for rows that are used instead of rownames.
#' @param labels_col similar to labels_row, but for columns.
#' @param filename file path where to save the picture. Filetype is decided by
#' the extension in the path. Currently following formats are supported: png, pdf, tiff,
#'  bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is
#' calculated so that the plot would fit there, unless specified otherwise.
#' @param width manual option for determining the output file width in inches.
#' @param height manual option for determining the output file height in inches.
#' @param silent do not draw the plot (useful when using the gtable output)
#' @param na_col specify the color of the NA cell in the matrix.
#' @param \dots graphical parameters for the text used in plot. Parameters passed to
#' \code{\link{grid.text}}, see \code{\link{gpar}}.
#'
#' @return
#' Invisibly a list of components
#' \itemize{
#'     \item \code{tree_row} the clustering of rows as \code{\link{hclust}} object
#'     \item \code{tree_col} the clustering of columns as \code{\link{hclust}} object
#'     \item \code{kmeans} the kmeans clustering of rows if parameter \code{kmeans_k} was
#' specified
#' }
#'
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#' # Create test matrix
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#'
#' # Draw heatmaps
#' pheatmap(test)
#' pheatmap(test, kmeans_k = 2)
#' pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
#' pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#' pheatmap(test, cluster_row = FALSE)
#' pheatmap(test, legend = FALSE)
#'
#' # Show text within cells
#' pheatmap(test, display_numbers = TRUE)
#' pheatmap(test, display_numbers = TRUE, number_format = "\%.1e")
#' pheatmap(test, display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test)))
#' pheatmap(test, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c("0",
#' "1e-4", "1e-3", "1e-2", "1e-1", "1"))
#'
#' # Fix cell sizes and save to file with correct size
#' pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap")
#' pheatmap(test, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")
#'
#' # Generate annotations for rows and columns
#' annotation_col = data.frame(
#'                     CellType = factor(rep(c("CT1", "CT2"), 5)),
#'                     Time = 1:5
#'                 )
#' rownames(annotation_col) = paste("Test", 1:10, sep = "")
#'
#' annotation_row = data.frame(
#'                     GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
#'                 )
#' rownames(annotation_row) = paste("Gene", 1:20, sep = "")
#'
#' # Display row and color annotations
#' pheatmap(test, annotation_col = annotation_col)
#' pheatmap(test, annotation_col = annotation_col, annotation_legend = FALSE)
#' pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row)
#'
#'
#' # Specify colors
#' ann_colors = list(
#'     Time = c("white", "firebrick"),
#'     CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
#'     GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
#' )
#'
#' pheatmap(test, annotation_col = annotation_col, annotation_colors = ann_colors, main = "Title")
#' pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row,
#'          annotation_colors = ann_colors)
#' pheatmap(test, annotation_col = annotation_col, annotation_colors = ann_colors[2])
#'
#' # Gaps in heatmaps
#' pheatmap(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14))
#' pheatmap(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14),
#'          cutree_col = 2)
#'
#' # Show custom strings as row/col names
#' labels_row = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
#' "", "", "Il10", "Il15", "Il1b")
#'
#' pheatmap(test, annotation_col = annotation_col, labels_row = labels_row)
#'
#' # Specifying clustering from distance matrix
#' drows = dist(test, method = "minkowski")
#' dcols = dist(t(test), method = "minkowski")
#' pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
#'
#' # Modify ordering of the clusters using clustering callback option
#' callback = function(hc, mat){
#'     sv = svd(t(mat))$v[,1]
#'     dend = reorder(as.dendrogram(hc), wts = sv)
#'     as.hclust(dend)
#' }
#'
#' pheatmap(test, clustering_callback = callback)
#'
#' \dontrun{
#' # Same using dendsort package
#' library(dendsort)
#'
#' callback = function(hc, ...){dendsort(hc)}
#' pheatmap(test, clustering_callback = callback)
#' }
#'
#' @export
pheatmap_msf = function(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), kmeans_k = NA, breaks = NA, border_color = "grey60", cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", clustering_callback = identity2, cutree_rows = NA, cutree_cols = NA,  treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows, 50, 0), treeheight_col = ifelse((class(cluster_cols) == "hclust") || cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA, legend_labels = NA, annotation_row = NA, annotation_col = NA, annotation = NA, annotation_colors = NA, annotation_legend = TRUE, annotation_names_row = TRUE, annotation_names_col = TRUE, drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8 * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL, labels_col = NULL, filename = NA, width = NA, height = NA, silent = FALSE, na_col = "#DDDDDD", ...){

  # Set labels
  if(is.null(labels_row)){
    labels_row = rownames(mat)
  }
  if(is.null(labels_col)){
    labels_col = colnames(mat)
  }

  # Preprocess matrix
  mat = as.matrix(mat)
  if(scale != "none"){
    mat = scale_mat(mat, scale)
    if(is.na2(breaks)){
      breaks = generate_breaks(mat, length(color), center = T)
    }
  }


  # Kmeans
  if(!is.na(kmeans_k)){
    # Cluster data
    km = kmeans(mat, kmeans_k, iter.max = 100)
    mat = km$centers

    # Compose rownames
    t = table(km$cluster)
    labels_row = sprintf("Cluster: %s Size: %d", names(t), t)
  }
  else{
    km = NA
  }

  # Format numbers to be displayed in cells
  if(is.matrix(display_numbers) | is.data.frame(display_numbers)){
    if(nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != ncol(mat)){
      stop("If display_numbers provided as matrix, its dimensions have to match with mat")
    }

    display_numbers = as.matrix(display_numbers)
    fmat = matrix(as.character(display_numbers), nrow = nrow(display_numbers), ncol = ncol(display_numbers))
    fmat_draw = TRUE

  }
  else{
    if(display_numbers){
      fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw = TRUE
    }
    else{
      fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw = FALSE
    }
  }

  # Do clustering
  if((class(cluster_rows) == "hclust") || cluster_rows){
    if(class(cluster_rows) == "hclust"){
      tree_row = cluster_rows
    } else {
      tree_row = cluster_mat(mat, distance = clustering_distance_rows, method = clustering_method)
      tree_row = clustering_callback(tree_row, mat)
    }
    mat = mat[tree_row$order, , drop = FALSE]
    fmat = fmat[tree_row$order, , drop = FALSE]
    labels_row = labels_row[tree_row$order]
    if(!is.na(cutree_rows)){
      gaps_row = find_gaps(tree_row, cutree_rows)
    }
    else{
      gaps_row = NULL
    }
  }
  else{
    tree_row = NA
    treeheight_row = 0
  }

  if((class(cluster_cols) == "hclust") || cluster_cols){
    if(class(cluster_cols) == "hclust"){
      tree_col = cluster_cols
    } else {
      tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, method = clustering_method)
      tree_col = clustering_callback(tree_col, t(mat))
    }
    mat = mat[, tree_col$order, drop = FALSE]
    fmat = fmat[, tree_col$order, drop = FALSE]
    labels_col = labels_col[tree_col$order]
    if(!is.na(cutree_cols)){
      gaps_col = find_gaps(tree_col, cutree_cols)
    }
    else{
      gaps_col = NULL
    }
  }
  else{
    tree_col = NA
    treeheight_col = 0
  }

  attr(fmat, "draw") = fmat_draw

  # Colors and scales
  if(!is.na2(legend_breaks) & !is.na2(legend_labels)){
    if(length(legend_breaks) != length(legend_labels)){
      stop("Lengths of legend_breaks and legend_labels must be the same")
    }
  }


  if(is.na2(breaks)){
    breaks = generate_breaks(as.vector(mat), length(color))
  }
  if (legend & is.na2(legend_breaks)) {
    legend = grid.pretty(range(as.vector(breaks)))
    names(legend) = legend
  }
  else if(legend & !is.na2(legend_breaks)){
    legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]

    if(!is.na2(legend_labels)){
      legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
      names(legend) = legend_labels
    }
    else{
      names(legend) = legend
    }
  }
  else {
    legend = NA
  }
  mat = scale_colours(mat, col = color, breaks = breaks, na_col = na_col)

  # Preparing annotations
  if(is.na2(annotation_col) & !is.na2(annotation)){
    annotation_col = annotation
  }
  # Select only the ones present in the matrix
  if(!is.na2(annotation_col)){
    annotation_col = annotation_col[colnames(mat), , drop = F]
  }

  if(!is.na2(annotation_row)){
    annotation_row = annotation_row[rownames(mat), , drop = F]
  }

  annotation = c(annotation_row, annotation_col)
  annotation = annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]
  if(length(annotation) != 0){
    annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
  }
  else{
    annotation_colors = NA
  }

  if(!show_rownames){
    labels_row = NULL
  }

  if(!show_colnames){
    labels_col = NULL
  }

  # Draw heatmap
  gt = heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend, annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, number_color = number_color, gaps_row = gaps_row, gaps_col = gaps_col, labels_row = labels_row, labels_col = labels_col, ...)

  if(is.na(filename) & !silent){
    grid.newpage()
    grid.draw(gt)
  }

  invisible(list(tree_row = tree_row, tree_col = tree_col, kmeans = km, gtable = gt))
}
