#' function for pheatmap
#' @description function for pheatmap
#' @export
kmeans_pheatmap = function(mat, k = min(nrow(mat), 150), sd_limit = NA, ...){
  # Filter data
  if(!is.na(sd_limit)){
    s = apply(mat, 1, sd)
    mat = mat[s > sd_limit, ]
  }

  # Cluster data
  set.seed(1245678)
  km = kmeans(mat, k, iter.max = 100)
  mat2 = km$centers

  # Compose rownames
  t = table(km$cluster)
  rownames(mat2) = sprintf("cl%s_size_%d", names(t), t)

  # Draw heatmap
  pheatmap(mat2, ...)
}

