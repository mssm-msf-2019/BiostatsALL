#' function for pheatmap
#' @description function for pheatmap
#' @export

convert_annotations = function(annotation, annotation_colors){
  new = annotation
  for(i in 1:ncol(annotation)){
    a = annotation[, i]
    b = annotation_colors[[colnames(annotation)[i]]]
    if(is.character(a) | is.factor(a)){
      a = as.character(a)

      if(length(setdiff(setdiff(a, NA), names(b))) > 0){
        stop(sprintf("Factor levels on variable %s do not match with annotation_colors", colnames(annotation)[i]))
      }
      new[, i] = b[a]
    }
    else{
      a = cut(a, breaks = 100)
      new[, i] = colorRampPalette(b)(100)[a]
    }
  }
  return(as.matrix(new))
}
