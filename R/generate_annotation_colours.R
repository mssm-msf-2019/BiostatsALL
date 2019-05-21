#' function for pheatmap
#' @description function for pheatmap
#' @export

generate_annotation_colours = function(annotation, annotation_colors, drop){
  if(is.na2(annotation_colors)){
    annotation_colors = list()
  }
  count = 0
  for(i in 1:length(annotation)){
    annotation[[i]] = annotation[[i]][!is.na(annotation[[i]])]
    if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
      if (is.factor(annotation[[i]]) & !drop){
        count = count + length(levels(annotation[[i]]))
      }
      else{
        count = count + length(unique(annotation[[i]]))
      }
    }
  }

  factor_colors = dscale(factor(1:count), hue_pal(l = 75))

  oldseed = NULL
  if (exists(".Random.seed"))
    oldseed = get(".Random.seed", pos=.GlobalEnv)

  set.seed(3453)

  cont_counter = 2
  for(i in 1:length(annotation)){
    if(!(names(annotation)[i] %in% names(annotation_colors))){
      if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
        n = length(unique(annotation[[i]]))
        if (is.factor(annotation[[i]]) & !drop){
          n = length(levels(annotation[[i]]))
        }
        ind = sample(1:length(factor_colors), n)
        annotation_colors[[names(annotation)[i]]] = factor_colors[ind]
        l = levels(as.factor(annotation[[i]]))
        l = l[l %in% unique(annotation[[i]])]
        if (is.factor(annotation[[i]]) & !drop){
          l = levels(annotation[[i]])
        }

        names(annotation_colors[[names(annotation)[i]]]) = l
        factor_colors = factor_colors[-ind]
      }
      else{
        annotation_colors[[names(annotation)[i]]] = brewer_pal("seq", cont_counter)(5)[1:4]
        cont_counter = cont_counter + 1
      }
    }
  }

  if(!is.null(oldseed)){
    assign(".Random.seed", oldseed, pos=.GlobalEnv)
  }
  else{
    remove(.Random.seed, pos=.GlobalEnv)
  }

  return(annotation_colors)
}
