#' function for pheatmap
#' @description function for pheatmap
#' @export

find_coordinates = function(n, gaps, m = 1:n){
  if(length(gaps) == 0){
    return(list(coord = unit(m / n, "npc"), size = unit(1 / n, "npc") ))
  }

  if(max(gaps) > n){
    stop("Gaps do not match with matrix size")
  }

  size = (1 / n) * (unit(1, "npc") - length(gaps) * unit("4", "bigpts"))

  gaps2 = apply(sapply(gaps, function(gap, x){x > gap}, m), 1, sum)
  coord = m * size + (gaps2 * unit("4", "bigpts"))

  return(list(coord = coord, size = size))
}

