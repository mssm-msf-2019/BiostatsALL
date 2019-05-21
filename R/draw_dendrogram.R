#' function for pheatmap
#' @description function for pheatmap
#' @export

draw_dendrogram = function(hc, gaps, horizontal = T){
  h = hc$height / max(hc$height) / 1.05
  m = hc$merge
  o = hc$order
  n = length(o)

  m[m > 0] = n + m[m > 0]
  m[m < 0] = abs(m[m < 0])

  dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y")))
  dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)

  for(i in 1:nrow(m)){
    dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
    dist[n + i, 2] = h[i]
  }


  x = rep(NA, nrow(m) * 4)
  y = rep(NA, nrow(m) * 4)
  id = rep(1:nrow(m), rep(4, nrow(m)))

  for(i in 1:nrow(m)){
    c = draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
    k = (i - 1) * 4 + 1
    x[k : (k + 3)] = c$x
    y[k : (k + 3)] = c$y
  }

  x = find_coordinates(n, gaps, x * n)$coord
  y = unit(y, "npc")

  if(!horizontal){
    a = x
    x = unit(1, "npc") - y
    y = unit(1, "npc") - a
  }
  res = polylineGrob(x = x, y = y, id = id)

  return(res)
}
