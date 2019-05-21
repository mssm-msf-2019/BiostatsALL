#' Function from John Fox that is called by scatter3d.Fox.
#' @description It returns a matriz with 3 columns  to be used to draw confidence ellipsoid in 3d.
#' @param center center of the ellipsoid
#' @param radius radius of the ellipsoid
#' @param shape shape of the ellipsoid
#' @examples 
#' ellipsoid.Fox()
#' 
ellipsoid.Fox <- function(center=c(0, 0, 0), radius=1, shape=diag(3),segments=51) {
  angles <- (0:segments)*2*pi/segments
  ecoord2 <- function(p) {
    c(cos(p[1])*sin(p[2]), sin(p[1])*sin(p[2]), cos(p[2])) }
  unit.sphere <- t(apply(expand.grid(angles, angles), 1, ecoord2))
  t(center + radius * t(unit.sphere %*% chol(shape))) 
}
