#' Plot data ellipses.
#'
#' @param level The confidence level at which to draw an ellipse (default is 0.75),
#'   or, if \code{type="euclid"}, the radius of the circle to be drawn.
#' @param type The type of ellipse.
#'   The default \code{"t"} assumes a multivariate t-distribution, and
#'   \code{"norm"} assumes a multivariate normal distribution.
#'   \code{"euclid"} draws a circle with the radius equal to \code{level},
#'   representing the euclidian distance from the center.
#'   This ellipse probably won't appear circular unless \code{coord_fixed()} is applied.
#' @param segments The number of segments to be used in drawing the ellipse.
#' @param na.rm If \code{FALSE} (the default), removes missing values with
#'   a warning.  If \code{TRUE} silently removes missing values.
#' @inheritParams stat_identity
#'
#' @details The method for calculating the ellipses has been modified from car::ellipse (Fox and Weisberg, 2011)
#'
#' @references
#' https://raw.githubusercontent.com/low-decarie/FAAV/master/r/stat-ellipse.R
#'
#' @export
#' @importFrom MASS cov.trob
#'
#' @examples
#' ggplot(faithful, aes(waiting, eruptions))+
#'   geom_point()+
#'   stat_ellipse()
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3))+
#'   geom_point()+
#'   stat_ellipse()
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3))+
#'   geom_point()+
#'   stat_ellipse(type = "norm", linetype = 2)+
#'   stat_ellipse(type = "t")
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3))+
#'   geom_point()+
#'   stat_ellipse(type = "norm", linetype = 2)+
#'   stat_ellipse(type = "euclid", level = 3)+
#'   coord_fixed()
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3))+
#'   stat_ellipse(geom = "polygon")

#
# stat_ellipse <- function(mapping=NULL, data=NULL, geom="path", position="identity", ...) {
#  # require(proto)
#   StatEllipse <- proto(ggplot2:::Stat,
# {
#   required_aes <- c("x", "y")
#   default_geom <- function(.) GeomPath
#   objname <- "ellipse"
#
#   calculate_groups <- function(., data, scales, ...){
#     .super$calculate_groups(., data, scales,...)
#   }
#   calculate <- function(., data, scales, level = 0.75, segments = 51,...){
#     dfn <- 2
#     dfd <- length(data$x) - 1
#     if (dfd < 3){
#       ellipse <- rbind(c(NA,NA))
#     } else {
#       require(MASS)
#       v <- cov.trob(cbind(data$x, data$y))
#       shape <- v$cov
#       center <- v$center
#       radius <- sqrt(dfn * qf(level, dfn, dfd))
#       angles <- (0:segments) * 2 * pi/segments
#       unit.circle <- cbind(cos(angles), sin(angles))
#       ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
#     }
#
#     ellipse <- as.data.frame(ellipse)
#     colnames(ellipse) <- c("x","y")
#     return(ellipse)
#   }
# }
#   )
#
#   StatEllipse$new(mapping=mapping, data=data, geom=geom, position=position, ...)
# }
