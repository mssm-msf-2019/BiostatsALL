#' Function that draws an scatterplot with three dimensions with confidence ellipsoid. From John Fox. This is an older version. Please use ellipsoid3d in package car 
#' @description  The scatter3d function uses the rgl package to draw 3D scatterplots with various regression surfaces. The function identify3d allows you to label points interactively with the mouse: Press the right mouse button (on a two-button mouse) or the centre button (on a three-button mouse), drag a rectangle around the points to be identified, and release the button. Repeat this procedure for each point or set of “nearby” points to be identified. To exit from point-identification mode, click the right (or centre) button in an empty region of the plot.
#' @param x,y,x vectors to be plotted
#' @param ellipsoid TRUE if a confidence ellipsoid is to be draw. Default to FALSE
#' @param xlab, ylab, zlab axis labels.
#' @param axis.scales if TRUE, label the values of the ends of the axes. Note: For identify3d to work properly, the value of this argument must be the same as in scatter3d.
#' @param revolutions number of full revolutions of the display.
#' @param bg.col background colour; one of "white", "black".
#' @param axis.col colours for axes; if axis.scales is FALSE, then the second colour is used for all three axes.
#' @param surface.col vector of colours for regression planes, used in the order specified by fit; for multi-group plots, the colours are used for the regression surfaces and points in the several groups.
#' @param surface.alpha transparency of regression surfaces, from 0.0 (fully transparent) to 1.0 (opaque); default is 0.5.
#' @param neg.res.col, pos.res.col colours for lines representing negative and positive residuals.
#' @param square.col colour to use to plot squared residuals.
#' @param point.col colour of points.
#' @param text.col colour of axis labels.
#' @param grid.col colour of grid lines on the regression surface(s).
#' @param fogtype type of fog effect; one of "exp2", "linear", "exp", "none".
#' @param residuals plot residuals if TRUE; if residuals="squares", then the squared residuals are shown as squares (using code adapted from Richard Heiberger). Residuals are available only when there is one surface plotted.
#' @param surface plot surface(s) (TRUE or FALSE).
#' @param fill fill the plotted surface(s) with colour (TRUE or FALSE).
#' @param grid plot grid lines on the regression surface(s) (TRUE or FALSE).
#' @param grid.lines number of lines (default, 26) forming the grid, in each of the x and z directions.
#' @param df.smooth degrees of freedom for the two-dimensional smooth regression surface; if NULL (the default), the gam function will select the degrees of freedom for a smoothing spline by generalized cross-validation; if a positive number, a fixed regression spline will be fit with the specified degrees of freedom.
#' @param df.additive degrees of freedom for each explanatory variable in an additive regression; if NULL (the default), the gam function will select degrees of freedom for the smoothing splines by generalized cross-validation; if a positive number or a vector of two positive numbers, fixed regression splines will be fit with the specified degrees of freedom for each term.
#' @param sphere.size general size of spheres representing points; the actual size is dependent on the number of observations.
#' @param radius relative radii of the spheres representing the points. This is normally a vector of the same length as the variables giving the coordinates of the points, and for the formula method, that must be the case or the argument may be omitted, in which case spheres are the same size; for the default method, the default for the argument, 1, produces spheres all of the same size. The radii are scaled so that their median is 1.
#' @param threshold if the actual size of the spheres is less than the threshold, points are plotted instead.
#' @param speed relative speed of revolution of the plot.
#' @param fov field of view (in degrees); controls degree of perspective.
#' @param fit one or more of "linear", "quadratic", "smooth", "additive"; to display fitted surface(s); partial matching is supported -- e.g., c("lin", "quad").
#' @param groups if NULL (the default), no groups are defined; if a factor, a different surface or set of surfaces is plotted for each level of the factor; in this event, the colours in surface.col are used successively for the points, surfaces, and residuals corresponding to each level of the factor.
#' @param parallel when plotting surfaces by groups, should the surfaces be constrained to be parallel? A logical value, with default TRUE.
#' @param ellipsoid plot concentration ellipsoid(s) (TRUE or FALSE).
#' @param level expected proportion of bivariate-normal observations included in the concentration ellipsoid(s); default is 0.5.
#' @param model.summary print summary or summaries of the model(s) fit (TRUE or FALSE). scatter3d rescales the three variables internally to fit in the unit cube; this rescaling will affect regression coefficients.
#' @examples 
#' x=rnorm(50); y=3*x+4+rnorm(50,0,.1); z=x+y+rnorm(50,0,.1)
#' scatter3d.Fox(x,y,x)
#' 
#' 
scatter3d.Fox <- function(x, y, z, xlab=deparse(substitute(x)),
                      ylab=deparse(substitute(y)),  zlab=deparse(substitute(z)), revolutions=0,
                      bg.col=c("white", "black"), axis.col=if (bg.col == "white") "black" else "white",
                      surface.col=c("blue", "green", "orange", "magenta","cyan", "red", "yellow", "gray"),
                      neg.res.col="red", pos.res.col="green",point.col="yellow",
                      text.col=axis.col, grid.col=if (bg.col == "white") "black" else "gray",
                      fogtype=c("exp2", "linear", "exp", "none"),
                      residuals=(length(fit) == 1), surface=TRUE, grid=TRUE, grid.lines=26,
                      df.smooth=NULL, df.additive=NULL,
                      sphere.size=1, threshold=0.01, speed=1, fov=60,
                      fit="linear", groups=NULL, parallel=TRUE, ellipsoid=FALSE, level=0.5, model.summary=FALSE){
  
  ###Function from John Fox
  
  require(rgl)
  require(mgcv)
  summaries <- list()
  if ((!is.null(groups)) && (nlevels(groups) > length(surface.col))) 
    stop(sprintf(gettextRcmdr("Number of groups (%d) exceeds number of
                              colors (%d)."),
                 nlevels(groups), length(surface.col)))
  if ((!is.null(groups)) && (!is.factor(groups)))
    stop(gettextRcmdr("groups variable must be a factor."))
  bg.col <- match.arg(bg.col)
  fogtype <- match.arg(fogtype)
  if ((length(fit) > 1) && residuals && surface)
    stop(gettextRcmdr("cannot plot both multiple surfaces and
                      residuals"))
  xlab  # cause these arguments to be evaluated
  ylab
  zlab
  rgl.clear()
  rgl.viewpoint(fov=fov)
  rgl.bg(col=bg.col, fogtype=fogtype)
  valid <- if (is.null(groups)) complete.cases(x, y, z)
  else complete.cases(x, y, z, groups)
  x <- x[valid]
  y <- y[valid]
  z <- z[valid]
  if (!is.null(groups)) groups <- groups[valid]
  x <- (x - min(x))/(max(x) - min(x))
  y <- (y - min(y))/(max(y) - min(y))
  z <- (z - min(z))/(max(z) - min(z))
  size <- sphere.size*((100/length(x))^(1/3))*0.015
  if (is.null(groups)){
    if (size > threshold) rgl.spheres(x, y, z, color=point.col,
                                      radius=size)
    else rgl.points(x, y, z, color=point.col)
  }
  else {
    if (size > threshold) rgl.spheres(x, y, z,
                                      color=surface.col[as.numeric(groups)], radius=size)
    else rgl.points(x, y, z, color=surface.col[as.numeric(groups)])
  }
  rgl.lines(c(0,1), c(0,0), c(0,0), color=axis.col)
  rgl.lines(c(0,0), c(0,1), c(0,0), color=axis.col)
  rgl.lines(c(0,0), c(0,0), c(0,1), color=axis.col)
  rgl.texts(1, 0, 0, xlab, adj=1, color=text.col)
  rgl.texts(0, 1, 0, ylab, adj=1, color=text.col)
  rgl.texts(0, 0, 1, zlab, adj=1, color=text.col)
  if (ellipsoid) {
    dfn <- 3
    if (is.null(groups)){
      dfd <- length(x) - 1
      radius <- sqrt(dfn * qf(level, dfn, dfd))
      ellips <- ellipsoid.Fox(center=c(mean(x), mean(y), mean(z)),
                              shape=cov(cbind(x,y,z)), radius=radius)
      quads3d(ellips[,1], ellips[,2], ellips[,3], front="lines",
              back="lines", alpha=.5, 
              lit=FALSE, col=surface.col[1])
    }
    else{
      levs <- levels(groups)
      for (j in 1:length(levs)){
        group <- levs[j]
        select.obs <- groups == group
        xx <- x[select.obs]
        yy <- y[select.obs]
        zz <- z[select.obs]
        dfd <- length(xx) - 1
        radius <- sqrt(dfn * qf(level, dfn, dfd))
        ellips <- ellipsoid(center=c(mean(xx), mean(yy), mean(zz)),
                            shape=cov(cbind(xx,yy,zz)), radius=radius)
        quads3d(ellips[,1], ellips[,2], ellips[,3], front="lines",
                back="lines", alpha=.5, 
                lit=FALSE, col=surface.col[j])
      }
    }
  }               
  if (surface){
    vals <- seq(0, 1, length=grid.lines)
    dat <- expand.grid(x=vals, z=vals)
    for (i in 1:length(fit)){
      f <- match.arg(fit[i], c("linear", "quadratic", "smooth",
                               "additive"))
      if (is.null(groups)){
        mod <- switch(f,
                      linear = lm(y ~ x + z),
                      quadratic = lm(y ~ (x + z)^2 + I(x^2) + I(z^2)),
                      smooth = if (is.null(df.smooth)) gam(y ~ s(x, z))
                      else gam(y ~ s(x, z, fx=TRUE, k=df.smooth)),
                      additive = if (is.null(df.additive)) gam(y ~ s(x) +s(z))
                      else gam(y ~ s(x, fx=TRUE, k=df.additive[1]+1) +
                                 s(z, fx=TRUE, k=(rev(df.additive+1)[1]+1)))
        )
        if (model.summary) summaries[[f]] <- summary(mod)
        yhat <- matrix(predict(mod, newdata=dat), grid.lines,grid.lines)
        rgl.surface(vals, vals, yhat, color=surface.col[i],alpha=0.5, lit=FALSE)
        if(grid) rgl.surface(vals, vals, yhat, color=grid.col,alpha=0.5, lit=FALSE, front="lines", back="lines")
        if (residuals){
          n <- length(y)
          fitted <- fitted(mod)
          colors <- ifelse(residuals(mod) > 0, pos.res.col,
                           neg.res.col)
          rgl.lines(as.vector(rbind(x,x)),
                    as.vector(rbind(y,fitted)), as.vector(rbind(z,z)),
                    color=as.vector(rbind(colors,colors)))
        }
      }
      else{
        if (parallel){
          mod <- switch(f,
                        linear = lm(y ~ x + z + groups),
                        quadratic = lm(y ~ (x + z)^2 + I(x^2) + I(z^2) +
                                         groups),
                        smooth = if (is.null(df.smooth)) gam(y ~ s(x, z) +
                                                               groups)
                        else gam(y ~ s(x, z, fx=TRUE, k=df.smooth) +
                                   groups),
                        additive = if (is.null(df.additive)) gam(y ~ s(x) +
                                                                   s(z) + groups)
                        else gam(y ~ s(x, fx=TRUE, k=df.additive[1]+1) +
                                   s(z, fx=TRUE, k=(rev(df.additive+1)[1]+1)) +
                                   groups)
          )
          if (model.summary) summaries[[f]] <- summary(mod)
          levs <- levels(groups)
          for (j in 1:length(levs)){
            group <- levs[j]
            select.obs <- groups == group
            yhat <- matrix(predict(mod, newdata=cbind(dat,
                                                      groups=group)), grid.lines, grid.lines)
            rgl.surface(vals, vals, yhat, color=surface.col[j],
                        alpha=0.5, lit=FALSE)
            if (grid) rgl.surface(vals, vals, yhat,
                                  color=grid.col, alpha=0.5, lit=FALSE, front="lines", back="lines")
            rgl.texts(0, predict(mod, newdata=data.frame(x=0,
                                                         z=0, groups=group)), 0,
                      paste(group, " "), adj=1, color=surface.col[j])
            if (residuals){
              yy <- y[select.obs]
              xx <- x[select.obs]
              zz <- z[select.obs]
              fitted <- fitted(mod)[select.obs]
              rgl.lines(as.vector(rbind(xx,xx)), as.vector(rbind(yy,fitted)), as.vector(rbind(zz,zz)),
                        col=surface.col[j])
            }
          }
        }
        else {
          levs <- levels(groups)
          for (j in 1:length(levs)){
            group <- levs[j]
            select.obs <- groups == group
            mod <- switch(f,
                          linear = lm(y ~ x + z, subset=select.obs),
                          quadratic = lm(y ~ (x + z)^2 + I(x^2) + I(z^2),subset=select.obs),
                          smooth = if (is.null(df.smooth)) gam(y ~ s(x,z), subset=select.obs)
                          else gam(y ~ s(x, z, fx=TRUE, k=df.smooth),subset=select.obs),
                          additive = if (is.null(df.additive)) gam(y ~s(x) + s(z), subset=select.obs)
                          else gam(y ~ s(x, fx=TRUE,k=df.additive[1]+1) +
                                     s(z, fx=TRUE,k=(rev(df.additive+1)[1]+1)), subset=select.obs)
            )
            if (model.summary) summaries[[paste(f, ".", group,sep="")]] <- summary(mod)
            yhat <- matrix(predict(mod, newdata=dat),grid.lines, grid.lines)
            rgl.surface(vals, vals, yhat, color=surface.col[j],alpha=0.5, lit=FALSE)
            rgl.surface(vals, vals, yhat, color=grid.col,alpha=0.5, lit=FALSE, front="lines", back="lines")
            rgl.texts(0, predict(mod, newdata=data.frame(x=0,z=0, groups=group)), 0,
                      paste(group, " "), adj=1, color=surface.col[j])
            if (residuals){
              yy <- y[select.obs]
              xx <- x[select.obs]
              zz <- z[select.obs]
              fitted <- fitted(mod)
              rgl.lines(as.vector(rbind(xx,xx)),as.vector(rbind(yy,fitted)), as.vector(rbind(zz,zz)),
                        col=surface.col[j])
            }
          }
        }
      }
    }
  }
  if (revolutions > 0) {
    for (i in 1:revolutions){
      for (angle in seq(1, 360, length=360/speed))
        rgl.viewpoint(-angle, fov=fov)
    }
    if (model.summary) return(summaries) else return(invisible(NULL))
  }
}
