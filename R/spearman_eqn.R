#' A function that returns an equation with the pearson correlation coef and the pvalue
#' @description This function calls with(db,cor(x,y,method='spearman')) and constructs an equation with the estimates. Use to add equation to scatterplot
#' @param df : data base
#' @param x : name of the independent variable
#' @param y : name of the dependent variable
#' @examples
#' spearman_eqn(data,'x','y')
#'
#'


spearman_eqn<-function(df,x,y) {
  cp <- with(data.frame(x = df[, x], y = df[, y]), cor.test(x, y, use = "p", method = "spearman"))
  db <- list(rc = format(cp$estimate, digits = 2),
             pc = format(cp$p.value, digits = 2),
             pc.exponent = floor(log10(cp$p.value)),
             pc.base = round(cp$p.value/(10^floor(log10(cp$p.value))), digits = 1))
  if (cp$p.value == 0)
    eq <- substitute(~italic(rho) ~ "=" ~ rc ~ ", " ~ italic(p) ~ "<" * 10^-16, db)
  else if (db[["pc.base"]] == 1)
    eq <- substitute(~italic(rho) ~ "=" ~ rc ~ ", " ~ italic(p) ~ "=" ~ 10^pc.exponent, db)
  else if (db[["pc.exponent"]] == (-1))
    eq <- substitute(~italic(rho) ~ "=" ~ rc ~ ", " ~ italic(p) ~ "=" ~ pc, db)
  else eq <- substitute(~italic(rho) ~ "=" ~ rc ~ ", " ~ italic(p) ~ "=" ~ pc.base %*% 10^pc.exponent, db)
  paste(as.character(as.expression(eq)))
}

# spearman_eqn = function(df,x,y){
#   cp<-with(data.frame(x=df[,x],y=df[,y]),cor.test(x,y,use='p',method='spearman'))
#   db<-list(rc = format(cp$estimate, digits = 2),
#            pc = format(cp$p.value, digits = 2),
#            pc.exponent=floor(log10(cp$p.value)),
#            pc.base =round(cp$p.value/(10^floor(log10(cp$p.value))), digits=1)
#   )
#   if  (cp$p.value==0)  eq <- substitute(~italic(r)~"="~rc~", "~italic(p)~"<"*10^-16, db)
#   else  if  (db[['pc.base']]==1)  eq <- substitute(~italic(rho)~"="~rc~", "~italic(p)~"="~10^pc.exponent, db)
#   else  if  (db[['pc.exponent']]==(-1))  eq <- substitute(~italic(r)~"="~rc~", "~italic(p)~"="~pc, db)
#   else   eq <- substitute(~italic(rho)~"="~rc~", "~italic(p)~"="~pc.base%*%10^pc.exponent, db)
#   paste(as.character(as.expression(eq)));
# }
