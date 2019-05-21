#' A function that returns the equation y=ax+b after estimating a and b
#' @description This function calls lm(y~x,db) and constructs an equation with the estimates. Use to add equation to scatterplot
#' @param df : data base
#' @param x : name of the independent variable
#' @param y : name of the dependent variable
#' @examples 
#' lm_eqn(data,'x','y')

lm_eqn = function(df,x,y){
  m = lm(y ~x, data.frame(x=df[,x],y=df[,y]));
  db<-list(a = format(coef(m)[1], digits = 2), 
           b = format(abs(coef(m)[2]), digits = 2))
  eq <- substitute(italic(y) == a + b %.% italic(x), db)
  if (coef(m)[2]<0) {
    eq <- substitute(italic(y) == a - b %.% italic(x), db)
  }
  paste(as.character(as.expression(eq)));                 
}
