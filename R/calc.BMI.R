#' BMI
#' @description Given height(cm) and weight(kg), calculatiing BMI 
#' @param height.cm: height(cm)
#' @param weight.kg: weight(kg)


calc.BMI<-function(height.cm,weight.kg){weight.kg/((height.cm/100)^2)}
