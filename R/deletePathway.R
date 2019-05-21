#' A function that deletes pathways whose frequency is greater than 5
#' @description a function that examines a frequency table and conserves the elements 
#' with frequency greater than min and less than max parameters
#' @param SET a data frame 
#' @param MyGroups is the data frame column that contains the names of pathways
#' @param max is the max acceptable frequency
#' @param min is the minimum acceptable frequency
#' @examples 
#' deletePathway(SET,MyGroups=mypathways,max=9999,min=0)

deletePathway<-function(SET,MyGroups=mypathways,max=9999,min=0){
  SET<-SET[(SET[,MyGroups]%in%names(table(SET))[(table(SET)>min & table(SET)<max)]),]
}

