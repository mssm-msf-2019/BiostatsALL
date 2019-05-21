#' A function that select a set of dissimilar models from caret package
#' @description this function extracts a list of models in caret that are dissimilar
#' @param method is the kind of model Classification / Regression
#' @param N the number of dissimilar models to be extracted
#' @param include a character to separate the names of the columns
#' @keywords caret , regression classification
#' @examples getModels_Caret('Regression',5,'')
#' 

getModels_Caret<-function(method="Regression",N=5,include=''){
  
  
  ###This function tells you what are most N disimilar models available in caret for a certain proble (tag)
  require(proxy)
  require(caret)
  
  tag <-read.csv("tag_data.csv", row.names = 1)
  tag <- as.matrix(tag)
  
  ## Select only models for regression
  tagModels <- tag[tag[,method] == 1,]
  
  all <- 1:nrow(tagModels)
  start <- grep(paste("(",include,")",sep=''), rownames(tagModels), fixed = TRUE)
  pool <- all[all != start]
  
  ## Select 4 model models by maximizing the Jaccard
  ## dissimilarity between sets of models
  nextMods <- maxDissim(tagModels[start,,drop = FALSE],
                        matrix(tagModels[pool, ],ncol=1),
                        method = "Jaccard",
                        n = N)
  
  rownames(tagModels)[c(start, nextMods)]
}
