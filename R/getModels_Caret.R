#' A function that tells you what are most N disimilar models available in caret for a certain problem (tag)
#' @description It selects the N most dissimilar models among those in caret suite by maximizing the Jaccard dissimilarity between sets of models.
#' @param method: type of methods. Options are Regression, Classification adn Dual use
#' @param N: number of methods to select
#' @param include: name of a method to be included in the set
#' @examples 
#' getModels_Caret<-function(method="Regression",N=5,include='cubist')


getModels_Caret<-function(method="Regression",N=5,include=''){
  ###This function tells you what are most N disimilar models available in caret for a certain proble (tag)
  
  tag<-read.csv(file='/Users/mayte/Rlibrary/caret-model-tags.csv', row.names = 1)
  tag <- as.matrix(tag)
  
  ## Select only models for regression
  tagModels <- tag[tag[,method] == 1,]
  
  all <- 1:nrow(tagModels)
  start <- grep(paste("(",include,")",sep=''), rownames(tagModels), fixed = TRUE)
  pool <- all[all != start]
  
  ## Select 4 model models by maximizing the Jaccard
  ## dissimilarity between sets of models
  nextMods <- maxDissim(tagModels[start,,drop = FALSE],
                        tagModels[pool, ],
                        method = "Jaccard",
                        n = N)
  
  rownames(tagModels)[c(start, nextMods)]
}
