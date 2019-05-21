#' A function that runs analyzeGeneSetCollections function
#' @description This function runs a function that takes a list of gene set collections, a named phenotype vector 
#'              (with names of the phenotype vector as the universe), a vector of hits (gene names only) 
#'              and returns the results of hypergeometric and gene set enrichment analyses for all of the gene set collections 
#'              (with multiple hypothesis testing corrections).
#' @param listOfGeneSetCollections a list of gene set collections (a 'gene set collection' is a list of gene sets). 
#'        Even if only one collection is being tested, it must be entered as an element of a 1-element list, 
#'        e.g. ListOfGeneSetCollections = list(YourOneGeneSetCollection). 
#'        Naming the elements of listOfGeneSetCollections will result in these names being associated 
#'        with the relevant data frames in the output (meaningful names are advised)
#' @param geneList a numeric or integer vector of phenotypes in descending or ascending order 
#'        with elements named by their EntrezIds (no duplicates nor NA values)
#' @param hits a character vector of the EntrezIds of hits, as determined by the user
#' @param pAdjustMethod a single character value specifying the p-value adjustment method to be used (see 'p.adjust' for details) 
#' @param pValueCutoff numeric value, p value of cut of
#' @param nPermutations numeric value, number of permutations
#' @param minGeneSetSize numeric value, minimun number of genes accepted
#' @param fname character string, filename

RunGeneSetCollection_MSF<- function(listOfGeneSetCollections, geneList, hits, pAdjustMethod="BH", pValueCutoff=1, 
                                    nPermutations=1000, minGeneSetSize, fname=NULL, doGSEA=TRUE, doGSOA=TRUE,entrez2symbs){
  print('DEG')
  Res<- HyperGeneSetCollection_MSF(listOfGeneSetCollections, geneList, hits, pAdjustMethod, pValueCutoff, nPermutations,
                                   minGeneSetSize, fname=paste(fname,'Diff',sep='.'), doGSEA=doGSEA, doGSOA=doGSOA,
                                   entrez2symbs=entrez2symbs)
  
  hitsUp<-hits[which(geneList[hits]>0)]
  print('UP')
  ResUp<-HyperGeneSetCollection_MSF(listOfGeneSetCollections, geneList, hitsUp, pAdjustMethod, pValueCutoff, 
                                    nPermutations=10, minGeneSetSize, fname=paste(fname,'Up',sep='.'),
                                    doGSEA=doGSEA, doGSOA=doGSOA, entrez2symbs=entrez2symbs)
  
  hitsDown<-hits[which(geneList[hits]<0)]
  #geneList_scrambled<-geneList+runif(length(geneList),-max(abs(geneList_scrambled)),max(abs(geneList_scrambled)))
  #geneList_scrambled<-geneList_scrambled[order(geneList_scrambled,decreasing=TRUE)]
  print('Down')
  ResDown<-HyperGeneSetCollection_MSF(listOfGeneSetCollections, geneList, hitsDown, pAdjustMethod, pValueCutoff, 
                                      nPermutations=10, minGeneSetSize, fname=paste(fname,'Down',sep='.'),
                                      doGSEA=doGSEA, doGSOA=doGSOA, entrez2symbs=entrez2symbs)
  
  TheRes<-c(Res,HyperGeo.results.Up=ResUp$HyperGeo.results, HyperGeo.results.Down=ResDown$HyperGeo.results)
  return(TheRes)
  
}