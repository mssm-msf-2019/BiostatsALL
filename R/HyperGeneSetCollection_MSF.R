#' A function that runs GeneSets analysis (ORA and GSEA)  over a set of gen sets and create outputs
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
#' @param doGSEA if TRUE do GSEA type of analysis
#' @param doGSOA if TRUE do Over representation Analysis using HyperGeom test

HyperGeneSetCollection_MSF<-function (listOfGeneSetCollections, geneList, hits, pAdjustMethod = "BH",
                                      pValueCutoff = 1, nPermutations = 1000, minGeneSetSize, fname = NULL,
                                      doGSEA = TRUE, doGSOA = TRUE, entrez2symbs = NULL) {
  library(GSEABase)
  library(HTSanalyzeR)
  library(GO.db)
  library(stringr)
  require(KEGG.db)
  require(plyr)
  Res <- analyzeGeneSetCollections_MSF(listOfGeneSetCollections = listOfGeneSetCollections,
                                       geneList = geneList, hits = hits, pAdjustMethod = pAdjustMethod,
                                       pValueCutoff = pValueCutoff, nPermutations = nPermutations,
                                       minGeneSetSize = minGeneSetSize, exponent = 1, verbose = TRUE,
                                       doGSOA = doGSOA, doGSEA = doGSEA)
  getKEGGNAMES <- function(kegg.ids) {
    kegg.names <- unlist(as.list(KEGGPATHID2NAME)[str_replace(kegg.ids,"hsa", replacement = "")])
    return(kegg.names)
  }
  separateHTTP <- function(x) {
    if (grepl("http",x,fixed=T)) {
      y <- sapply(strsplit(gsub("http",">>http",as.character(x)), " >>",fixed=T),"[[", 2)
    } else {
      y <- NULL
    }
  }
  if (is.null(names(listOfGeneSetCollections))) {
    names.lists <- paste("List", 1:length(listOfGeneSetCollections),
                         sep = "_")
  } else names.lists <- names(listOfGeneSetCollections)
  aux.fcn <- function(g, hits, mapgenes) {
    gns <- uniquenona(sort(intersect(g, hits)))
    if ((!is.null(mapgenes)) & (length(gns) > 0)) {
      gns <- paste(gns, "(", mapgenes[gns], ")", sep = "")
    }
    return(gsub(" ", ", ", do.call("paste", as.list(gns))))
  }
  for (l in c(1:length(listOfGeneSetCollections))) {
    print(paste("preparing 3 save outputs", names.lists[l]))
    if (doGSOA) {
      PTWNAMES <- rownames(Res$HyperGeo.results[[l]])
      if (length(intersect(PTWNAMES, keys(GO.db))) > 0)
        PTWNAMES <- select(GO.db, keys = rownames(Res$HyperGeo.results[[l]]),
                           columns = c("TERM", "ONTOLOGY"))
      if (length(intersect(PTWNAMES, as.list(KEGGPATHID2NAME))) >
          0)
        PTWNAMES <- getKEGGNAMES(rownames(Res$HyperGeo.results[[l]]))
      genesbyPTW <- sapply(listOfGeneSetCollections[[l]],
                           aux.fcn, hits, entrez2symbs)
      print(head(genesbyPTW[1:3]))
      Res$HyperGeo.results[[l]] <- cbind.data.frame(GeneSet = PTWNAMES,
                                                    Res$HyperGeo.results[[l]], Genes = genesbyPTW[rownames(Res$HyperGeo.results[[l]])])
    }
    if (doGSEA) {
      PTWNAMES <- rownames(Res$GSEA.results[[l]])
      if (length(intersect(PTWNAMES, keys(GO.db))) > 0)
        PTWNAMES <- select(GO.db, keys = rownames(Res$HyperGeo.results[[l]]),
                           columns = c("TERM", "ONTOLOGY"))
      if (length(intersect(PTWNAMES, as.list(KEGGPATHID2NAME))) >
          0)
        PTWNAMES <- getKEGGNAMES(rownames(Res$HyperGeo.results[[l]]))
      Res$GSEA.results[[l]] <- cbind.data.frame(GeneSet = PTWNAMES,
                                                Res$GSEA.results[[l]])
    }
  }
  if (!is.null(fname)) {
    print(fname)
    for (l in c(1:length(listOfGeneSetCollections))) {
      if (doGSOA) {
        db <- cbind.data.frame(#Set = rownames(Res$HyperGeo.results[[l]]),
                               Res$HyperGeo.results[[l]])
        db[,"Link"] <- sapply(db$GeneSet, separateHTTP)
        db$GeneSet <- gsub("http(.*)","",db$GeneSet)
        write.csv(file = paste(names.lists[l], fname,
                               "Hyper.csv", sep = "_"), x = plyr::arrange(subset(db,
                                                                                 Pvalue <= pValueCutoff), Pvalue), row.names = TRUE)
      }
      if (doGSEA) {
        db <- cbind.data.frame(#Set = rownames(Res$GSEA.results[[l]]),
                               Res$GSEA.results[[l]])
        db[,"Link"] <- sapply(db$GeneSet, separateHTTP)
        db$GeneSet <- gsub("http(.*)","",db$GeneSet)
        write.csv(file = paste(names.lists[l], fname,
                               "GSEA.csv", sep = "_"), x = plyr::arrange(subset(db,
                                                                                Pvalue <= pValueCutoff), Pvalue), row.names = TRUE)
      }
    }
  }
  return(Res)
}
