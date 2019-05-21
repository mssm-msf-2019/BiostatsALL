#' A function called from HyperGeneSetCollection_MSF
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
#' @param doGSEA if TRUE do GSEA type of analysis
#' @param doGSOA if TRUE do Over representation Analysis using HyperGeom test
#'
analyzeGeneSetCollections_MSF<-function (listOfGeneSetCollections, geneList, hits, pAdjustMethod = "BH",
                                         pValueCutoff = 0.05, nPermutations = 1000, minGeneSetSize = 15,
                                         exponent = 1, verbose = TRUE, doGSOA = TRUE, doGSEA = TRUE) {
  require(HTSanalyzeR)
  if ((!doGSOA) && (!doGSEA))
    stop("Please choose to perform hypergeometric tests and/or GSEA by specifying 'doGSOA' and 'doGSEA'!\n")
  if (doGSOA || doGSEA) {
  }
  if (doGSOA) {
  }
  if (doGSEA) {
  }
  numGeneSetCollections <- length(listOfGeneSetCollections)
  max.size <- 0
  for (l in 1:numGeneSetCollections) {
    gs.size <- unlist(mclapply(mclapply(listOfGeneSetCollections[[l]],
                                        intersect, y = names(geneList), mc.cores = 4), length,
                               mc.cores = 4))
    max.size <- max(max(gs.size), max.size)
    gs.id <- which(gs.size >= minGeneSetSize)
    n.gs.discarded <- length(listOfGeneSetCollections[[l]]) -
      length(gs.id)
    listOfGeneSetCollections[[l]] <- listOfGeneSetCollections[[l]][gs.id]
    if (verbose && n.gs.discarded > 0)
      cat(paste("--", n.gs.discarded, " gene sets don't have >= ",
                minGeneSetSize, " overlapped genes with universe in gene",
                " set collection named ", names(listOfGeneSetCollections)[l],
                "!\n", sep = ""))
  }
  if (all(unlist(mclapply(listOfGeneSetCollections, length,
                          mc.cores = 4)) == 0))
    stop(paste("No gene set has >= ", minGeneSetSize, " overlapped ",
               "genes with universe!\n The largest number of overlapped ",
               "genes between gene sets and universe is: ", max.size,
               sep = ""))
  result.names <- names(listOfGeneSetCollections)
  if (doGSOA) {
    HGTresults <- mclapply(listOfGeneSetCollections, function(L) {
      if (length(L) > 0)
        out <- multiHyperGeoTest(L, universe = names(geneList),
                                 hits = hits, minGeneSetSize = minGeneSetSize,
                                 pAdjustMethod = pAdjustMethod, verbose = verbose)
      else {
        out <- matrix(, nrow = 0, ncol = 7)
        colnames(out) <- c("Universe Size", "Gene Set Size",
                           "Total Hits", "Expected Hits", "Observed Hits",
                           "Pvalue", "Adjusted.Pvalue")
      }
      out<-cbind(out[,1:5],OR=out[,"Observed Hits"]/out[, "Expected Hits"],out[,6:7])
      return(out)
    }, mc.cores = 4)
    pvals <- NULL
    sapply(1:numGeneSetCollections, function(i) {
      if (nrow(HGTresults[[i]]) > 0) {
        pv <- HGTresults[[i]][, "Pvalue"]
        names(pv) <- rownames(HGTresults[[i]])
        pvals <<- c(pvals, pv)
      }
    })
    HGTpvals <- p.adjust(pvals, method = pAdjustMethod)
    sapply(1:numGeneSetCollections, function(i) {
      if (nrow(HGTresults[[i]]) > 0) {
        ind <- match(rownames(HGTresults[[i]]), names(HGTpvals))
        Adjusted.Pvalue <- HGTpvals[ind]
        HGTresults[[i]][, "Adjusted.Pvalue"] <<- Adjusted.Pvalue
      }
    })
    names(HGTresults) <- result.names
    cat("-Hypergeometric analysis complete\n\n")
    sign.hgt <- mclapply(HGTresults, function(x) {
      if (nrow(x) > 0) {
        a <- which(x[, "Pvalue"] < pValueCutoff)
        x <- x[a, , drop = FALSE]
        if (length(a) > 1) {
          x <- x[order(x[, "Pvalue"]), , drop = FALSE]
        }
        return(x)
      }
    }, mc.cores = 4)
    sign.hgt.adj <- mclapply(HGTresults, function(x) {
      if (nrow(x) > 0) {
        a <- which(x[, "Adjusted.Pvalue"] < pValueCutoff)
        x <- x[a, , drop = FALSE]
        if (length(a) > 1) {
          x <- x[order(x[, "Adjusted.Pvalue"]), , drop = FALSE]
        }
        return(x)
      }
    }, mc.cores = 4)
  }
  else {
    HGTresults = NULL
  }
  if (doGSEA) {
    GSEA.results.list <- list()
    cat("-Performing gene set enrichment analysis ...\n")
    test.collection <- mclapply(listOfGeneSetCollections,
                                function(L) {
                                  if (verbose) {
                                    cat("--For", "names(listOfGeneSetCollections)[i]",
                                        "\n")
                                  }
                                  if (length(L) > 0)
                                    collectionGsea(L, geneList = geneList, exponent = exponent,
                                                   nPermutations = nPermutations, minGeneSetSize = minGeneSetSize,
                                                   verbose = verbose)
                                  else {
                                    test.collection[[i]] <<- list(Observed.scores = NULL,
                                                                  Permutation.scores = NULL)
                                  }
                                }, mc.cores = 4)
    obs.scores <- NULL
    prm.scores <- NULL
    sapply(1:numGeneSetCollections, function(i) {
      obs.scores <<- c(obs.scores, test.collection[[i]]$Observed.scores)
      prm.scores <<- rbind(prm.scores, test.collection[[i]]$Permutation.scores)
    })
    if (length(obs.scores) > 0) {
      total.test.collect <- list(Observed.scores = obs.scores,
                                 Permutation.scores = prm.scores)
      test.FDR.collection <- FDRcollectionGsea(permScores = total.test.collect$Permutation.scores,
                                               dataScores = total.test.collect$Observed.scores)
      test.pvalues.collection <- permutationPvalueCollectionGsea(permScores = total.test.collect$Permutation.scores,
                                                                 dataScores = (total.test.collect$Observed.scores))
      gsea.adjust.pval <- p.adjust(test.pvalues.collection,
                                   method = pAdjustMethod)
      test.GSEA.results <- cbind(total.test.collect$Observed.scores,
                                 test.pvalues.collection, gsea.adjust.pval, test.FDR.collection)
      colnames(test.GSEA.results) <- c("Observed.score",
                                       "Pvalue", "Adjusted.Pvalue", "FDR")
      sapply(1:numGeneSetCollections, function(i) {
        if (!is.null(test.collection[[i]])) {
          match.ind <- match(names(listOfGeneSetCollections[[i]]),
                             rownames(test.GSEA.results))
          match.ind <- match.ind[which(!is.na(match.ind))]
          GSEA.res.mat <- test.GSEA.results[match.ind,
                                            , drop = FALSE]
          GSEA.res.mat <- GSEA.res.mat[order(GSEA.res.mat[,
                                                          "Adjusted.Pvalue"]), , drop = FALSE]
          GSEA.results.list[[i]] <<- GSEA.res.mat
        }
      })
      names(GSEA.results.list) <- result.names
    }
    else {
      sapply(1:numGeneSetCollections, function(i) {
        GSEA.results.list[[i]] <<- matrix(, nrow = 0,
                                          ncol = 4)
        colnames(GSEA.results.list[[i]]) <<- c("Observed.score",
                                               "Pvalue", "Adjusted.Pvalue", "FDR")
      })
    }
    sign.gsea <- lapply(GSEA.results.list, function(x) {
      if (nrow(x) > 0) {
        a <- which(x[, "Pvalue"] < pValueCutoff)
        x <- x[a, , drop = FALSE]
        if (length(a) > 1) {
          x <- x[order(x[, "Pvalue"]), , drop = FALSE]
        }
        return(x)
      }
    })
    sign.gsea.adj <- lapply(GSEA.results.list, function(x) {
      if (nrow(x) > 0) {
        a <- which(x[, "Adjusted.Pvalue"] <= pValueCutoff)
        x <- x[a, , drop = FALSE]
        if (length(a) > 1) {
          x <- x[order(x[, "Adjusted.Pvalue"]), , drop = FALSE]
        }
        return(x)
      }
    })
  }
  else {
    GSEA.results.list = NULL
  }
  if (doGSOA && doGSEA) {
    overlap <- list()
    overlap.adj <- list()
    sapply(1:numGeneSetCollections, function(i) {
      a1 <- intersect(rownames(sign.gsea[[i]]), rownames(sign.hgt[[i]]))
      a2 <- intersect(rownames(sign.gsea.adj[[i]]), rownames(sign.hgt.adj[[i]]))
      Hypergeometric.Pvalue <- HGTresults[[i]][a1, "Pvalue", drop = FALSE]
      Hypergeometric.Adj.Pvalue <- HGTresults[[i]][a2, "Adjusted.Pvalue", drop = FALSE]
      GSEA.Pvalue <- GSEA.results.list[[i]][a1, "Pvalue", drop = FALSE]
      GSEA.Adj.Pvalue <- GSEA.results.list[[i]][a2, "Adjusted.Pvalue", drop = FALSE]
      overlap[[i]] <<- cbind(Hypergeometric.Pvalue, GSEA.Pvalue)
      colnames(overlap[[i]]) <<- c("HyperGeo.Pvalue", "GSEA.Pvalue")
      overlap.adj[[i]] <<- cbind(Hypergeometric.Adj.Pvalue,
                                 GSEA.Adj.Pvalue)
      colnames(overlap.adj[[i]]) <<- c("HyperGeo.Adj.Pvalue",
                                       "GSEA.Adj.Pvalue")
    })
    names(overlap) <- result.names
    names(overlap.adj) <- result.names
  }
  else {
    overlap = NULL
    overlap.adj = NULL
  }
  cat("-Gene set enrichment analysis complete \n")
  final.results <- list(HyperGeo.results = HGTresults, GSEA.results = GSEA.results.list,
                        Sig.pvals.in.both = overlap, Sig.adj.pvals.in.both = overlap.adj)
  print("Finish Running GSOA/GSEA")
  return(final.results)
}
