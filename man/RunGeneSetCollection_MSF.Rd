% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/RunGeneSetCollection_MSF.R
\name{RunGeneSetCollection_MSF}
\alias{RunGeneSetCollection_MSF}
\title{A function that runs analyzeGeneSetCollections function}
\usage{
RunGeneSetCollection_MSF(listOfGeneSetCollections, geneList, hits,
  pAdjustMethod = "BH", pValueCutoff = 1, nPermutations = 1000,
  minGeneSetSize, fname = NULL)
}
\arguments{
\item{listOfGeneSetCollections}{a list of gene set collections (a 'gene set collection' is a list of gene sets).
Even if only one collection is being tested, it must be entered as an element of a 1-element list,
e.g. ListOfGeneSetCollections = list(YourOneGeneSetCollection).
Naming the elements of listOfGeneSetCollections will result in these names being associated
with the relevant data frames in the output (meaningful names are advised)}

\item{geneList}{a numeric or integer vector of phenotypes in descending or ascending order
with elements named by their EntrezIds (no duplicates nor NA values)}

\item{hits}{a character vector of the EntrezIds of hits, as determined by the user}

\item{pAdjustMethod}{a single character value specifying the p-value adjustment method to be used (see 'p.adjust' for details)}

\item{pValueCutoff}{numeric value, p value of cut of}

\item{nPermutations}{numeric value, number of permutations}

\item{minGeneSetSize}{numeric value, minimun number of genes accepted}

\item{fname}{character string, filename}
}
\description{
This function runs a function that takes a list of gene set collections, a named phenotype vector
             (with names of the phenotype vector as the universe), a vector of hits (gene names only)
             and returns the results of hypergeometric and gene set enrichment analyses for all of the gene set collections
             (with multiple hypothesis testing corrections).
}

