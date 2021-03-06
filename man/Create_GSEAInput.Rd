% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/Create_GSEAInput.R
\name{Create_GSEAInput}
\alias{Create_GSEAInput}
\title{A function that creates ranked files for GSEA}
\usage{
Create_GSEAInput(coefs, ann, fname)
}
\arguments{
\item{coefs}{a matrix with coefficients of a linear model}

\item{ann}{is the package with annotation for the platform that is being used}

\item{fname}{is the prefix for the files to be output}
}
\description{
a function that eliminates ambiguity in gene symbols and prepares preRanked Gsea files for different contrasts
depends on getMat and myfindLargest
}
\examples{
prepareGSEArnk(coefs = coefs, Annpkg = hgu133plus2.db, fname='mycontrasts',symbs='NULL')
}

