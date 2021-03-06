% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/processGMT.R
\name{processGMT}
\alias{processGMT}
\title{Function for parsing a GMT gene set list}
\usage{
processGMT(gene_sets, trim_snam = F)
}
\arguments{
\item{gene_sets}{Name of geneset list that was created by "readGmt"}

\item{trim_snam}{TRUE/FALSE if the geneset.names should be trimmed using "str_trim" from the "stringr" package. Defaults to FALSE.}
}
\description{
This function takes a GMT gene set list (from "readGmt") as input, and outputs a list of gene sets, with names being the geneset.names, and values the corresponding genesets. Note this function also works on a GSA .gmt file object.
}
\examples{
jkr.symbol <- processGMT(jkr.list)
}
\keyword{gene}
\keyword{gmt}
\keyword{sets,}

