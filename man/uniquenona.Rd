% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/uniquenona.R
\name{uniquenona}
\alias{uniquenona}
\title{A function that given a vector returns only the unique elements that are not NA or epty characters ("").
 function unique does keep NA in the output}
\usage{
uniquenona(x)
}
\arguments{
\item{x:}{vector of characters}
}
\description{
A function that given a vector returns only the unique elements that are not NA or epty characters (""). Function unique does keep NA in the output
}
\examples{
uniquenona(c(LETTERS[1:10],"",NA))
}

