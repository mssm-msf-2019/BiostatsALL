% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/normDataWithin.R
\name{normDataWithin}
\alias{normDataWithin}
\title{Summarizes data.
Norms the data within specified groups in a data frame;}
\usage{
normDataWithin(data = NULL, idvar, measurevar, betweenvars = NULL,
  na.rm = FALSE, .drop = TRUE)
}
\arguments{
\item{data:}{a data frame.}

\item{idvar:}{the name of a column that identifies each subject (or matched subjects)}

\item{measurevar:}{the name of a column that contains the variable to be summariezed}

\item{betweenvars:}{a vector containing names of columns that are between-subjects variables}

\item{na.rm:}{a boolean that indicates whether to ignore NA's}
}
\description{
it normalizes each subject (identified by idvar) so that they have the same mean, within each groups pecified by betweenvars.
}

