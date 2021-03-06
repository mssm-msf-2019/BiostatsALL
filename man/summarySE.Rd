% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/summarySE.R
\name{summarySE}
\alias{summarySE}
\title{Summarizes data.}
\usage{
summarySE(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
  conf.interval = 0.95, .drop = TRUE)
}
\arguments{
\item{data:}{a data frame.}

\item{measurevar:}{the name of a column that contains the variable to be summariezed}

\item{groupvars:}{a vector containing names of columns that contain grouping variables}

\item{na.rm:}{a boolean that indicates whether to ignore NA's}

\item{conf.interval:}{the percent range of the confidence interval (default is 95%)}
}
\description{
Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
}

