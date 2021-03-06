% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{do.the.whole.thing}
\alias{do.the.whole.thing}
\title{do.the.whole.thing}
\usage{
do.the.whole.thing(cent, y.col.name = col.name, x.col.name = "AGE",
  group.col.name = "SEX", data, x.pred = NULL, fam = c("BCPE", "BCT",
  "BCCG", "BCPEo", "BCTo", "BCCGo"), n.cyc = 30, refit = F)
}
\arguments{
\item{cent}{the centiles that should be calculated}

\item{y.col.name}{the dependent variable}

\item{x.col.name}{the independent variable}

\item{group.col.name}{group, the whole thing is done for every group seperately}

\item{data}{data frame with at least a dependent an indepenent and a group variable}

\item{x.pred}{x variable for prediction}

\item{n.cyc}{number of iterations for lms}
}
\value{
a list with everything
}
\description{
does what it does
}
