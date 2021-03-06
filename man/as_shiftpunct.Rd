% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiftpunct.R
\name{as_shiftpunct}
\alias{as_shiftpunct}
\alias{print.shiftpunct}
\alias{plot.shiftpunct}
\title{'shiftpunct' object}
\usage{
as_shiftpunct(listopt, tree, zscores, lambda, alphaOU)

\method{print}{shiftpunct}(x, digits = 3, ...)

\method{plot}{shiftpunct}(x, digits = 3, ...)
}
\arguments{
\item{listopt}{an output of 'estimate_shifts()'}

\item{tree}{tree}

\item{zscores}{vector of observed z-scores}

\item{lambda}{regularization parameter}

\item{alphaOU}{parameter}

\item{x}{a 'shiftpunct' object.}

\item{digits}{Number of digits to round to.}

\item{...}{further arguments to be passed to or from other methods.}
}
\value{
a 'shiftpunct' object
}
\description{
'shiftpunct' object
}
