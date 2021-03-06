% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint_global.R
\name{estimate_confint}
\alias{estimate_confint}
\title{Estimate confidence interval from scaled lasso.}
\usage{
estimate_confint(
  shiftpunct,
  alpha_conf = 0.05,
  method = c("score system", "colwise inverse"),
  ...
)
}
\arguments{
\item{shiftpunct}{a \code{shiftpunct} object with method \code{scaledlasso}.}

\item{alpha_conf}{the confidence level.}

\item{method}{\code{score system} or \code{colwise inverse}.}

\item{...}{further argument...}
}
\value{
a \code{shiftconf} object.
}
\description{
Estimate confidence interval from scaled lasso.
}
