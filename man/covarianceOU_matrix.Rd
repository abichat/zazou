% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariance_matrix.R
\name{covarianceOU_matrix}
\alias{covarianceOU_matrix}
\title{Covariance matrix}
\usage{
covarianceOU_matrix(tree, alphaOU)
}
\arguments{
\item{tree}{a phylo object}

\item{alphaOU}{hyperparameter}
}
\value{
the covariance matrix of the OU process on the tree
}
\description{
Compute normalized covariance matrix
}
\details{
The coefficient is proportional to
\eqn{(exp(- \alpha * d_{i,j}) - exp(-2 * \alpha * h)) / (2 * \alpha)}
}
\examples{
tree <- ape::rtree(3)
covarianceOU_matrix(tree, alpha = 1)
}
