% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_helpers.R
\name{fitch}
\alias{fitch}
\title{Fitch parsimony score}
\usage{
fitch(phy, states)
}
\arguments{
\item{phy}{tree}

\item{states}{named state vector (typically the estimated z-scores)}
}
\value{
Fitch parsimony score
}
\description{
Fitch parsimony score
}
\examples{
nwk <- "(((t4:8,(t6:5,t5:5):3):66,t2:74):40,(t7:39,(t1:10,t3:10):29):75);"
tree <- ape::read.tree(text = nwk)
states <- c(t1 = "A", t2 = "B", t3 = "A", t4 = "C",
            t5 = "C", t6 = "C", t7 = "A")
fitch(tree, states) ## 2
}
