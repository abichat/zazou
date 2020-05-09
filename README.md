
<!-- README.md is generated from README.Rmd. Please edit that file -->

# zazou

<!-- badges: start -->

<!-- [![Last-changedate](https://img.shields.io/badge/Last%20change-Dec-yellowgreen.svg)]() -->

[![packageversion](https://img.shields.io/badge/Package%20version-0.0.0.9009-orange.svg)]()
[![R build
status](https://github.com/abichat/zazou/workflows/R-CMD-check/badge.svg)](https://github.com/abichat/zazou/actions)
<!-- badges: end -->

The goal of **zazou** is to model the evolution of z-scores along a tree
to do differential analysis.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("abichat/zazou")
```

## Notations

To keep consistent notations in the package:

  - the considered tree has `m` leafs and `n` internal branches,

  - the vector of z-scores `Z` or `zscores` is size `m`,

  - the vector of shift `Delta` or `shifts` is size `n+m`,

  - the incidence matrix `mat_incidence` is size `m*(n+m)`,

  - the OU covariance matrix `mat_covarOU` and its square root `R` are
    size `m*m`,

  - `Y = R %*% Z` is a vector of size `m`,

  - `X = R %*% mat_incidence` is a matrix of size size `m*(n+m)`.
