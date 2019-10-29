#' Datasets
#'
#' \code{alcohol} and \code{chlamydiae} datasets
#'
#' The \code{alcohol} data set was taken from a study of the dietary
#' association with the gut microbiome (Wu, 2011). Here we include the
#' alcohol intake as the primary phenotype. The OTU abundance data consists
#' of 98 samples with 949 OTUs with prevalence > 10%. The raw counts are
#' normalized to account for variable library sizes.
#'
#' The \code{chlamydia} data set is taken from a 16s rRNA diversity study
#' narrowed to the Chlamydiae bacteria taxon, originally appearing in PNAS.
#'
#' @format A list containing
#' \describe{
#'   \item{\code{X}}{the normalized OTU data,}
#'   \item{\code{Y}}{the alcohol intake phenotype,}
#'   \item{\code{tree}}{the corresponding phylogenetic tree,}
#'   \item{\code{taxonomy}}{the taxonomic lineages of the OTUs.}
#' }
#' @source Wu, Gary D., et al. "Linking long-term dietary patterns with gut
#' microbial enterotypes." Science 334.6052 (2011): 105-108.
#' @source Caporaso, J. G., et al. (2011). "Global patterns of 16S rRNA
#' diversity at a depth of millions of sequences per sample." PNAS,
#' 108, 4516-4522. PMCID: PMC3063599
"alcohol"

#' @rdname alcohol
"chlamydiae"
