#' Open documentation site of the package
#'
#' @importFrom utils browseURL
#'
#' @export
open_doc <- function() {
  guide_path <- system.file('docs/index.html', package = 'zazou')

  if (guide_path == "") stop('There is no pkgdown site in ', 'docs/index.html')

  browseURL(paste0('file://', guide_path))
}
