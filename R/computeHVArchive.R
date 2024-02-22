#' Hypervolume for all prefixes of an archive of points
#' 
#' Implements an efficient (\eqn{n \log n}) algorithm to compute the hypervolume
#' for all prefixes of an archive of objective vectors.
#' 
#' Currently only supports **bi-objective problems**.
#' Assumes **minimization**.
#' 
#' The "historical" hypervolume of an archive is particularly interesting for
#' benchmarking purposes.
#'
#' @param archive Numeric matrix with one pair of objective values per column,
#' and two rows
#' @param ref_point Numeric vector of length two containing the reference point.
#'
#' @return A vector of length `ncol(archive)` containing the hypervolume of all
#' prefixes of the archive.
#' @export
#'
computeHVArchive <- function(archive, ref_point) {
  checkmate::assert_matrix(archive, any.missing = FALSE, nrows = 2L)
  checkmate::assert_numeric(ref_point, any.missing = FALSE, len = 2L)
  
  computeHVArchiveCPP(archive, ref_point)
}