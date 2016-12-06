#' EfficientMaxEigenpair: A package for computating the maximal eigenpair for a matrix.
#'
#' The EfficientMaxEigenpair package provides some auxillary functions and
#'  four categories of important functions:
#'  \code{\link{tridiag}}, \code{\link{ray.quot}}
#' \code{\link{eff.ini.maxeig.tri}}, \code{\link{eff.ini.maxeig.general}},
#' \code{\link{eff.ini.seceig.tri}} and \code{\link{eff.ini.seceig.general}}.
#'
#' @section EfficientMaxEigenpair functions:
#' \code{\link{tridiag}}: Generate tridiagonal matrix Q based on three input vectors.
#'
#' \code{\link{ray.quot}}: Rayleigh quotient iteration algorithm to computing the maximal eigenpair of
#' matrix Q.
#'
#' \code{\link{eff.ini.maxeig.tri}}: Calculate the maximal eigenpair for the tridiagonal matrix.
#'
#' \code{\link{eff.ini.maxeig.general}}: Calculate the maximal eigenpair for the general matrix.
#'
#' \code{\link{eff.ini.seceig.tri}}: Calculate the next to maximal eigenpair for the tridiagonal matrix whose sums of each row should be 0.
#'
#' \code{\link{eff.ini.seceig.general}}: Calculate the next to maximal eigenpair for the general conservative matrix.
#'
#' @docType package
#' @name EfficientMaxEigenpair
NULL
