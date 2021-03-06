#' EfficientMaxEigenpair: A package for computating the maximal eigenpair for a matrix.
#'
#' The EfficientMaxEigenpair package provides some auxillary functions and
#' five categories of important functions:
#' \code{\link{tridiag}}, \code{\link{tri.sol}}, \code{\link{find_deltak}},
#' \code{\link{ray.quot.tri}}, \code{\link{shift.inv.tri}},
#' \code{\link{ray.quot.seceig.tri}}, \code{\link{ray.quot.general}}, \code{\link{ray.quot.seceig.general}},
#' \code{\link{eff.ini.maxeig.tri}}, \code{\link{eff.ini.maxeig.shift.inv.tri}},
#' \code{\link{eff.ini.maxeig.general}}, \code{\link{eff.ini.seceig.tri}} 
#' and \code{\link{eff.ini.seceig.general}}.
#'
#' @section EfficientMaxEigenpair functions:
#' \code{\link{tridiag}}: generate tridiagonal matrix Q based on three input vectors.
#'
#' \code{\link{tri.sol}}: construct the solution of linear equation (-Q-zI)w=v.
#'
#' \code{\link{find_deltak}}: compute \eqn{\delta_k} for given vector \eqn{v} and matrix \eqn{Q}.
#'
#' \code{\link{ray.quot.tri}}: rayleigh quotient iteration algorithm to computing the maximal eigenpair of
#' tridiagonal matrix \eqn{Q}.
#' 
#' \code{\link{shift.inv.tri}}: shifted inverse iteration algorithm to computing the maximal eigenpair of
#' tridiagonal matrix \eqn{Q}.
#' 
#' \code{\link{ray.quot.seceig.tri}}: rayleigh quotient iteration algorithm to computing the next to 
#' maximal eigenpair of tridiagonal matrix \eqn{Q}.
#'
#' \code{\link{ray.quot.general}}: rayleigh quotient iteration algorithm to computing the maximal eigenpair of
#' general matrix \eqn{A}.
#' 
#' \code{\link{ray.quot.seceig.general}}: rayleigh quotient iteration algorithm to computing the  next to maximal eigenpair of
#' general matrix \eqn{A}.
#' 
#' \code{\link{eff.ini.maxeig.tri}}: calculate the maximal eigenpair for the tridiagonal matrix by
#' rayleigh quotient iteration algorithm.
#' 
#' \code{\link{eff.ini.maxeig.shift.inv.tri}}: calculate the maximal eigenpair for the tridiagonal matrix by
#' shifted inverse iteration algorithm.
#' 
#' \code{\link{eff.ini.maxeig.general}}: calculate the maximal eigenpair for the general matrix.
#'
#' \code{\link{eff.ini.seceig.tri}}: calculate the next to maximal eigenpair for the tridiagonal matrix whose sums of each row should be 0.
#'
#' \code{\link{eff.ini.seceig.general}}: calculate the next to maximal eigenpair for the general conservative matrix.
#'
#' @docType package
#' @name EfficientMaxEigenpair
NULL
