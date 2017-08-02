#' @title Rayleigh quotient iteration
#' @description Rayleigh quotient iteration algorithm to computing the maximal eigenpair of
#' matrix Q.
#'
#' @param Q The input matrix to find the maximal eigenpair.
#' @param mu A vector.
#' @param v0_tilde The unnormalized initial vector \eqn{\tilde{v0}}.
#' @param zstart The initial \eqn{z_0} as an approximation of \eqn{\rho(Q)}.
#' @param digit.thresh The precise level of output results.
#' @return A list of eigenpair object are returned, with components \eqn{z}, \eqn{v} and \eqn{iter}.
#' \item{z}{The approximating sequence of the maximal eigenvalue.}
#' \item{v}{The approximating sequence of the corresponding eigenvector.}
#' \item{iter}{The number of iterations.}
#'
#' @examples
#' Q = matrix(c(1, 1, 3, 2, 2, 2, 3, 1, 1), 3, 3)
#' ray.quot.seceig.general(Q, mu=rep(1,dim(Q)[1]), v0_tilde=rep(1,dim(Q)[1]), zstart=6,
#'  digit.thresh = 6)

#' @export
ray.quot.seceig.general = function(Q, mu, v0_tilde, zstart, digit.thresh = 6) {
    z = list()
    rz = list()
    v = list()
    w = list()

    v[[1]] = v0_tilde/sqrt(sum(v0_tilde^2 * mu))

    z[[1]] = zstart
    rz[[1]] = round(zstart, digit.thresh)

    ratio = 1
    iter = 0
    while (ratio >= 10^(-digit.thresh)) {
        iter = iter + 1
        w = append(w, list(solve(-Q - z[[iter]] * diag(1, length(v0_tilde)),
            v[[iter]], tol = 1e-100)))
        v = append(v, list(w[[iter]]/sqrt(sum(w[[iter]]^2 * mu))))

        z = append(z, list(sum(v[[iter + 1]] * (-Q %*% v[[iter +
            1]]) * mu)))

        ratio = abs(round(z[[iter + 1]], digit.thresh) - round(z[[iter]],
            digit.thresh))

        rz[[iter + 1]] = round(z[[iter + 1]], digit.thresh)
    }

    if (ratio == 0) {
        v = v[-(iter + 1)]
        rz = rz[-(iter + 1)]
        
        iter = iter - 1
    }

    return(list(v = v, z = rz, iter = iter))
}
