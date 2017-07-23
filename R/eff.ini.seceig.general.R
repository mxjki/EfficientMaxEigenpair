getzstart = function(Q, mu, x, r, N) {
    v0_tilde = c(r, sqrt(1 - x[2:N]))
    v0_bar = v0_tilde - sum(mu * v0_tilde)/sum(mu)
    v0 = v0_bar/sqrt(sum(v0_bar^2 * mu))
    zstart = sum(v0_bar * (-Q %*% v0_tilde) * mu)/sum(v0_bar^2 * mu)
    return(zstart)
}

#' @title General conservative matrix maximal eigenpair
#'
#' @description Calculate the next to maximal eigenpair for the general conservative matrix.
#'
#' @param Q The input general matrix.
#' @param z0 The type of initial \eqn{z_0} used to calculate the approximation of \eqn{\rho(Q)}.
#'        There are two types: 'fixed' and 'Auto' corresponding to two choices
#'        of \eqn{z_0} in paper.
#' @param c1 A large constant.
#' @param digit.thresh The precise level of output results.
#'
#' @return A list of eigenpair object are returned, with components \eqn{z}, \eqn{v} and \eqn{iter}.
#' \item{z}{The approximating sequence of the maximal eigenvalue.}
#' \item{v}{The approximating sequence of the corresponding eigenvector.}
#' \item{iter}{The number of iterations.}
#' 
#' @note The conservativity of matrix \eqn{Q=(q_{ij})} means that the sums of each row of
#' matrix \eqn{Q} are all 0.
#'
#' @seealso \code{\link{eff.ini.seceig.tri}} for the tridiagonal matrix next to the maximal eigenpair.
#'
#' @examples
#' Q = matrix(c(-30, 1/5, 11/28, 55/3291, 30, -17, 275/42, 330/1097,
#' 0, 84/5, -20, 588/1097, 0, 0, 1097/84, -2809/3291), 4, 4)
#' eff.ini.seceig.general(Q, z0 = 'Auto', digit.thresh = 5)
#' eff.ini.seceig.general(Q, z0 = 'fixed', digit.thresh = 5)

#' @importFrom stats

#' @export
eff.ini.seceig.general = function(Q, z0 = NULL, c1 = 1000, digit.thresh = 6) {
    if (sum(rowSums(Q) > 1e-10) > 0)
        stop("Input matrix Q should be a conservative matrix!")

    N = dim(Q)[1]

    Q1 = Q
    Q1[N, N] = c1 * Q[N, N]
    q1 = 1/diag(Q1)
    P = diag(q1, N) %*% Q1 + diag(1, N)

    if (z0 == "fixed") {
        lambda0Q1 = -eff.ini.maxeig.general(Q1, digit.thresh = digit.thresh,
            improved = T, xi = 1)$z
    }

    x = rep(NA, N)
    x[1] = 1
    x[2:N] = solve(diag(1, N - 1) - P[-1, -1], P[2:N, 1])

    mu = rep(1, N)
    mu[2:N] = solve(t(Q)[-N, -1], -t(Q)[1:(N - 1), 1])

    z0func = function(r) {
        return(getzstart(Q, mu, x, r, N))
    }
    ropt = optim(0.9, z0func, method = c("L-BFGS-B"), lower = 0,
        upper = 1)$par

    v0_tilde = c(ropt, sqrt(1 - x[2:N]))
    v0_bar = v0_tilde - sum(mu * v0_tilde)/sum(mu)
    v0 = v0_bar/sqrt(sum(v0_bar^2 * mu))

    zstart = switch(z0, fixed = lambda0Q1[length(lambda0Q1)], Auto = sum(v0_bar *
        (-Q %*% v0_tilde) * mu)/sum(v0_bar^2 * mu))
    ray = ray.quot(Q = Q, mu = mu, v0_tilde = v0_bar, zstart = zstart,
        digit.thresh = digit.thresh)

    return(list(z = unlist(ray$z), v = ray$v, iter = ray$iter))
}
