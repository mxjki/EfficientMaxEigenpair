#' @title Tridiagonal matrix next to the maximal eigenpair
#'
#' @description Calculate the next to maximal eigenpair for the tridiagonal matrix
#' whose sums of each row should be 0.
#'
#' @param a The lower diagonal vector.
#' @param b The upper diagonal vector.
#' @param xi The coefficient used in the improved initials to form
#'        the convex combination of \eqn{\delta_1^{-1}} and \eqn{(v_0,-Q*v_0)_\mu},
#'        it should between 0 and 1.
#' @param digit.thresh The precise level of output results.
#'
#' @return A list of eigenpair object are returned, with components \eqn{z} and \eqn{v}.
#' \item{z}{The approximating sequence of the maximal eigenvalue.}
#' \item{v}{The approximating sequence of the corresponding eigenvector.}
#'
#' @note The sums of each row of the input tridiagonal matrix should be 0.
#'
#' @seealso \code{\link{eff.ini.seceig.general}} for the general conservative matrix next to the maximal eigenpair.
#'
#' @examples
#' a = c(1:7)^2
#' b = c(1:7)^2
#'
#' eff.ini.seceig.tri(a, b, xi = 0)
#' eff.ini.seceig.tri(a, b, xi = 1)
#' eff.ini.seceig.tri(a, b, xi = 2/5)

#' @export
eff.ini.seceig.tri = function(a, b, xi = 1, digit.thresh = 6) {
    N = length(a)
    Q = tridiag(b, a, -c(b[1], a[1:N - 1] + b[2:N], a[N]))

    # check input vectors
    if (sum(a <= 0) > 0)
        stop("Input vector a should be all postive!")
    if (sum(b <= 0) > 0)
        stop("Input vector b should be all postive!")
    if (xi < 0 | xi > 1)
        stop("The coefficient xi should between 0 and 1!")

    mu = rep(1, N + 1)
    mu[2:(N + 1)] = cumprod(b/a)

    pivec = mu/sum(mu)

    phi = rep(0, N + 1)
    phi[2:(N + 1)] = cumsum(1/(mu[1:N] * b))

    v0_tilde = sqrt(phi)
    v0_bar = v0_tilde - sum(pivec * v0_tilde)
    v0 = v0_bar/sqrt(sum(v0_bar^2 * mu))

    etapart = mu * v0_bar
    etapartsum = rev(cumsum(rev(etapart)))
    eta1 = max(1/(mu[1:N] * b[1:N] * diff(v0_tilde)) * etapartsum[2:(N +
        1)])

    zstart = xi * 1/eta1 + (1 - xi) * sum(v0_bar * (-Q %*% v0_tilde) *
        mu)/(sum(v0_bar^2 * mu))

    ray = ray.quot(Q = Q, mu = mu, v0_tilde = v0_bar, zstart = zstart,
        digit.thresh = digit.thresh)

    return(list(z = unlist(ray$z), v = ray$v))
}
