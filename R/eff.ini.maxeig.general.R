#' @title General matrix maximal eigenpair
#'
#' @description Calculate the maximal eigenpair for the general matrix.
#'
#' @param A The input general matrix.
#' @param v0_tilde The unnormalized initial vector \eqn{\tilde{v0}}.
#' @param z0 The type of initial \eqn{z_0} used to calculate the approximation of \eqn{\rho(Q)}.
#'        There are three types: 'fixed', 'Auto' and 'numeric' corresponding to three choices
#'        of \eqn{z_0} in paper.
#' @param z0numeric The numerical value assigned to initial \eqn{z_0} as an approximation of
#'        \eqn{\rho(Q)} when z_0='numeric'.
#' @param improved With improved=T, the improved algorithm to calculate initial approximating
#'        eigenpairs is used.
#' @param xi The coefficient used to form the convex combination of \eqn{\delta_1^{-1}} and
#'        \eqn{(v_0,-Q*v_0)_\mu}, it should between 0 and 1.
#' @param digit.thresh The precise level of output results.
#'
#' @return A list of eigenpair object are returned, with components \eqn{z} and \eqn{v}.
#' \item{z}{The approximating sequence of the maximal eigenvalue.}
#' \item{v}{The approximating sequence of the corresponding eigenvector.}
#'
#' @seealso \code{\link{eff.ini.maxeig.tri}} for the tridiagonal matrix maximal eigenpair.

#' @examples A = matrix(c(1, 1, 3, 2, 2, 2, 3, 1, 1), 3, 3)
#' eff.ini.maxeig.general(A, v0_tilde = rep(1, dim(A)[1]), z0 = 'fixed')
#'
#' A = matrix(c(1, 1, 3, 2, 2, 2, 3, 1, 1), 3, 3)
#' eff.ini.maxeig.general(A, v0_tilde = rep(1, dim(A)[1]), z0 = 'Auto')
#'
#' ##Symmetrizing A converge to second largest eigenvalue
#' A = matrix(c(1, 3, 9, 5, 2, 14, 10, 6, 0, 11, 11, 7, 0, 0, 1, 8), 4, 4)
#' S = (t(A) + A)/2
#' N = dim(S)[1]
#' a = diag(S[-1, -N])
#' b = diag(S[-N, -1])
#' c = rep(NA, N)
#' c[1] = -diag(S)[1] - b[1]
#' c[2:(N - 1)] = -diag(S)[2:(N - 1)] - b[2:(N - 1)] - a[1:(N - 2)]
#' c[N] = -diag(S)[N] - a[N - 1]
#'
#' z0ini = eff.ini.maxeig.tri(a, b, c, xi = 7/8, improved = TRUE)$z[1]
#' eff.ini.maxeig.general(A, v0_tilde = rep(1, dim(A)[1]), z0 = 'numeric',
#' z0numeric = 28 - z0ini)

#' @export
eff.ini.maxeig.general = function(A, v0_tilde = NULL, z0 = NULL,
    z0numeric, improved = F, xi = 1, digit.thresh = 6) {
    if (xi < 0 | xi > 1)
        stop("The coefficient xi should between 0 and 1!")

    N = dim(A)[1]
    h = rep(NA, N)
    m = 0

    # check input matrix Q
    if (max(rowSums(A)) <= 1e-10) {
        Q = A
    } else {
        #print("Adjust the diagonal of matrix by Q=A-mI")
        m = max(rowSums(A))
        print(paste0("m=", max(rowSums(A))))
        Q = A - m * diag(1, N)
    }

    if ((sum(rowSums(Q)[1:(N - 1)] < 0) == 0) & (rowSums(Q)[N] <
        0)) {
        h = rep(1, N)
    } else {
        h[1] = 1
        h[2:N] = solve(Q[-N, -1], -Q[1:(N - 1), 1])
    }

    q = -diag(Q)

    P = matrix(rep(1/(q * h), N), N, N) * Q * matrix(rep(h, each = N),
        N, N) + diag(rep(1, N))

    if (-Q[N, N] - sum(Q[N, 1:(N - 1)] * h[1:(N - 1)])/h[N] < sum(Q[N,
        1:(N - 1)] * h[1:(N - 1)])/h[N]/100) {
        x = rep(1, N)
    } else {
        x = rep(NA, N)
        x[1] = 1
        x[2:N] = solve(diag(1, N - 1) - P[-1, -1], P[2:N, 1])
    }

    if (is.null(v0_tilde)) {
        v0_tilde = h * sqrt(x)
    }

    if (improved) {
        Qtilde = diag(1/h, N) %*% Q %*% diag(h, N)
        Q0tilde = Qtilde
        Q0tilde[N, N] = -sum(Q0tilde[N, 1:(N - 1)])

        mu = rep(1, N)
        mu[2:N] = solve(t(Qtilde)[-N, -1], -t(Qtilde)[1:(N - 1), 1])

        deltapart1 = mu * sqrt(x)
        deltapart1sum = cumsum(deltapart1)
        deltapart2 = mu * x^(3/2)
        deltapart2sum = rev(cumsum(rev(deltapart2)))

        delta1 = max(c(sqrt(x[1:(N - 1)]) * deltapart1sum[1:(N -
            1)] + 1/sqrt(x[1:(N - 1)]) * deltapart2sum[2:N], sqrt(x[N]) *
            deltapart1sum[N]))/(1 - x[2])

        v0 = v0_tilde/sqrt(sum(v0_tilde^2 * mu))
        zstart = xi * 1/delta1 + (1 - xi) * sum(v0 * (-Q %*% v0) *
            mu)
        ray = ray.quot(Q = Q, mu = mu, v0_tilde = v0_tilde, zstart = zstart,
            digit.thresh = digit.thresh)
    } else {
        v0 = v0_tilde/sqrt(sum(v0_tilde^2))
        zstart = switch(z0, fixed = 0, Auto = sum(v0 * (-Q %*% v0)),
            numeric = z0numeric)

        ray = ray.quot(Q = Q, mu = rep(1, N), v0_tilde = v0_tilde,
            zstart = zstart, digit.thresh = digit.thresh)
    }

    return(list(z = m - unlist(ray$z), v = ray$v))
}
