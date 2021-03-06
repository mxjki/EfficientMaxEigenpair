#' @title Tridiagonal matrix maximal eigenpair
#' @description Calculate the maximal eigenpair for the tridiagonal matrix by 
#' shifted inverse iteration algorithm.
#'
#' @param a The lower diagonal vector.
#' @param b The upper diagonal vector.
#' @param c The shifted main diagonal vector. The corresponding unshift diagonal vector
#'        is -c(b[1] + c[1], a[1:N - 1] + b[2:N] + c[2:N], a[N] + c[N + 1]) where N+1
#'        is the dimension of matrix.
#' @param xi The coefficient used to form the convex combination of \eqn{\delta_1^{-1}} and
#'        \eqn{(v_0,-Q*v_0)_\mu}, it should between 0 and 1.
#' @param digit.thresh The precise level of output results.
#'
#' @return A list of eigenpair object are returned, with components \eqn{z}, \eqn{v} and \eqn{iter}.
#' \item{z}{The approximating sequence of the maximal eigenvalue.}
#' \item{v}{The approximating eigenfunction of the corresponding eigenvector.}
#' \item{iter}{The number of iterations.}
#' 
#' @seealso \code{\link{eff.ini.maxeig.tri}} for the tridiagonal matrix maximal 
#' eigenpair by rayleigh quotient iteration algorithm.
#' \code{\link{eff.ini.maxeig.general}} for the general matrix maximal eigenpair.
#'
#' @examples
#' a = c(1:7)^2
#' b = c(1:7)^2
#' c = rep(0, length(a) + 1)
#' c[length(a) + 1] = 8^2
#' eff.ini.maxeig.shift.inv.tri(a, b, c, xi = 1)

#' @export
eff.ini.maxeig.shift.inv.tri = function(a, b, c, xi = 1, digit.thresh = 6) {
    N = length(a)
    Q = tridiag(b, a, -c(b[1] + c[1], a[1:N - 1] + b[2:N] + c[2:N], a[N] + c[N + 
        1]))
    m = 0
    
    # check input vectors
    if (sum(a <= 0) > 0) 
        stop("Input vector a should be all positive!")
    if (sum(b <= 0) > 0) 
        stop("Input vector b should be all positive!")
    if (sum(c < 0) > 0) {
        print("Input vector c should be all nonnegative!")
        # print('Adjust the diagonal of matrix by Q=A-mI')
        m = max(-c)
        print(paste0("m=", m))
        c = c + m
        
        Q = tridiag(b, a, -c(b[1] + c[1], a[1:N - 1] + b[2:N] + c[2:N], a[N] + c[N + 
            1]))
    }
    if (sum(c != 0) == 0) 
        stop("Input vector c cannot be all 0!")
    if (xi < 0 | xi > 1) 
        stop("The coefficient xi should between 0 and 1!")
    
    mu = rep(1, N + 1)
    mu[2:(N + 1)] = cumprod(b/a)
    
    phi = rep(NA, N + 1)
    r = rep(NA, N)
    h = rep(NA, N + 2)
    
    if (sum(c[1:N] != 0) == 0) {
        
        b[N + 1] = c[N + 1]
        
        # phi = rev(cumsum(rev(1/(mu * b))))
        
        r = rep(1, N)
        
        h = rep(1, N + 2)
        h[N + 2] = c[N + 1]
        
    } else {
        r[1] = 1 + c[1]/b[1]
        
        for (i in 2:N) {
            r[i] = 1 + (a[i - 1] + c[i])/b[i] - a[i - 1]/(b[i] * r[i - 1])
        }
        
        h[1] = 1
        h[2:(N + 1)] = cumprod(r)
        h[N + 2] = c[N + 1] * h[N + 1] + a[N] * (h[N + 1] - h[N])
        
        b[N + 1] = 1
        # phi = rev(cumsum(rev(1/(h[1:(N + 1)] * h[2:(N + 2)] * mu * b))))
    }
    
    Q_tilde = diag(1/h[1:(N + 1)], N + 1) %*% Q %*% diag(h[1:(N + 1)], N + 1)
    mu_tilde = h[1:(N + 1)]^2 * mu
    b_tilde = c(Q_tilde[cbind(1:N, 2:(N + 1))], -Q_tilde[N + 1, N + 1] - Q_tilde[N + 
        1, N])
    phi_tilde = rev(cumsum(rev(1/(mu_tilde * b_tilde))))
    
    # deltapart1 = mu * h[-(N + 2)]^2 * sqrt(phi) deltapart1sum = cumsum(deltapart1)
    # deltapart2 = mu * h[-(N + 2)]^2 * phi^(3/2) deltapart2sum =
    # rev(cumsum(rev(deltapart2)))
    
    # delta1 = max(c(sqrt(phi[1:N]) * deltapart1sum[1:N] + 1/sqrt(phi[1:N]) *
    # deltapart2sum[2:(N + 1)], sqrt(phi[N + 1]) * deltapart1sum[N + 1]))
    
    deltapart1 = mu_tilde * sqrt(phi_tilde)
    deltapart1sum = cumsum(deltapart1)
    deltapart2 = mu_tilde * phi_tilde^(3/2)
    deltapart2sum = rev(cumsum(rev(deltapart2)))
    
    delta1 = max(c(sqrt(phi_tilde[1:N]) * deltapart1sum[1:N] + 1/sqrt(phi_tilde[1:N]) * 
        deltapart2sum[2:(N + 1)], sqrt(phi_tilde[N + 1]) * deltapart1sum[N + 1]))
    
    # v0_tilde = h[-(N + 2)] * sqrt(phi)
    
    v0_tilde = sqrt(phi_tilde)
    
    if (sum(b[1:N] != a) != 0) {
        v0 = v0_tilde/sqrt(sum(v0_tilde^2 * mu_tilde))
        zstart = xi * 1/delta1 + (1 - xi) * sum(v0 * (-Q_tilde %*% v0) * mu_tilde)
        ray = shift.inv.tri(Q = Q_tilde, mu = mu_tilde, v0_tilde = v0_tilde, zstart = zstart, 
            digit.thresh = digit.thresh)
    } else {
        v0 = v0_tilde/sqrt(sum(v0_tilde^2))
        zstart = xi * 1/delta1 + (1 - xi) * sum(v0 * (-Q_tilde %*% v0))
        ray = shift.inv.tri(Q = Q_tilde, mu = rep(1, (N + 1)), v0_tilde = v0_tilde, 
            zstart = zstart, digit.thresh = digit.thresh)
    }
    
    return(list(z = round(m - unlist(ray$z), digit.thresh), v = ray$v, iter = ray$iter))
}
