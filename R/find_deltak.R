#' @title Compute \eqn{\delta_k}
#' @description Compute \eqn{\delta_k} for given vector \eqn{v} and matrix \eqn{Q}.
#'
#' @param Q The given tridiagonal matrix.
#' @param v The column vector on the right hand of equation.
#' @return A list of \eqn{\delta_k} for given vector \eqn{v} and matrix \eqn{Q}.
#' 
#' @examples
#' a = c(1:7)^2
#' b = c(1:7)^2
#' c = rep(0, length(a) + 1)
#' c[length(a) + 1] = 8^2
#' N = length(a)
#' Q = tridiag(b, a, -c(b[1] + c[1], a[1:N - 1] + b[2:N] + c[2:N], a[N] + c[N + 1]))
#' find_deltak(Q, v=rep(1,dim(Q)[1]))

#' @export
find_deltak = function(Q, v) {
    
    N = dim(Q)[1] - 1
    indx <- seq.int(N)
    a = Q[cbind(indx + 1, indx)]
    b = c(Q[cbind(indx, indx + 1)], -Q[N + 1, N + 1] - Q[N + 1, N])
    
    if (sum(b[1:N] != a) == 0) {
        mu = rep(1, N + 1)
    } else {
        mu = rep(1, N + 1)
        mu[2:(N + 1)] = cumprod(b[1:N]/a)
    }
    
    phi = rev(cumsum(rev(1/(mu * b))))
    
    deltapart1 = mu * v
    deltapart1sum = cumsum(deltapart1)
    deltapart2 = mu * phi * v
    deltapart2sum = rev(cumsum(rev(deltapart2)))
    
    delta1 = rep(NA, N + 1)
    delta1 = phi[1:N] * deltapart1sum[1:N] + deltapart2sum[2:(N + 1)]
    delta1[N + 1] = phi[N + 1] * deltapart1sum[N + 1]
    
    deltak = max(delta1/v)
    deltak_prime = min(delta1/v)
    
    return(list(deltak = deltak, deltak_prime = deltak_prime))
}
