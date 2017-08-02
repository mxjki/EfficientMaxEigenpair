#' @title Solve the linear equation (-Q-zI)w=v.
#' @description Construct the solution of linear equation (-Q-zI)w=v.
#'
#' @param Q The given tridiagonal matrix.
#' @param v The column vector on the right hand of  equation.
#' @param z The Rayleigh shift.
#' @return A solution sequence \eqn{w} to the equation (-Q-zI)w=v.
#' 
#' @examples
#' a = c(1:7)^2
#' b = c(1:7)^2
#' c = rep(0, length(a) + 1)
#' c[length(a) + 1] = 8^2
#' N = length(a)
#' zstart = 6
#' Q = tridiag(b, a, -c(b[1] + c[1], a[1:N - 1] + b[2:N] + c[2:N], a[N] + c[N + 1]))
#' tri.sol(Q, z=zstart, v=rep(1,dim(Q)[1]))

#' @export
tri.sol = function(Q, z, v) {
    
    N = dim(Q)[1] - 1
    indx <- seq.int(N)
    a = Q[cbind(indx + 1, indx)]
    b = c(Q[cbind(indx, indx + 1)], -sum(Q[N + 1, ]))
    
    if (sum(b[1:N] != a) == 0) {
        mu0 = rep(1, N + 1)
    } else {
        mu0 = rep(1, N + 1)
        mu0[2:(N + 1)] = cumprod(b[-(N + 1)]/a)
    }
    
    M = matrix(0, N + 1, N + 1)
    nv = 1/(mu0 * b)
    
    for (s in 1:(N + 1)) {
        
        w0 = rep(0, N + 1)
        w0[1:s] = rev(cumsum(rev(nv[1:s])))
        M[s, ] = mu0 * w0
        
    }
    
    A = rep(1, N + 1)
    A[1] = 0
    B = rep(1, N + 1)
    B[1] = 1
    
    for (s in 2:(N + 1)) {
        A[s] = -sum(M[s - 1, ] * (v + z * A))
        B[s] = 1 - z * sum(M[s - 1, ] * B)
    }
    
    x = (sum(mu0 * (v + z * A)) - mu0[N + 1] * b[N + 1] * A[N + 1])/(mu0[N + 1] * 
        b[N + 1] * B[N + 1] - z * sum(mu0 * B))
    
    w = A + x * B
    
    return(w)
}
