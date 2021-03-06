#' @title Shifted inverse iteration algorithm for Tridiagonal matrix
#' @description Shifted inverse iteration algorithm algorithm to computing the maximal eigenpair of
#' tridiagonal matrix \eqn{Q}.
#'
#' @param Q The input matrix to find the maximal eigenpair.
#' @param mu A vector.
#' @param v0_tilde The unnormalized initial vector \eqn{\tilde{v0}}.
#' @param zstart The initial \eqn{z_0} as an approximation of \eqn{\rho(Q)}.
#' @param digit.thresh The precise level of output results.
#' @return A list of eigenpair object are returned, with components \eqn{z}, \eqn{v} and \eqn{iter}.
#' \item{z}{The approximating sequence of the maximal eigenvalue.}
#' \item{v}{The approximating eigenfunction of the corresponding eigenvector.}
#' \item{iter}{The number of iterations.}
#'
#' @examples
#' a = c(1:7)^2
#' b = c(1:7)^2
#' c = rep(0, length(a) + 1)
#' c[length(a) + 1] = 8^2
#' N = length(a)
#' Q = tridiag(b, a, -c(b[1] + c[1], a[1:N - 1] + b[2:N] + c[2:N], a[N] + c[N + 1]))
#' shift.inv.tri(Q, mu=rep(1,dim(Q)[1]), v0_tilde=rep(1,dim(Q)[1]), zstart=6,
#'  digit.thresh = 6)

#' @export
shift.inv.tri = function(Q, mu, v0_tilde, zstart, digit.thresh = 6) {
    z = list()
    rz = list()
    v = list()
    w = list()
    deltak = list()
    
    v[[1]] = v0_tilde/sqrt(sum(v0_tilde^2 * mu))
    
    z[[1]] = zstart
    rz[[1]] = round(zstart, digit.thresh)
    
    ratio_z = 1
    iter = 0
    
    while (ratio_z >= 10^(-digit.thresh)) {
        iter = iter + 1
        
        wk = tri.sol(Q, z[[iter]], v[[iter]])
        
        if (sum(sign(wk)) == length(wk) | sum(sign(wk)) == -length(wk)) {
            w = append(w, list(wk))
        } else {
            w = append(w, list(tri.sol(Q, z[[1]], v[[iter]])))
        }
        
        v = append(v, list(w[[iter]]/sqrt(sum(w[[iter]]^2 * mu))))
        
        delta1 = find_deltak(Q, v[[iter + 1]])
        
        deltak = append(deltak, list(delta1$deltak))
        
        z = append(z, list(1/deltak[[iter]]))
        
        ratio_z = abs(round(z[[iter + 1]], digit.thresh) - round(z[[iter]], digit.thresh))
        
        rz[[iter + 1]] = round(z[[iter + 1]], digit.thresh)
        
    }
    
    if (ratio_z == 0) {
        v = v[-(iter + 1)]
        rz = rz[-(iter + 1)]
        
        iter = iter - 1
    }
    
    return(list(v = v[iter + 1], z = rz, iter = iter))
}
