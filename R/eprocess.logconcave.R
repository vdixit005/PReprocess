#' Calculate PReprocess for testing logconcavity
#'
#' This function takes in data and gives the user an e-value for testing logconcavity
#'
#' @param X Univariate data as a vector
#' @param d Parametric kernel function of form d(x,U)
#' @param U Mixing density support, where U is an equispaced and odd length vector as the support for univariate f
#' @param f A vector of initial guess for f or final f based on old data
#' @param n Total sample size already used (default is 0)
#' @param L Value of PR log-likelihood based on old data (default is 0)
#' @return A scalar log(e-value) based on PR, PR object f_n and PR log-likelihood L (to use if further data gets incorporated)
#' @examples
#' # Generate data from a bimodal density (not log-concave)
#' datagen = function(n){
#' ans = numeric(n)
#' u = runif(n)
#' for(i in 1:n){
#'  if(u[i]<0.75){
#'    ans[i] = rnorm(1,0,sqrt(2))
#'  } else {
#'    ans[i] = rnorm(1,10,sqrt(2))
#'  }
#' }
#' return(ans)
#' }
#' X = datagen(1000)
#' # PR initialization
#' U = cbind(seq(-10, 20, length.out = 101), seq(0.01, 3, length.out = 101))
#' f0 = rep(1, 101*101)
#' d1 = function(x,u){mapply(function(mean,sd) dnorm(x,mean,sd), mean=u[,1], sd=u[,2])}
#' ans = eprocess.logconcave(X,d=d1,U=U,f=f0,n=0,L=0)
#' # If new data is received, update the process
#' Y = datagen(1000)
#' X = c(X, Y)
#' ans2 = eprocess.logconcave(X,d=d1,U=U,f=ans$f_n$f,n=1000, L=ans$`PR-loglik`)
#' # More data is received
#' Y = datagen(1000)
#' X = c(X, Y)
#' ans3 = eprocess.logconcave(X,d=d1,U=U,f=ans2$f_n$f,n=2000, L = ans2$`PR-loglik`)
#'
#' @source The PR functions are based on code from \url{https://www4.stat.ncsu.edu/~rmartin/Codes/pr.R}
#' @references
#' Dixit, Vaidehi, and Ryan Martin. "Anytime valid and asymptotically optimal statistical inference
#' driven by predictive recursion" arXiv preprint arXiv:2309.13441 (2023)
#' ########################################################################################
#' @export
eprocess.logconcave = function(X,d,U,f,n=0,L=0){
  w <- function(i) {1 / (i + n + 1)^0.67}
  # PR
  f.pr = pr(X = X[(n+1):length(X)], d = d, U = U, f0 = f, w = w)
  L = L + f.pr$L
  # Get log-concave MLE on X
  npmle <- logcondens::logConDens(x = X)

  # Evaluate log-concave MLE on X1
  eval_X <- logcondens::evaluateLogConDens(
    xs = X, res = npmle, which = 2)[, 3]

  # log-PRevalue
  log_evalue = (L - sum(log(eval_X)))

  return(list("logeprocess" = log_evalue, "f_n" = f.pr, "PR-loglik" = L))
}

int <- function(f, x, tol=1e-10) {

  n <- length(x)
  if(n == 1) return(0)
  simp <- (sum(range(x[-1] - x[-n]) * c(-1, 1)) < tol)
  if(!simp) out <- sum((f[-1] + f[-n]) * (x[-1] - x[-n])) / 2
  else out <- ((x[2] - x[1]) / 3) * (sum(f) + sum(f[2:(n-1)]) + 2 * sum(f[seq(2, n-1, 2)]) )
  return(out)

}

simp.int2 <- function(x, y, Fxy, tol=1e-10) {

  n <- length(x); if(n %% 2 == 0) stop("need an odd number of grid points")
  m <- length(y); if(m %% 2 == 0) stop("need an odd number of grid points")

  check <- (sum(range(x[-1]-x[-n])*c(-1,1)) < tol) && (sum(range(y[-1]-y[-m])*c(-1,1)) < tol)
  if(!check) stop("xy-grid must be equi-spaced")

  h <- x[2] - x[1]  # spacings for x (assumed to be equi-spaced)
  k <- y[2] - y[1]  # spacings for y

  u <- c(1, rep(c(4,2), (n-3)/2), 4, 1)
  v <- c(1, rep(c(4,2), (m-3)/2), 4, 1)
  M <- outer(u, v, '*')

  return(h * k * sum(M * Fxy) / 9)

}

pr <- function(X, d, U, f0, w, nperm = 1, perm = NULL,...) {
  X = as.matrix(X)
  U = as.matrix(U)
  n <- nrow(X)
  t <- nrow(U)
  du <- ncol(U)
  if(missing(f0)) f0 <- 1 + 0 * U[,1]
  # if(missing(w)) w <- function(i) 1 / (i + 1)
  N <- nperm
  if(N == 1 && is.null(perm)) {
    perm <- matrix(0, n, N)
    perm[,1] <- 1:n}
  else if(N > 1 && is.null(perm)) {
    perm <- matrix(0, n, N)
    perm[,1] <- 1:n
    for(j in 2:N) perm[,j] <- sample(n)

  }
  f.avg = 0 * f0
  L.avg = 0
  if(du==1){
    f0 = f0 / int(f0, U)
    for(j in 1:N) {
      f = f0
      L = 0
      x <- X[perm[,j],]
      for(i in 1:n) {
        num <- d(x[i], U,...) * f
        den <- int(num, U)
        L <- L + log(den)
        f <- (1 - w(i)) * f + w(i) * num / den
      }
      f.avg <- (j - 1) * f.avg / j + f / j
      L.avg <- (j - 1) * L.avg / j + L / j
    }
    ans = list(U = U, f=f.avg, L=L.avg)
    class(ans) = "pr"
    return(ans)
  }
  else if(du==2){
    f0 = f0 / simp.int2(U[,2], U[,1], matrix(f0, t, t, byrow=TRUE))
    U.l = as.matrix(expand.grid(U[,1],U[,2]))
    for(j in 1:N) {
      f = f0
      L = 0
      x <- as.matrix(X[perm[,j], ])
      for(i in 1:n){
        num = d(x[i,], U.l,...) * f
        f_matrix = matrix(num, nrow = t, ncol = t, byrow=TRUE)
        den = simp.int2(U[,2], U[,1], f_matrix)
        L <- L + log(den)
        f <- (1 - w(i)) * f + w(i) * num / den
      }
      f.avg <- (j - 1) * f.avg / j + f / j
      L.avg <- (j - 1) * L.avg / j + L / j
    }
    ans = list(U = U, f=f.avg, L=L.avg)
    class(ans) = "pr"
    return(ans)
  }
  else {
    w.vec = w(1:n)
    cpp.final = pr_cpp(f0, U, X, d, N, w.vec)
    ans = list(U = U, f=cpp.final$f, L=-cpp.final$L, D = cpp.final$D)
    class(ans) = "pr"
    return(ans)
  }
}
