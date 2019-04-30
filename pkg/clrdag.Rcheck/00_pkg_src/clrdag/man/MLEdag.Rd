\name{clrdag}
\alias{clrdag}
\alias{MLEdag}
\title{
MLEdag
}
\description{
A function computes the MLE/LRT of a Gaussian directed acyclic graph with specified constraints.
}
\usage{
MLEdag(X, A = NULL, Lambda = NULL, D = NULL, tau, mu, rho, 
        tol_abs = 1e-04, tol_rel = 1e-04, 
        dc_max_iter = 20, admm_max_iter = 1000, trace_obj = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    An n by p data matrix, where n is the number of observations and p is the dimension. 
}
  \item{A, Lambda}{
  Initial estimate. \code{A} is a p by p adjacency matrix, 
  \code{Lambda} is a p by p dual matrix in acyclicity condition. 
  \code{A} must be a DAG! 
  If \code{A} is NULL (default), 
  the initial estimate is provided automatically (Be careful!). 
}
  \item{D}{
    A p by p matrix indicating hypothesized edges. 
    For the entries equal to 1, no sparse penalty is imposed. 
}
  \item{tau}{
    A positive real number. \code{tau} is the threshold parameter in TLP.
}
  \item{mu}{
    A positive real number. \code{mu} is the sparsity parameter. 
}
  \item{rho}{
    A positive real number. \code{rho} is the ADMM dual parameter. 
}
  \item{tol_abs, tol_rel}{
    Positive real. The absolute and relative tolerance. 
}
  \item{dc_max_iter, admm_max_iter}{
    Positive integer. The maximum iteration number of DC and ADMM. 
}
%  \item{test_path}{
%    Logical. If TRUE, the path test is used. 
%}
  \item{trace_obj}{
    Logical. If TRUE, the objective values are printed after each iteration. 
}
}

\value{
  The function returns a LIST containing the following components. 
  \item{X}{The input data matrix.}
  \item{A}{The final estimate of adjacency matrix.}
  \item{Lambda}{The final estimate of dual variables in the acyclicity condition.}
  \item{mu}{The input sparsity parameter.}
  \item{tau}{The input threshold parameter in TLP.}
}
\references{
  Li, C., Shen, X., and Pan, W. (2019). Likelihood ratio tests of a large directed acyclic graph.
  Submitted. 
}
\author{
Chunlin Li
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## library(mvtnorm)
##
## Example 1: random graph
##
set.seed(2019)
p<-10
n<-1000
## random graph: randomly generate adjacency matrix A, A lower triangular
sparsity <- 2/p
A <- matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A, diag = TRUE)] <- 0
#Sigma <- solve(diag(p) - A)
#Sigma <- Sigma \%*\% t(Sigma)
#X <- rmvnorm(n,mean=rep(0,p), sigma=Sigma, method="chol") 
X <- matrix(rnorm(n*p), n, p) \%*\% t(solve(diag(p) - A))
out <- MLEdag(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE) # compute the MLE
B <- out$A
B <- ifelse(abs(B)>0.3,1,0)
all(B == abs(A))
##
## Example 2: hub graph
##
set.seed(2019)
p<-10
n<-1000
## hub graph: randomly generate adjacency matrix A, A lower triangular
A <- matrix(0,p,p)
A[,1] <- sign(runif(p,min=-1,max=1))
A[1,1] <- 0
#Sigma <- solve(diag(p) - A)
#Sigma <- Sigma \%*\% t(Sigma)
#X <- rmvnorm(n,mean=rep(0,p), sigma=Sigma, method="chol") 
X <- matrix(rnorm(n*p), n, p) \%*\% t(solve(diag(p) - A))
out <- MLEdag(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE) # compute the MLE
B <- out$A
B <- ifelse(abs(B)>0.3,1,0)
all(B == abs(A))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MLE}
\keyword{Likelihood ratio test}
\keyword{Directed acyclic graph}