\name{MLEdag}
\alias{clrdag}
\alias{MLEdag}
\title{
MLE/LRT of a Gaussian directed acyclic graph
}
\description{
Computes the MLE/LRT of a Gaussian directed acyclic graph using difference convex programming and alternating direction method of multipliers.
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
    For the entries equal to one, no sparse penalty is imposed. 
    If \code{D} is not provided, or if all off-diagonal entries of \code{D} are zero, 
    no test is performed. 
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
  \item{A}{The final estimate of adjacency matrix. 
            Returned if no test is performed.}
  \item{A.H1, A.H0}{The final estimates of adjacency matrix under alternative and null.
                    Returned if a test is performed.}
  \item{Lambda}{The final estimate of dual variables in the acyclicity condition.
                Returned if no test is performed.}
  \item{Lambda.H1}{The final estimate of dual variables in the acyclicity condition under alternative.
                Returned if a test is performed.}
  \item{D}{A matrix indicating hypothesized edges. Returned if a test if performed.}
  \item{tau}{The input threshold parameter in TLP.} 
  \item{mu}{The input sparsity parameter.}
  \item{lrt}{(2\eqn{\times}log-likelihood ratio) of alternative over null.
             Returned if a test is performed.}
  \item{df}{Degrees of freedom of the test. Returned if a test is performed.}
  \item{pval}{The \eqn{p}-value of the test. Returned if a test is performed.}
}

\references{
  Li, C., Shen, X., and Pan, W. (2019). Likelihood ratio tests for a large directed acyclic graph.
  \emph{Journal of the American Statistical Association}. Accepted. <doi:10.1080/01621459.2019.1623042>.
}

\author{
Chunlin Li <li000007@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##
## Example: random graph
##
library(clrdag)
set.seed(2019)
p<-10
n<-1000
## random graph: randomly generate adjacency matrix A, A lower triangular
sparsity <- 2/p
A <- matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A,diag=TRUE)] <- 0
X <- matrix(rnorm(n*p),n,p) \%*\% t(solve(diag(p)-A))
out <- MLEdag(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE) # compute the MLE
sum(abs((out$A!=0)-(A!=0))) # Hamming distance to the truth graph

# test edge 1 --> 2
D <- matrix(0,p,p)
D[2,1] <- 1
out <- MLEdag(X=X,D=D,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE) # compute the MLE
out$pval
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Directed acyclic graph}
\keyword{Likelihood ratio test}
\keyword{Maximum likelihood estimation}
