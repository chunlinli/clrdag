pkgname <- "clrdag"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('clrdag')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cmleDAG")
### * cmleDAG

flush(stderr()); flush(stdout())

### Name: clrdag
### Title: cmleDAG
### Aliases: clrdag cmleDAG
### Keywords: MLE Likelihood ratio test Directed acyclic graph

### ** Examples

library(mvtnorm)
##
## Example 1: random graph
##
set.seed(2019)
p<-50
n<-1000
## random graph: randomly generate adjacency matrix A, A lower triangular
sparsity <- 2/p
A <- matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A, diag = TRUE)] <- 0
Sigma <- solve(diag(p) - A)
Sigma <- Sigma %*% t(Sigma)
X <- rmvnorm(n,mean=rep(0,p), sigma=Sigma, method="chol") 
out <- cmleDAG(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE) # compute the MLE
B <- out$A
B <- ifelse(abs(B)>0.3,1,0)
B == abs(A)
##
## Example 2: hub graph
##
set.seed(2019)
p<-50
n<-1000
## hub graph: randomly generate adjacency matrix A, A lower triangular
A <- matrix(0,p,p)
A[,1] <- sign(runif(p,min=-1,max=1))
A[1,1] <- 0
Sigma <- solve(diag(p) - A)
Sigma <- Sigma %*% t(Sigma)
X <- rmvnorm(n,mean=rep(0,p), sigma=Sigma, method="chol") 
out <- cmleDAG(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE) # compute the MLE
B <- out$A
B <- ifelse(abs(B)>0.3,1,0)
B == abs(A)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
