library(clrdag)
set.seed(2018)
p=250
n=200
# random graph: randomly generate adjacency matrix A, A lower triangular
sparsity = 1/p
A = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A, diag = T)] = 0
# X ~ Normal(0,Sigma)
X <- matrix(rnorm(n*p), n, p) %*% t(solve(diag(p) - A))
t1 = proc.time()
out <- MLEdag(X=X,tau=0.35,mu=1,rho=1.2,trace_obj=FALSE) # compute the MLE
B <- out$A
B <- ifelse(abs(B)>0.35,1,0)
proc.time() - t1