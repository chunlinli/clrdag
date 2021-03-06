library(clrdag)
set.seed(2018)
p <- 250
n <- 200
A <- matrix(0,p,p)
A[,1] <- sign(runif(p,min=-1,max=1))
A[1,1] <- 0
D <- matrix(0, p, p)
D[p,2] = 1
pval <- rep(NA, 1000)
for(i in 1:1000) {
    X <- matrix(rnorm(n*p), n, p) %*% t(solve(diag(p) - A))
    out <- MLEdag(X=X,D=D,tau=0.35,mu=1,rho=1.2,trace_obj=FALSE)
    pval[i] <- out$pval
    print(out$pval)
}