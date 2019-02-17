#
# R version of cmleDAG
# 
# @reference: Li, Shen, Pan "Likelihood ratio tests for a large directed acyclic graph"
# @author: Chunlin Li (li000007@umn.edu), Ziyue Zhu (zhux0502@umn.edu)
# @version: 0.19.02.17
# 

library(mvtnorm)
source("clrdag.R")

#
# Example 1: random graph
#

set.seed(2018)

p=250
n=200
# random graph: randomly generate adjacency matrix A, A lower triangular
sparsity = 1/p
A = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A, diag = T)] = 0
Sigma = solve(diag(p) - A)
Sigma = Sigma%*%t(Sigma)
# X ~ Normal(0,Sigma)
X = rmvnorm(n,mean=rep(0,p), sigma=Sigma, method = "chol") 

t1 = proc.time()
out = cmleDAG(X=X,tau=0.35,mu=1,rho=1.2)
Eout = ifelse(out$A!=0,1,0)
E = ifelse(A!=0,1,0)
sum(abs(E - Eout))
sum(E - Eout > 0)
proc.time() - t1



set.seed(2018)

p=250
n=200
# hub graph: randomly generate adjacency matrix A, A lower triangular
A = matrix(0,p,p)
A[,1] = sign(runif(p,min=-1,max=1))
A[1,1] = 0
Sigma = solve(diag(p) - A)
Sigma = Sigma%*%t(Sigma)
# X ~ Normal(0,Sigma)
X = rmvnorm(n,mean=rep(0,p), sigma=Sigma, method = "chol") 

t1 = proc.time()
out = cmleDAG(X=X,tau=0.35,mu=1,rho=1.2)
Eout = ifelse(out$A!=0,1,0)
E = ifelse(A!=0,1,0)
sum(abs(E - Eout))
sum(E - Eout > 0)
proc.time() - t1







