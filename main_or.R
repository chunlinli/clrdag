# 
# @reference: Li, Shen, Pan "Likelihood ratio tests for a large directed acyclic graph"
# @author: Chunlin Li (li000007@umn.edu), Ziyue Zhu (zhux0502@umn.edu)
# @version: 0.18.12
# 
library(mvtnorm)
source("clrdag.R")

#
# Example 1: random graph
#
set.seed(2019)
p=50
n=1000
# random graph: randomly generate adjacency matrix A, A lower triangular
sparsity = 2/p
A = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A, diag = T)] = 0
Sigma = solve(diag(p) - A)
Sigma = Sigma%*%t(Sigma)
OR_random = rep(0,100)
for(ii in 1:100)
{
  # X ~ Normal(0,Sigma)
  X = rmvnorm(n,mean=rep(0,p), sigma=Sigma, method = "chol") 
  out = cmleDAG(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE)
  B = out$A
  B = ifelse(abs(B)>0.3,1,0)
  if(sum(abs(B-abs(A)))==0) OR_random[ii] = 1
}
mean(OR_random) # 0.97



#
# Example 2: hub graph
#
set.seed(2019)
p=50
n=1000
# hub graph: randomly generate adjacency matrix A, A lower triangular
A = matrix(0,p,p)
A[,1] = sign(runif(p,min=-1,max=1))
A[1,1] = 0
Sigma = solve(diag(p) - A)
Sigma = Sigma%*%t(Sigma)
OR_hub = rep(0,100)
for(ii in 1:100)
{
  # X ~ Normal(0,Sigma)
  X = rmvnorm(n,mean=rep(0,p), sigma=Sigma, method = "chol") 
  out = cmleDAG(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE)
  B = out$A
  B = ifelse(abs(B)>0.3,1,0)
  if(sum(abs(B-abs(A)))==0) OR_hub[ii] = 1
}
mean(OR_hub) # 1.00





#
# Example 3: random graph
#
set.seed(2019)
p=30
n=600
# random graph: randomly generate adjacency matrix A, A lower triangular
sparsity = 2/p
A = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A, diag = T)] = 0
Sigma = solve(diag(p) - A)
Sigma = Sigma%*%t(Sigma)
OR_random = rep(0,100)
for(ii in 1:100)
{
  # X ~ Normal(0,Sigma)
  X = rmvnorm(n,mean=rep(0,p), sigma=Sigma, method = "chol") 
  out = cmleDAG(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE)
  B = out$A
  B = ifelse(abs(B)>0.3,1,0)
  if(sum(abs(B-abs(A)))==0) OR_random[ii] = 1
}
mean(OR_random) # 1.00






#
# Example 4: hub graph
#
set.seed(2019)
p=30
n=600
# hub graph: randomly generate adjacency matrix A, A lower triangular
A = matrix(0,p,p)
A[,1] = sign(runif(p,min=-1,max=1))
A[1,1] = 0
Sigma = solve(diag(p) - A)
Sigma = Sigma%*%t(Sigma)
OR_hub = rep(0,100)
for(ii in 1:100)
{
  # X ~ Normal(0,Sigma)
  X = rmvnorm(n,mean=rep(0,p), sigma=Sigma, method = "chol") 
  out = cmleDAG(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE)
  B = out$A
  B = ifelse(abs(B)>0.3,1,0)
  if(sum(abs(B-abs(A)))==0) OR_hub[ii] = 1
}
mean(OR_hub) # 1.00






#
# Example 5: random graph
#
set.seed(2019)
p=10
n=200
# random graph: randomly generate adjacency matrix A, A lower triangular
sparsity = 2/p
A = matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A, diag = T)] = 0
Sigma = solve(diag(p) - A)
Sigma = Sigma%*%t(Sigma)
OR_random = rep(0,100)
for(ii in 1:100)
{
  # X ~ Normal(0,Sigma)
  X = rmvnorm(n,mean=rep(0,p), sigma=Sigma, method = "chol") 
  out = cmleDAG(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE)
  B = out$A
  B = ifelse(abs(B)>0.3,1,0)
  if(sum(abs(B-abs(A)))==0) OR_random[ii] = 1
}
mean(OR_random) # 0.99



#
# Example 6: hub graph
#
set.seed(2019)
p=10
n=200
# hub graph: randomly generate adjacency matrix A, A lower triangular
A = matrix(0,p,p)
A[,1] = sign(runif(p,min=-1,max=1))
A[1,1] = 0
Sigma = solve(diag(p) - A)
Sigma = Sigma%*%t(Sigma)
OR_hub = rep(0,100)
for(ii in 1:100)
{
  # X ~ Normal(0,Sigma)
  X = rmvnorm(n,mean=rep(0,p), sigma=Sigma, method = "chol") 
  out = cmleDAG(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE)
  B = out$A
  B = ifelse(abs(B)>0.3,1,0)
  if(sum(abs(B-abs(A)))==0) OR_hub[ii] = 1
}
mean(OR_hub) # 1.00


