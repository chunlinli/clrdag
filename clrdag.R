## 
## reference: Li, Shen, Pan "Likelihood ratio tests for a large directed acyclic graph"
## author: Chunlin Li (li000007@umn.edu), Ziyue Zhu (zhux0502@umn.edu)
## version: 2019.01: cpp functions are updated
## 

library(Rcpp)
library(RcppArmadillo)
sourceCpp("dcadmm.cc")

## cmleDAG: estimate a gaussian directed acyclic graph with specified constraints 
## X: data matrix
## A, Lambda: initial estimate, A must be a DAG !!!
## D: index of hypothesized edges, no penalty is imposed
## tau: threshold par in TLP
## mu: sparsity par
## rho: ADMM par
## tol, rel_tol: absolute/relative tolerance
## dc_max_iter, admm_max_iter: DC/ADMM max iteration
## trace_obj: boolean value, indicate whether to print the obj function after each DC iter

cmleDAG = function(X, A=NULL, Lambda=NULL, D=NULL, tau, mu, rho, tol_abs=1e-4, 
                  tol_rel=1e-4, dc_max_iter=20, admm_max_iter=500, trace_obj=TRUE) {
  p = ncol(X)
  
  ## initialization of (A, Lambda, D)
  if (is.null(A)) { A = matrix(0,p,p); Lambda = matrix(1,p,p) } 
  else { 
    if(!all.equal(dim(A0),c(p,p))) stop("Invalid input A!") 
    if(!all.equal(dim(lambda0),c(p,p))) stop("Invalid input Lambda!")
  }
  if(is.null(D)) { D = matrix(FALSE,p,p) }
  else { if(!all.equal(dim(D),c(p,p))) stop("Invalid input D!") }
  
  ## call DCADMM function
  out = DC_ADMM(X,A,Lambda,D,tau,mu,rho,tol_abs,tol_rel,dc_max_iter,admm_max_iter,trace_obj) 
  A = out$A
  Lambda = out$Lambda
  
  ## return values
  return(list(X=X, A=A, Lambda=Lambda, tau=tau, mu=mu))
}
