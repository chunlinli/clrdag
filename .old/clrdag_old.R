## 
## reference: Li, Shen, Pan "Likelihood ratio tests for a large directed acyclic graph"
## author: Chunlin Li (li000007@umn.edu), Ziyue Zhu (zhux0502@umn.edu)
## version: 0.19.02.09: cpp functions are updated, auto initialization added
## 

## depends on libraries: Rcpp, RcppArmadillo
## recommend to link R to openblas 
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

cmleDAG = function(X, A=NULL, Lambda=NULL, D=NULL, tau, mu, rho, tol_abs=1e-4, tol_rel=1e-4, 
                    dc_max_iter=20, admm_max_iter=1000, test_path=FALSE, trace_obj=TRUE) {
    p = ncol(X)
  
    ## initialization of (A, Lambda, D), may need more work
    if (is.null(A)) 
    {
        auto_init = TRUE
        A = matrix(0,p,p)
        Lambda = matrix(1,p,p) 
    } 
    else 
    { 
        auto_init = FALSE
        if (!all.equal(dim(A),c(p,p))) 
        {
            stop("Invalid input A!")
        }
        if (!all.equal(dim(Lambda),c(p,p))) 
        {
            stop("Invalid input Lambda!")
        }
    }
    if (is.null(D)) 
    { 
        test_link = FALSE
        test_path = FALSE
        D = matrix(FALSE,p,p) 
    } 
    else 
    { 
        test_link = ifelse(test_path, FALSE, TRUE)
        if (!all.equal(dim(D),c(p,p))) 
        {
            stop("Invalid input D!")
        } 
    }
  
    ## call DC_ADMM cpp function
    out = DC_ADMM(X, A, Lambda, D, tau, mu, rho, tol_abs, tol_rel, dc_max_iter, admm_max_iter, auto_init, trace_obj) 
    A = out$A
    Lambda = out$Lambda
    
    ## testing




    ## return values
    return(list(X=X, A=A, Lambda=Lambda, tau=tau, mu=mu))
}
