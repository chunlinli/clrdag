## 
## reference: Li, Shen, Pan "Likelihood ratio tests for a large directed acyclic graph"
## author: Chunlin Li (li000007@umn.edu)
## version: 0.99.99: cpp functions are parallelized
## 

## depends on libraries: Rcpp, RcppArmadillo
## recommend to link R to openblas or Intel MKL

## MLEdag: estimate a gaussian directed acyclic graph with specified constraints 
## X: data matrix
## A, Lambda: initial estimate, A must be a DAG !!!
## D: index of hypothesized edges, no penalty is imposed
## tau: threshold par in TLP
## mu: sparsity par
## rho: ADMM par
## tol, rel_tol: absolute/relative tolerance
## dc_max_iter, admm_max_iter: DC/ADMM max iteration
## trace_obj: boolean value, indicate whether to print the obj function after each DC iter

MLEdag <- function(X, A = NULL, Lambda = NULL, D = NULL, tau, mu, rho, 
                    tol_abs = 1e-4, tol_rel = 1e-4, 
                    dc_max_iter = 20, admm_max_iter = 1000, trace_obj = TRUE) {
    p <- ncol(X)
  
    ## initialization of (A, Lambda, D), may need more work
    if (is.null(A)) 
    {
        auto_init <- TRUE
        A <- matrix(0,p,p)
        Lambda <- matrix(1,p,p) 
    } 
    else 
    { 
        auto_init <- FALSE
        if (!all.equal(dim(A),c(p,p))) 
        {
            stop("Invalid input A!")
        }
        if (!all.equal(dim(Lambda),c(p,p))) 
        {
            stop("Invalid input Lambda!")
        }
        diag(A) <- 0
    }
    if (is.null(D)) 
    { 
        D <- matrix(FALSE,p,p) 
    } 
    else 
    { 
        if (!all.equal(dim(D),c(p,p))) 
        {
            stop("Invalid input D!")
        } 
        diag(D) <- FALSE
    }
  
    ## call DC_ADMM cpp function
    out <- DC_ADMM(X, A, Lambda, D, tau, mu, rho, tol_abs, tol_rel, dc_max_iter, admm_max_iter, auto_init, trace_obj) 
    A <- out$A
    Lambda <- out$Lambda
    
    ## hypothesis testing 
    if (sum(D) == 0) # no test
    {
        return(list(X = X, A = A, Lambda = Lambda, tau = tau, mu = mu))
    }
    else # likelihood ratio test
    {
        n <- nrow(X)
        XTX <- crossprod(X)
        NZ <- (A != 0)
        df <- sum(NZ * D != 0) # degrees of freedom 
        A0 <- matrix(0, p, p)
        for (i in 1:p)
        {
            idx <- which(A[i,] != 0 & D[i,] == 0)
            if (length(idx) > 0)
            {
                A0[i,idx] <- solve(XTX[idx,idx]) %*% XTX[idx,i]
            }
        }
        lrt <- max((p*n - sum(A))*(1 - sum((X %*% (diag(p) - t(A)))^2)/sum((X %*% (diag(p) - t(A0)))^2)),0)
        pval <- ifelse(df < 30, pchisq(lrt, df, lower.tail = FALSE),
                                pnorm((lrt - df) / sqrt(2*df), lower.tail = FALSE))
        return(list(X = X, A = A, Lambda = Lambda, D = D, tau = tau, mu = mu, lrt = lrt, df = df, pval = pval))
    }
}
