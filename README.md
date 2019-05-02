## Installation 

This package requires compilation to install from source. The following instructions are for new R users. 

### Windows users: 

0. Make sure you have R (>= 3.5.0) and Rtools35 or later (available [here](https://cran.r-project.org/bin/windows/Rtools/)) properly installed. In particular, add R and Rtools to your PATH variable. In R, install the packages ```Rcpp``` and ```RcppArmadillo```.

1. Download ```clrdag_*.tar.gz``` from ```pkg``` directory.

2. In Command Prompt, 
```
cd "path\\to\\downloaded\\package\\" && R CMD INSTALL clrdag_*.tar.gz
```
3. Delete ```clrdag_*.tar.gz```.

### Unix(-like) users:

0. Make sure you have R (>= 3.5.0) and building tools installed. 
In R, install the packages ```Rcpp``` and ```RcppArmadillo```.

1. Download ```clrdag_*.tar.gz``` from ```pkg``` directory.

2. In terminal (Bash), 
```
cd "path/to/downloaded/package/" && R CMD INSTALL clrdag_*.tar.gz
```
3. Remove ```clrdag_*.tar.gz```.

## Examples

Try the following example in R!

```
library(clrdag)
set.seed(2018)
p <- 50; n <- 1000; sparsity <- 2/p
A <- matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1))*runif(p*p,min = 0.7, max = 1),p,p)
A[upper.tri(A, diag = TRUE)] <- 0
idx <- sample(1:p)
A <- A[idx,idx]
X <- matrix(rnorm(n*p), n, p) %*% t(solve(diag(p) - A))

t <- proc.time()
out <- MLEdag(X=X,tau=0.3,mu=1,rho=1.2)
proc.time() - t

sum((out$A - A)^2)
sum(abs((out$A != 0) - (A!=0)))
```

## To-dos in future versions

Cross-validation functions, Diagnostic plots. 

## About ```clrdag``` 

This is a C++/R implementation of the constrained likelihood ratio tests of a DAG. The reference for the details is 

Li, C., Shen, X., and Pan, W. (2019). Likelihood ratio tests of a large directed acyclic graph. Submitted. 

The program was originally written in July 2018. Then Ziyue Zhu (thanks to her!) and I refactored the R code version (available at [here](https://github.umn.edu/li000007/clrdag_r/)) in December 2018 as part of the final project for the optimization course EE 5239. As of April 2019, this R package can be accessed at UMN. It still in development and be cautious.

To report an issue, please file an issue at [here](https://github.com/chunlinli/clrdag/issues).
