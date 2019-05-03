## Installation 

This package requires compilation to install from source. The following instructions are for new R users. 

### Windows users: 

0. Make sure you have [R](https://www.r-project.org/) (>= 3.5.0) and [Rtools35](https://cran.r-project.org/bin/windows/Rtools/) or later properly installed. In particular, add R and Rtools to your PATH variable. In R, install the packages ```Rcpp``` and ```RcppArmadillo```.

1. Download ```clrdag.zip``` from ```pkg``` directory. Unzip ```clrdag.zip```.

2. In Command Prompt, 
```
cd "path\\to\\clrdag.zip"
R CMD build clrdag
R CMD INSTALL clrdag_*.tar.gz
```
3. Delete ```clrdag``` folder, ```clrdag.zip```, and ```clrdag_*.tar.gz```.

### \*nix users:

0. Make sure you have [R](https://www.r-project.org/) (>= 3.5.0) and building tools installed. 
In R, install the packages ```Rcpp``` and ```RcppArmadillo```.

1. Download ```clrdag.zip``` from ```pkg``` directory. Unzip ```clrdag.zip```.

2. In a terminal app (Bash), 
```
cd "path/to/clrdag.zip"
R CMD build clrdag
R CMD INSTALL clrdag_*.tar.gz
```
3. Remove ```clrdag``` directory, ```clrdag.zip```, and ```clrdag_*.tar.gz```.

## Example(s)

Try the following example in R!

```
library(clrdag)
set.seed(2018)
p <- 50; n <- 1000; sparsity <- 2/p
A <- matrix(rbinom(p*p, 1, sparsity)*sign(runif(p*p, -1, 1))*runif(p*p, 0.7, 1), p, p)
A[upper.tri(A, diag=TRUE)] <- 0
idx <- sample(1:p)
A <- A[idx,idx]
X <- matrix(rnorm(n*p), n, p) %*% t(solve(diag(p) - A))

t <- proc.time()
out <- MLEdag(X=X, tau=0.3, mu=1, rho=1.2)
proc.time() - t

sum((out$A - A)^2)
sum(abs((out$A != 0) - (A != 0)))
```

## To-dos in future versions

Cross-validation functions, diagnostic plots. 

## About ```clrdag``` 

This is an R package for likelihood estimation and inference of a Gaussian directed acyclic graph. The reference for the details is 

Li, C., Shen, X., and Pan, W. (2019). Likelihood ratio tests of a large directed acyclic graph. Submitted. 

The program was originally written in July 2018. Ziyue Zhu (thanks to her!) and I refactored the R code version in December 2018 as part of the final project for the optimization course EE 5239. As of May 1, 2019, this R package is publicly available. It is still in development and be cautious.

To report an issue, please file it at [here](https://github.com/chunlinli/clrdag/issues).
