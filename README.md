## Introduction

This is an R package for likelihood estimation and inference of a Gaussian directed acyclic graph. 
See [```clrdag-slides.pdf```](https://github.com/chunlinli/clrdag/blob/master/clrdag-slides.pdf) for the introduction slides to the package. 
The complete reference for the methodology is 

Li, C., Shen, X., and Pan, W. (2019). Likelihood ratio tests of a large directed acyclic graph. Submitted. 

See also the package vignette by ```vignette('clrdag')``` in R.

## Installation 

This package requires compilation to install from source. The following instructions are for new R users. To achieve best performance, you are recommended to install the latest version of [OpenBLAS](https://github.com/xianyi/OpenBLAS) or [Intel Math Kernel Library](https://software.intel.com/mkl). 

### Windows users: 

0. Make sure you have [R](https://www.r-project.org/) (>= 3.5.0) and [Rtools35](https://cran.r-project.org/bin/windows/Rtools/) or later properly installed. In particular, add R and Rtools to your PATH variable. In R, install the packages ```Rcpp``` and ```RcppArmadillo```.

1. Download [```clrdag.zip```](https://github.com/chunlinli/clrdag/blob/master/pkg/clrdag.zip) to local ```Downloads``` folder. Unzip ```clrdag.zip```.

2. In Command Prompt, 
```
cd "path\\to\\downloads"
R CMD build clrdag
R CMD INSTALL clrdag_*.tar.gz
```
3. Delete ```clrdag``` folder, ```clrdag.zip```, and ```clrdag_*.tar.gz``` in ```Downloads```.

### \*nix users:

0. Make sure you have [R](https://www.r-project.org/) (>= 3.5.0) and building tools installed. 
In R, install the packages ```Rcpp``` and ```RcppArmadillo```.

1. Download [```clrdag.zip```](https://github.com/chunlinli/clrdag/blob/master/pkg/clrdag.zip) to local ```Downloads``` directory. Unzip ```clrdag.zip```.

2. In a terminal app (Bash), 
```
cd "path/to/downloads"
R CMD build clrdag
R CMD INSTALL clrdag_*.tar.gz
```
3. Remove ```clrdag``` directory, ```clrdag.zip```, and ```clrdag_*.tar.gz``` in ```Downloads```.

## Example(s)

Try an example in R.

```
library(clrdag)
set.seed(2018)
p <- 50; n <- 1000; sparsity <- 2/p

## generate a random lower triangular adjacnecy matrix
A <- matrix(rbinom(p*p, 1, sparsity)*sign(runif(p*p, -1, 1))*runif(p*p, 0.7, 1), p, p)
A[upper.tri(A, diag=TRUE)] <- 0

## permute the order of adjacency matrix
idx <- sample(1:p)
A <- A[idx,idx]

## num of edges in A
sum(A != 0)

## data matrix
X <- matrix(rnorm(n*p), n, p) %*% t(solve(diag(p) - A))

## estimate the graph
t <- proc.time()
out <- MLEdag(X=X, tau=0.3, mu=1, rho=1.2)
proc.time() - t

## Frobenius distance to the truth adjacency matrix
sum((out$A - A)^2)

## Hamming distance to the truth graph
sum(abs((out$A != 0) - (A != 0)))

## test edge 1 --> 2
D <- matrix(0, p, p)
D[2,1] <- 1
out <- MLEdag(X=X, D=D, tau=0.3, mu=1, rho=1.2)
out$pval

## test edge 1 --> 3
D <- matrix(0, p, p)
D[3,1] <- 1
out <- MLEdag(X=X, D=D, tau=0.3, mu=1, rho=1.2)
out$pval
```

The [```examples```](https://github.com/chunlinli/clrdag/tree/master/examples) directory contains Examples 1-2 in the paper. 
You can reproduce the other examples in the paper by changing the variables ```n```, ```p```, ```sparsity```, and ```D```. 

## To-dos in future versions

Cross-validation functions, information criterion functions (AIC, BIC), diagnostic plots. 

## More about ```clrdag``` 

The program was originally written in July 2018. 
Ziyue Zhu (thanks to her!) and I refactored the R code version in December 2018 as part of the final project for the optimization course EE 5239. 
As of May 1, 2019, this R package is publicly available. 
It is still in development and be cautious.

To report an issue, please file it at [here](https://github.com/chunlinli/clrdag/issues).
