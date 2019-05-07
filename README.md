## Introduction

<img src="https://www.r-pkg.org/badges/version/clrdag">
<img src="https://cranlogs.r-pkg.org/badges/grand-total/clrdag">

`clrdag` is an R package for likelihood estimation and inference of a Gaussian directed acyclic graph. 
See [`clrdag-slides.pdf`](https://github.com/chunlinli/clrdag/blob/master/clrdag-slides.pdf) for the introduction slides to the package. 
The complete reference for the methodology is 

Li, C., Shen, X., and Pan, W. (2019). Likelihood ratio tests for a large directed acyclic graph. *Journal of the American Statistical Association*. Accepted. 

~~See also the package vignette by `vignette('clrdag')` in R.~~

## Installation 

~~It is strongly recommended to install `clrdag` from [CRAN](https://cran.r-project.org/).~~ Alternatively, you may follow the instructions below. 

To achieve best performance, you may want to install the latest version of [OpenBLAS](https://github.com/xianyi/OpenBLAS) or [Intel Math Kernel Library](https://software.intel.com/mkl). You are also recommended to use a compiler that fully supports [OpenMP](https://www.openmp.org/). 

### Windows users: 

0. Make sure you have [R](https://www.r-project.org/) (>= 3.5.0) and [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (>= 35) installed. In particular, add R and Rtools to your PATH variable. In R, install the packages `Rcpp` and `RcppArmadillo`.

1. Download [`clrdag.zip`](https://github.com/chunlinli/clrdag/blob/master/pkg/clrdag.zip) to local `Downloads` folder. Unzip `clrdag.zip` in `Downloads`.

2. In Command Prompt, 
```
cd "path\\to\\Downloads"
R CMD build clrdag
R CMD INSTALL clrdag_*.tar.gz
```
3. Delete `clrdag` folder, `clrdag.zip`, and `clrdag_*.tar.gz` in `Downloads`.

### \*nix users:

0. Make sure you have [R](https://www.r-project.org/) (>= 3.5.0) and building tools installed. 
In R, install the packages `Rcpp` and `RcppArmadillo`.

1. Download [`clrdag.zip`](https://github.com/chunlinli/clrdag/blob/master/pkg/clrdag.zip) to local `Downloads` directory. Unzip `clrdag.zip` in `Downloads`.

2. In a terminal app (Bash), 
```
cd "path/to/Downloads"
R CMD build clrdag
R CMD INSTALL clrdag_*.tar.gz
```
3. Remove `clrdag` directory, `clrdag.zip`, and `clrdag_*.tar.gz` in `Downloads`.

## Example(s)

In this example, a 50 by 50 random adjacency matrix is generated, 
where there are 43 edges in the graph. 
For each edge, the strength is uniformly distributed in [-1,-0.7] U [0.7,1]. The random error is standard normal. 

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
sum(A != 0) # 43

## data matrix
X <- matrix(rnorm(n*p), n, p) %*% t(solve(diag(p) - A))

## estimate the graph
t <- proc.time()
out <- MLEdag(X=X, tau=0.3, mu=1, rho=1.2)
proc.time() - t 

## Frobenius distance to the truth adjacency matrix
sum((out$A - A)^2) # 0.0285789

## Hamming distance to the truth graph
sum(abs((out$A != 0) - (A != 0))) # 0

## test edge 1 --> 3
D <- matrix(0, p, p)
D[3,1] <- 1
out <- MLEdag(X=X, D=D, tau=0.3, mu=1, rho=1.2)
out$pval # 0.7496623

## test edge 7 --> 4
D <- matrix(0, p, p)
D[4,7] <- 1
out <- MLEdag(X=X, D=D, tau=0.3, mu=1, rho=1.2)
out$pval # 8.827349e-155
```

The [`examples`](https://github.com/chunlinli/clrdag/tree/master/examples) directory contains Examples 1-2 in Li et al. (2019). 
You can reproduce the other examples in the paper by changing the variables `n`, `p`, `sparsity`, `A`, and `D`, accordingly. 

## To-dos

Cross-validation functions, information criterion functions (AIC, BIC), and diagnostic plots. 

## More about `clrdag` 

The program was originally written in July 2018. 
Ziyue Zhu and I refactored the R code and wrote a vignette in December 2018 as part of the final project for the optimization course EE 5239. 
As of May 1, 2019, this R package is publicly available. 
It is still in development and be cautious.

Special thanks to Prof Charles Geyer and Prof Mingyi Hong for their comments.

To report an issue, please file one at [here](https://github.com/chunlinli/clrdag/issues).
