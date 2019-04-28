## Installation 

This package requires compilation to install from source. The following instructions are for new R users. 

### Windows users: 

0. Make sure you have R (>= 3.5.3) and Rtools35 or later (available [here](https://cran.r-project.org/bin/windows/Rtools/)) properly installed. In particular, add R and Rtools to your PATH variable. In R, install the packages ```Rcpp``` and ```RcppArmadillo```.

1. Download ```clrdag_*.tar.gz``` from ```pkg``` directory.

2. In Command Prompt, 
```
cd "path\\to\\downloaded\\package\\" && R CMD INSTALL clrdag_*.tar.gz
```
3. Delete ```clrdag_*.tar.gz```.

### Unix(-like) users:

0. Make sure you have R (>= 3.5.3) and building tools installed. 
In R, install the packages ```Rcpp``` and ```RcppArmadillo```.

1. Download ```clrdag_*.tar.gz``` from ```pkg``` directory.

2. In terminal (Bash), 
```
cd "path/to/downloaded/package/" && R CMD INSTALL clrdag_*.tar.gz
```
3. Remove ```clrdag_*.tar.gz```.

## To-dos in future versions

Cross-validation functions, Diagnostic plots. 

## About ```clrdag``` 

This is a C++/R implementation of the constrained likelihood ratio tests of a DAG.
The reference for the details is 

Li, C., Shen, X., and Pan, W. (2019). Likelihood ratio tests of a large directed acyclic graph. Submitted. 

The program was originally written by me in Summer 2018. Then Ziyue Zhu (thanks to her!) and I refactored the R code version (available [here](https://github.umn.edu/li000007/clrdag_r/)) in December 2018 as part of the final project for the optimization course EE5239. As of April 2019, this R package can be accessed by UMN students. It is not yet published and is still in development.

To report an issue, please file an issue at [here](https://github.umn.edu/li000007/clrdag/issues).
