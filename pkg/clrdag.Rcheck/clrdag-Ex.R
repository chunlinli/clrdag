pkgname <- "clrdag"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "clrdag-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('clrdag')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MLEdag")
### * MLEdag

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MLEdag
### Title: MLE/LRT of a Gaussian directed acyclic graph
### Aliases: clrdag MLEdag
### Keywords: Directed acyclic graph Likelihood ratio test Maximum
###   likelihood estimation

### ** Examples

##
## Example: random graph
##
library(clrdag)
set.seed(2019)
p<-10
n<-1000
## random graph: randomly generate adjacency matrix A, A lower triangular
sparsity <- 2/p
A <- matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1)),p,p)
A[upper.tri(A,diag=TRUE)] <- 0
X <- matrix(rnorm(n*p),n,p) 
out <- MLEdag(X=X,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE) # compute the MLE
sum(abs((out$A!=0)-(A!=0))) # Hamming distance to the truth graph

# test edge 1 --> 2
D <- matrix(0,p,p)
D[2,1] <- 1
out <- MLEdag(X=X,D=D,tau=0.3,mu=1,rho=1.2,trace_obj=FALSE) # compute the MLE
out$pval



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MLEdag", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
