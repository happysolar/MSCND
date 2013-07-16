## load the previous SaRa program code
dyn.load("diagnosticValue.so")
source("llce.r")

## Multi-SaRa starts here
multi.llce <- function(x, h, sd.est)
{
    if(is.vector(x)) x <- matrix(x, nrow = 1)
    N <- nrow(x)
    T <- ncol(x)
    diag.stat <- matrix(NA, nrow = N, ncol = T)
    for(i in 1:N) {
        diag.stat[i, ] <- llce(x[i, ], h)
    }
    return(diag.stat/sd.est * sqrt(h/2))
}

## Test multi.llce
##xx <- matrix(rnorm(10000 * 2000), 2000, 10000)
##system.time(stat.xx <- multi.llce(xx, 10))
##   user  system elapsed 
##  4.062   0.398   4.459 

zmat.pval <- function(x, sd.est = 1)
{
    if(is.vector(x)) x <- matrix(x, nrow = 1)
    N <- nrow(x)
    T <- ncol(x)
    p <- pnorm(abs(x), mean = 0, sd = sd.est, lower.tail = FALSE)
}

## Test llce.pval
##pval.xx <- llce.pval(stat.xx, 1)

##estimateSigma <-
##function(Y,h=10){                   #constant case
##  n     = length(Y)                                 # can we make it faster?
##  YBar  = rep(0,n)
##  for (i in 1:n) {
##     a       = min(n,i+h)
##     b       = max(1,i-h)
##     YBar[i] = mean(Y[b:a])
##  }
##  return(sd(Y-YBar))
##}

estimateSigma <- function(x, h = 10) {
  T <- length(x)
  sums <- cumsum(c(rep(0, h + 1), x, rep(0, h)))
  sum.loc <- sums[(1:T) + 2 * h + 1] - sums[1:T]
  ns <- c(h + (1:h), rep(2 * h + 1, T - (2 * h)), h + (h:1))
  xbar <- sum.loc/ns
  res <- x - xbar
  var0 <- var(res) * (2 * h + 1)/(2 * h)
  return(sqrt(var0))
}

estimate.sd <- function(x, h = 10) {
  sds <- apply(x, 1, estimateSigma)
  return(sds)
}

