## load the previous SaRa program code
#dyn.load("diagnosticValue64.dll")
#source("llce.r")

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
    p <- 2 * pnorm(abs(x), mean = 0, sd = sd.est, lower.tail = FALSE)
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



##estimate.sd <- function(x, h = 20) {
##  sds <- apply(x, 1, estimateSigma, h = h)
##  return(sds)
##}

estimate.sd <- function(x) {
  sds <- apply(x, 1, sd)
  return(sds)
}


## Multiple U kernel starts here
multi.scanU <- function(x, l, sd.est) {
  T <- ncol(x)
  N <- nrow(x)
  res <- matrix(NA, N, T - l + 1)
  for(i in 1:N) {
    sums <- cumsum(c(0, x[i, ]))
    res[i, ] <- sums[(1 + l):(T + 1)] - sums[1:(T - l + 1)]
  }
  means <- rowMeans(x)
  res <- (res - means * l)/sd.est/sqrt(l * (1 - l/T))
  return(res)
}

multi.cumsum <- function(x){
  T <- ncol(x)
  N <- nrow(x)
  sums1 <- matrix(0, nrow = N, ncol = T)
  for(i in 1:N) {
    sums1[i, ] <- cumsum(x[i, ])
  }
  sums <- cbind(0, sums1)
  return(sums)
}

multi.scanU.2 <- function(sums, l, sd.est, means = NULL) {
  T <- ncol(sums) - 1
  N <- nrow(sums)
  local.sums <- sums[, (l + 1):(T + 1)] - sums[, 1:(T - l + 1)]
  if (is.null(means)) means <- sums[, T + 1]/T
  res <- (local.sums - means * l)/sd.est/sqrt(l * (1 - l/T))
  return(res)
}


## Higher Criticism ##
higher.criticism <- function(p,alpha0)
{
  N <- nrow(p)
  T <- ncol(p)
  HC <- rep(0, T)
  for (i in 1:T) {
    pp <- sort(p[, i])
    W <- ((1:N)/N - pp)/sqrt(pp * (1 - pp)) * sqrt(N)
    HC[i] <- max(W[alpha0:floor(N/2)])
  }
  return(HC)
}


## Weighted sum of squared statistics ##
weighted.sum <- function (U, p0 = 1){
  U.dim1 <- dim(U)[1]
  U.dim2 <- dim(U)[2]
  U.squared <- U^2
  weights <- exp(U.squared/2)/((1 - p0)/p0 + exp(U.squared/2))
  res <- colSums(weights * U.squared/2)
  return(res)
}

## My function to estimate the trend in data ##
mytrend <- function(x, h = 100) {
  T <- length(x)
  sums <- cumsum(c(rep(0, h + 1), x, rep(0, h)))
  sum.loc <- sums[(1:T) + 2 * h + 1] - sums[1:T]
  ns <- c(h + (1:h), rep(2 * h + 1, T - (2 * h)), h + (h:1))
  xbar <- sum.loc/ns
  return(xbar)
}

detrend <- function(x, ...) {
    N <- nrow(x)
    T <- ncol(x)
    x.detrend <- matrix(0, N, T)
    for(i in 1:N) {
        x.detrend[i, ] <- x[i, ] - mytrend(x[i, ], ...)
    }
    return(x.detrend)
}


