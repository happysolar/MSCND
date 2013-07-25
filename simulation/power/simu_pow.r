#!/usr/local/cluster/software/installation/R/R-2.15.1/bin/Rscript
source("multi_Scan.r")
source("llce.r")
dyn.load("diagnosticValue.so")

args <- commandArgs(TRUE)

ind <- as.numeric(args[[1]])
set.seed(ind)

T <- 100

NN <- c(100, 500, 1000)
Delta <- c(0.25, 0.5, 0.75, 1, 1.5)
P <- (0:25)/100

param <- expand.grid(N = NN, delta = Delta, p = P)
N <- param[ind, 1]
delta <- param[ind, 2]
p <- param[ind, 3]

print(paste(ind, N, delta, p))

H <- c(10, 20)

##load(paste("../../simu_thres/summary/thres_D_result_N", N, "_h", h, ".rdata", sep = ""))

B <- 1000


res <- list()

Sum.max <- WSum.max <- HC.max <- matrix(NA, B, length(H))
for(i in 1:B) {
    print(i)
    n <- N * p
    carrier <- sample(1:N, n)
    x <- matrix(rnorm(N * T), N, T)
    x[carrier, 1:(T/2)] <- x[carrier, 1:(T/2)] + delta
    Sum <- WSum <- HC <- matrix(NA, length(H), T)
    for(j in 1:length(H)) {
        h <- H[j]
        Ds <- multi.llce(x, h, 1)
        Ds.pval <- zmat.pval(Ds)
        Sum[j, ] <- weighted.sum(Ds)
        WSum[j, ] <- weighted.sum(Ds, p0 = 0.01)
        HC[j, ] <- higher.criticism(Ds.pval, alpha0 = 4)
    }
    Sum.lm <- localmax.D(Sum, H, rep(0, length(H)))
    WSum.lm <- localmax.D(WSum, H, rep(0, length(H)))
    HC.lm <- localmax.D(HC, H, rep(0, length(H)))
    for(j in 1:length(H)) {
        Sum.max[i, j] <- max(Sum[j, Sum.lm[[j]]])
        WSum.max[i, j] <- max(WSum[j, WSum.lm[[j]]])
        HC.max[i, j] <- max(HC[j, HC.lm[[j]]])
    }
}

save(Sum.max, WSum.max, HC.max, file = paste("output/pow_N", N, "_delta", delta, "_p", p, ".rdata", sep = ""))
