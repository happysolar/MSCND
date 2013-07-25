#!/usr/local/cluster/software/installation/R/R-2.15.1/bin/Rscript
library(DNAcopy)
library(mscp)
source("llce.r")
source("multi_Scan.r")
dyn.load("diagnosticValue.so")

args <- commandArgs(TRUE)
seed <- as.numeric(args[[1]])
out <- args[[2]]

##seed <- 1
##out <- "simu_1"
set.seed(seed)

N <- 1000
T <- 500

cp <- c(27, 54, 115, 130, 221, 306)

delta <- c(1.29, -0.96, 0.87, -0.75)

p <- c(0.02, 0.05, 0.1)

H <- c(10, 20, 30)

alpha <- 0.005
cutoff.Sum <- cutoff.WSum <- cutoff.HC <- rep(NA, 3)
for(j in 1:length(H)) {
    h <- H[j]
    load(paste("../../simu_thres/summary/thres_D_result_N", N, "_h", h, "_hh", h, ".rdata", sep = ""))
    cutoff.Sum[j] <- quantile.D.Sum[ceiling(1000 * (1 - alpha))]
    cutoff.WSum[j] <- quantile.D.WSum[ceiling(1000 * (1 - alpha))]
    cutoff.HC[j] <- quantile.D.HC[ceiling(1000 * (1 - alpha))]
}




sigma <- 0.5
tau <- 2
##sigma1 <- 0
##sigma2 <- 0
a1 <- 0.025
a2 <- 0.01


mu <- matrix(0, N, T)
for(i in 1:N) {
    phase <- runif(1, 0, 2 * pi)
    mu[i, ] <- 0.1 * sigma * sin(a1 * pi * (1:T)) + 0.1 * tau * sigma * sin(a2 * pi * (1:T) + phase)
}


n <- N * p
carr <- list()
for(i in 1:length(p)) {
    carrier <- sort(sample(1:N, n[i]))
    mu[carrier, (cp[2 * i - 1] + 1):cp[2 * i]] <- mu[carrier, (cp[2 * i - 1] + 1):cp[2 * i]] + delta[i]
    carr[[i]] <- carrier
}




x <- mu + rnorm(N * T, sd = sigma)

DDs <- list()

sigma.est <- estimate.sd(x)

Sum <- WSum <- HC <- matrix(NA, length(H), T)
for(i in 1:length(H)) {
    h <- H[i]
    Ds <- multi.llce(x, h, sigma.est)
    DDs[[i]] <- Ds
    Sum[i, ] <- weighted.sum(Ds)
    WSum[i, ] <- weighted.sum(Ds, p0 = 0.01)
    Ds.p <- zmat.pval(Ds)
    HC[i, ] <- higher.criticism(Ds.p, alpha0 = 4)
}

cp.Sum <- localmax.D(Sum, H, cutoff.Sum)
cp.WSum <- localmax.D(WSum, H, cutoff.WSum)
cp.HC <- localmax.D(HC, H, cutoff.HC)

##BIC.Sum <- get.all.BICs(x, cp.Sum)
##BIC.WSum <- get.all.BICs(x, cp.WSum)
##BIC.HC <- get.all.BICs(x, cp.HC)
##BIC.Sum.best <- BIC.WSum.best <- BIC.HC.best <- list() 
##for(i in 1:N) {
##    BIC.Sum.best[[i]] <- iterative.remove.cp(BIC.Sum[[i]])
##    BIC.WSum.best[[i]] <- iterative.remove.cp(BIC.WSum[[i]])
##    BIC.HC.best[[i]] <- iterative.remove.cp(BIC.HC[[i]])
##}


##s <- 2
##cp.Sum.unique <- NULL
##for (i in 1:length(H)) {
##    cp.Sum.unique <- sort(c(cp.Sum.unique, setdiff(cp.Sum[[i]], rep(cp.Sum.unique, each = 2 * s + 1) + (-s:s))))
##}

##cp.WSum.unique <- NULL
##for (i in 1:length(H)) {
##    cp.WSum.unique <- sort(c(cp.WSum.unique, setdiff(cp.WSum[[i]], rep(cp.WSum.unique, each = 2 * s + 1) + (-s:s))))
##}

##cp.HC.unique <- NULL
##for (i in 1:length(H)) {
##    cp.HC.unique <- sort(c(cp.HC.unique, setdiff(cp.HC[[i]], rep(cp.HC.unique, each = 2 * s + 1) + (-s:s))))
##}



##cp.Sum.keep <- soft.thresholding.mu(x, cp.Sum.unique, 0.7)
cp.Sum.keep <- thresholding.by.h(x, cp.Sum, 2 * sigma.est * sqrt(2/min(H)), 2)
cp.WSum.keep <- thresholding.by.h(x, cp.WSum, 2 * sigma.est * sqrt(2/min(H)), 2)
cp.HC.keep <- thresholding.by.h(x, cp.HC, 2 * sigma.est * sqrt(2/min(H)), 2)

cp.Sum.all <- sort(unique(unlist(cp.Sum.keep)))
cp.WSum.all <- sort(unique(unlist(cp.WSum.keep)))
cp.HC.all <- sort(unique(unlist(cp.HC.keep)))


### CBS
cna <- CNA(t(x), maploc = 1:T, chrom = rep(1, ncol(x)))
cbs <- segment(cna)


### mSaRa
cut <- 2 * sigma * sqrt(2/H)
cp.msara <- list()
for(i in 1:N) {
    DD <- matrix(NA, length(H), T)
    for(j in 1:length(H)) {
        DD[j, ] <- DDs[[j]][i, ]
    }
    temp <- sort(unique(unlist(localmax.D(DD, H, cut))))
    BIC <- get.one.BIC(x[i, ], temp)
    cp.msara[[i]] <- iterative.remove.cp(BIC)$cp
}


### MSCBS
mscbs.result <- mscbs(t(x), MIN.SNPs = 10)

save(x, carr, cp.Sum.all, cp.WSum.all, cp.HC.all, cp.Sum, cp.WSum, cp.HC, cp.Sum.keep, cp.WSum.keep, cp.HC.keep, cbs, cp.msara, mscbs.result, file = paste("result/", out, "/simu_res_", seed, ".rdata", sep = ""))
