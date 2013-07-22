## Load data
load("../chr_1.rdata")

## Load the C function for the SaRa kernel
dyn.load("diagnosticValue64.dll")

## Load the SaRa function
source("../MSCND/mSARA/llce.r")

## Load the multi-Scan functions
source("../MSCND/mSARA/multi_Scan.r")

probe.info <- data.sort[, 1:3]
data.intensity <- t(data.sort[, -(1:3)])
nna <- apply(!is.na(data.intensity), 2, all)
probe <- probe.info[nna, ]
data.int <- data.intensity[, nna]
rm(data.sort)
rm(probe.info)
rm(data.intensity)
gc()

system.time(sd.est <- estimate.sd(data.int))
h <- c(10, 20, 30)
system.time(Ds1 <- multi.llce(data.int, 10, sd.est))
##system.time(Us1 <- multi.scanU(data.int, 10, sd.est))

##system.time(sums <- t(rbind(0, apply(data.int, 1, cumsum))))
T <- ncol(data.int)
N <- nrow(data.int)
sums1 <- matrix(nrow = N, ncol = T)
for(i in 1:N) {
    sums1[i, ] <- cumsum(data.int[i, ])
}
sums <- cbind(0, sums1)

##sd.est2 <- apply(data.int, 1, sd)
l <- 10
system.time(Us1 <- sums[, (l + 1):(T + 1)] - sums[, 1:(T - l + 1)])
means <- rowMeans(data.int)
Us <- (Us1 - means * l)/sd.est/sqrt(l * (1 - l/T))

data.trend <- matrix(0, N, T)
for(i in 1:10) {
    print(i)
    data.trend[i, ] <- ksmooth(1:T, data.int[i, ], kernel = "box", bandwidth = T/10)$y
}

data.temp <- data.int[1:10, ] - data.trend[1:10, ]
sd22 <- estimate.sd(data.temp)
U22 <- multi.scanU(data.temp, 10, sd22)

par(mfrow = c(3, 1))

#################
#### Detrend ####
#################



## batch 1: chr1~4
## batch 2: chr5~11
## batch 3: chr12~22

for(i in 22:12) {
    print(i)
    load(paste("../chr_", i, ".rdata", sep = ""))
    probe.info <- data.sort[, 1:3]
    data.intensity <- t(data.sort[, -(1:3)])
    nna <- apply(!is.na(data.intensity), 2, all)
    probe <- probe.info[nna, ]
    data.int <- data.intensity[, nna]
    rm(data.sort)
    rm(probe.info)
    rm(data.intensity)
    gc()
    ##N <- nrow(data.int)
    ##T <- ncol(data.int)
    ##data.detrend2 <- matrix(0, N, T)
    ##for(j in 1:N) {
    ##    data.detrend2[j, ] <- data.int[j, ] - mytrend(data.int[j, ])
    ##}
    data.detrend <- detrend(data.int)
    save(data.detrend, probe, file = paste("../detrend/chr_", i, "_detrend.rdata", sep = ""))
    gc()
##}

##################
#### Analysis ####
##################

##for (i in 22:10) {
    ##load(file = paste("../detrend/chr_", i, "_detrend.rdata", sep = ""))
    N <- nrow(data.detrend)
    T <- ncol(data.detrend)
    sd.est <- estimate.sd(data.detrend)
    hh <- c(10, 20, 30)
    D.HC <- matrix(0, nrow=3, ncol=T)
    D.Sum <- matrix(0, nrow=3, ncol=T)
    D.WSum <- matrix(0, nrow=3, ncol=T)
    for (j in 1:3){
        Ds <- multi.llce(data.detrend, hh[j], sd.est)
        p.Ds <- zmat.pval(Ds)
        D.HC[j, ] <- higher.criticism(p.Ds, alpha0=4)
        D.Sum[j, ] <- weighted.sum(Ds)
        D.WSum[j, ] <- weighted.sum(Ds, p0=0.01)
    }
    U.HC <- list()
    U.Sum <- list()
    U.WSum <- list()
    L <- 50
    ss <- multi.cumsum(data.detrend)
    for (ll in 1:L) {
        #print(ll)
        Us <- multi.scanU.2(ss, ll, sd.est)
        p.Us <- zmat.pval(Us)
        U.HC[[ll]] <- higher.criticism(p.Us, alpha0=4)
        U.Sum[[ll]] <- weighted.sum(Us)
        U.WSum[[ll]] <- weighted.sum(Us, p0=0.01)
    }
    save(D.HC, D.Sum, D.WSum, U.HC, U.Sum, U.WSum, file = paste("scan_results/chr_", i, "_scan_results.rdata", sep = ""))
    gc()
}
