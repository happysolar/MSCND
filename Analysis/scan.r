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

mytrend <- function(x, h = floor(length(x)/20)) {
  T <- length(x)
  sums <- cumsum(c(rep(0, h + 1), x, rep(0, h)))
  sum.loc <- sums[(1:T) + 2 * h + 1] - sums[1:T]
  ns <- c(h + (1:h), rep(2 * h + 1, T - (2 * h)), h + (h:1))
  xbar <- sum.loc/ns
  return(xbar)
}


for(i in 1:22) {
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
    N <- nrow(data.int)
    T <- ncol(data.int)
    data.detrend <- matrix(0, N, T)
    for(j in 1:N) {
        data.detrend[j, ] <- data.int[j, ] - mytrend(data.int[j, ])
    }
    save(data.detrend, probe, file = paste("../detrend/chr_", i, "_detrend.rdata", sep = ""))
}
