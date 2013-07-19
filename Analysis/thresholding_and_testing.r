source("../MSCND/mSARA/multi_Scan.r")
source("../MSCND/mSARA/getcutoff_func.r")



i <- 22
load(paste("scan_results/chr_", i, "_scan_results.rdata", sep = ""))
load(paste("../detrend/chr_", i, "_detrend.rdata", sep = ""))
load(paste("../chr_", i, ".rdata", sep = ""))
sample.names <- colnames(data.sort)[-(1:3)]
rm(data.sort)
gc()

N <- nrow(data.detrend)
T <- ncol(data.detrend)


th.Sum <- U.WSum.cutoff(0.01, T, 1, 50, N, 1)
ch.regions.Sum <- thres.U(U.Sum, th.Sum)
U.test.Sum <- calc.U.region(data.detrend, ch.regions.Sum)
rownames(U.test.Sum) <- sample.names

th.WSum <- U.WSum.cutoff(0.01, T, 1, 50, N, 0.01)
ch.regions.WSum <- thres.U(U.WSum, th.WSum)
U.test.WSum <- calc.U.region(data.detrend, ch.regions.WSum)
rownames(U.test.WSum) <- sample.names

#### Find the threshold such that the number of 
#### HC picked are the same as that for WSum and Sum. 
WSum.over <- sum(U.WSum > th.WSum)
Sum.over <- sum(U.Sum > th.Sum)
HC.over <- max(WSum.over, Sum.over)
th.HC <- mean(sort(U.HC, decreasing = TRUE)[HC.over+c(0,1)])

#### TODO: Plug in threshold here
ch.regions.HC <- thres.U(U.HC, th.HC)
U.test.HC <- calc.U.region(data.detrend, ch.regions.HC)
rownames(U.test.HC) <- sample.names


#### TODO (tomorrow): test U using CNVtools



chp.max <- localmax.D(D.Sum, c(20, 40, 60))
#### TODO: thresholding or FDR chp.max
#### TODO: merge the change points using different window sizes
#### TODO calculate D according to the change points


#### TODO (tomorrow): test D using CNVtools
