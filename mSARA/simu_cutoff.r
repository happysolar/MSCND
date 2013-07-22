source("multi_Scan.r")
dyn.load("diagnosticValue.so")
source("llce.r")
### Simulate the cutoff for D.HC, D.Sum, D.WSum
N <- 1894
T <- 60*100


seed <- 1
#ptm <- proc.time()
set.seed(seed)
data.simu <- matrix(rnorm(N*T), nrow = N)
hh <- c(10, 20, 30)
D.HC <- matrix(0, nrow=3, ncol=T)
D.Sum <- matrix(0, nrow=3, ncol=T)
D.WSum <- matrix(0, nrow=3, ncol=T)
for (j in 1:3){
    Ds <- multi.llce(data.simu, hh[j], sd.est=1)
    D.Sum[j, ] <- weighted.sum(Ds)
    D.WSum[j, ] <- weighted.sum(Ds, p0=0.01)
    p.Ds <- zmat.pval(Ds)
    rm(Ds)
    gc()
    D.HC[j, ] <- higher.criticism(p.Ds, alpha0=4)
    rm(p.Ds)
    gc()
}
rm(data.simu)
gc()

chp.max.Sum <- localmax.D(D.Sum, c(20, 40, 60))
D.Sum.Max <- list()
quantile.D.Sum <- matrix(NA, 3, 1000)
for (i in 1:3){
  D.Sum.Max[[i]] <- D.Sum[i, chp.max.Sum[[i]]]
#  quantile.D.Sum[i, ] <- quantile(D.Sum.Max[[i]], c(1:1000)/1000)
}

chp.max.WSum <- localmax.D(D.WSum, c(20, 40, 60))
D.WSum.Max=list()
quantile.D.WSum <- matrix(NA, 3, 1000)
for (i in 1:3){
  D.WSum.Max[[i]] <- D.WSum[i, chp.max.WSum[[i]]]
#  quantile.D.WSum[i, ] <- quantile(D.WSum.Max[[i]], c(1:1000)/1000)
}

chp.max.HC <- localmax.D(D.HC, c(20, 40, 60))
D.HC.Max=list()
quantile.D.HC <- matrix(NA, 3, 1000)
for (i in 1:3){
  D.HC.Max[[i]] <- D.HC[i, chp.max.HC[[i]]]
#  quantile.D.HC[i, ] <- quantile(D.HC.Max[[i]], c(1:1000)/1000)
}
#ttt <- proc.time() - ptm
#ttt



save(D.Sum.Max, D.WSum.Max, D.HC.Max, file = paste("../MSCND/simu_threshold_", seed,".rdata", sep = ""))

#save(D.Sum.Max, D.WSum.Max, D.HC.Max, quantile.D.Sum, quantile.D.WSum, quantile.D.HC, file='../MSCND/simu_threshold.rdata')
