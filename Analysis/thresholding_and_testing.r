source("../MSCND/mSARA/multi_Scan.r")
source("../MSCND/mSARA/getcutoff_func.r")



for (i in 22:1) {
  print(i)
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
  #### HC picked are the same as those for WSum and Sum.
  WSum.over <- sum(unlist(U.WSum)>th.WSum, na.rm = TRUE)
  Sum.over <- sum(unlist(U.Sum)>th.Sum, na.rm = TRUE)
  HC.over <- max(WSum.over, Sum.over)
  th.HC <- mean(sort(unlist(U.HC), decreasing = TRUE)[HC.over+c(0,1)])

  #### TODO: Plug in threshold here
  ch.regions.HC <- thres.U(U.HC, th.HC)
  U.test.HC <- calc.U.region(data.detrend, ch.regions.HC)
  rownames(U.test.HC) <- sample.names
  save(th.Sum, ch.regions.Sum, U.test.Sum, th.WSum, ch.regions.WSum, U.test.WSum, th.HC, ch.regions.HC, U.test.HC, file=paste("scan_results/chr_", i, "_region_test.rdata",sep=""))
}

  #### TODO (tomorrow): test U using CNVtools


source("../MSCND/mSARA/multi_Scan.r")
load("thres_D_result.rdata")
load(paste("../chr_", 22, ".rdata", sep = ""))
sample.names <- colnames(data.sort)[-(1:3)]
rm(data.sort)
gc()


for(chr in 1:22) {
    print(chr)
    load(paste("scan_results/chr_", chr, "_scan_results.rdata", sep = ""))
    load(paste("../detrend/chr_", chr, "_detrend.rdata", sep = ""))
    ##load(paste("chr_", i, "_scan_results.rdata", sep = ""))
    ##load(paste("chr_", i, "_detrend.rdata", sep = ""))
    ##load("thres_D_result.rdata")


    ##N <- nrow(data.detrend)
    ##T <- ncol(data.detrend)

    chp.max.Sum <- localmax.D(D.Sum, c(20, 40, 60), quantile.D.Sum[, 990])
    D.test.Sum <- calc.D.chp(data.detrend, chp.max.Sum, c(10, 20, 30))
    for (i in 1:3){
        rownames(D.test.Sum[[i]]) <- sample.names
    }

    chp.max.WSum <- localmax.D(D.WSum, c(20, 40, 60), quantile.D.WSum[, 990])
    D.test.WSum <- calc.D.chp(data.detrend, chp.max.WSum, c(10, 20, 30))
    for (i in 1:3){
        rownames(D.test.WSum[[i]]) <- sample.names
    }

    chp.max.HC <- localmax.D(D.HC, c(20, 40, 60), quantile.D.HC[, 990])
    D.test.HC <- calc.D.chp(data.detrend, chp.max.HC, c(10, 20, 30))
    for (i in 1:3){
        rownames(D.test.HC[[i]]) <- sample.names
    }
    save(chp.max.Sum, D.test.Sum, chp.max.WSum, D.test.WSum, chp.max.HC, D.test.HC, file=paste("test_D_results/chr_", chr, "_chp_test.rdata",sep=""))
    rm(data.detrend)
    gc()
}








  #### TODO: thresholding or FDR chp.max
  #### TODO: merge the change points using different window sizes
  #### TODO calculate D according to the change points


  #### TODO (tomorrow): test D using CNVtools



#######################
## Test U statistics ##
#######################

library(doSNOW)
source("clustering.r")

for (i in 22:1) {
    cat(i, "\n")

    load(paste("scan_results/chr_", i, "_region_test.rdata", sep = ""))

    pheno <- read.table("pheno_GPN_PBR_PLINK.txt", header = TRUE, as.is = TRUE)
    rownames(pheno) <- pheno$IID

    pheno <- pheno[rownames(U.test.Sum), ]

    ##batch.info <- read.delim("batch_info.txt", header = TRUE, as.is = TRUE, row.names = 1)

    ##batch1 <- factor(batch.info[, 1])
    ##names(batch1) <- rownames(batch.info)

    ##batch <- batch1[rownames(U.test.Sum)]
    ##pheno <- cbind(pheno, batch)

    y <- pheno$case_control - 1
    z <- pheno[, -(1:3)]


    cat("Sum\n")
    cl <- makeCluster(11, "SOCK")
    registerDoSNOW(cl)
    n.test.Sum <- ncol(U.test.Sum)
                                        ##ptm <- proc.time()
    test.result.Sum <- foreach(j = 1:n.test.Sum, .packages = c("mclust", "aod")) %dopar% {
        x <- U.test.Sum[, j]
        tryCatch(test.binary(y, x, z, G = 1:5, modelNames = "V"), error = function(ex) NULL)
    }
    ##ttt <- proc.time() - ptm
    ##print(ttt)
    stopCluster(cl)

    cat("WSum\n")
    cl <- makeCluster(11, "SOCK")
    registerDoSNOW(cl)
    n.test.WSum <- ncol(U.test.WSum)
    ##ptm <- proc.time()
    test.result.WSum <- foreach(j = 1:n.test.WSum, .packages = c("mclust", "aod")) %dopar% {
        x <- U.test.WSum[, j]
        tryCatch(test.binary(y, x, z, G = 1:5, modelNames = "V"), error = function(ex) NULL)
    }
                                        ##ttt <- proc.time() - ptm
                                        ##print(ttt)
    stopCluster(cl)

    cat("HC\n")
    cl <- makeCluster(11, "SOCK")
    registerDoSNOW(cl)
    n.test.HC <- ncol(U.test.HC)
    ##ptm <- proc.time()
    test.result.HC <- foreach(j = 1:n.test.HC, .packages = c("mclust", "aod")) %dopar% {
        x <- U.test.HC[, j]
        tryCatch(test.binary(y, x, z, G = 1:5, modelNames = "V"), error = function(ex) NULL)
    }
    ##ttt <- proc.time() - ptm
    ##print(ttt)
    stopCluster(cl)

    save(test.result.Sum, test.result.WSum, test.result.HC, file = paste("test_results/chr_", i, "_test_results_U.rdata", sep = ""))
    gc()
}


#######################
## Test D statistics ##
#######################

library(doSNOW)
source("clustering.r")

for (i in 22:1) {
    cat(i, "\n")

    load(paste("scan_D_results/chr_", i, "_chp_test.rdata", sep = ""))

    pheno <- read.table("pheno_GPN_PBR_PLINK.txt", header = TRUE, as.is = TRUE)
    rownames(pheno) <- pheno$IID

    pheno <- pheno[rownames(D.test.Sum[[1]]), ]

    ##batch.info <- read.delim("batch_info.txt", header = TRUE, as.is = TRUE, row.names = 1)

    ##batch1 <- factor(batch.info[, 1])
    ##names(batch1) <- rownames(batch.info)

    ##batch <- batch1[rownames(U.test.Sum)]
    ##pheno <- cbind(pheno, batch)

    y <- pheno$case_control - 1
    z <- pheno[, -(1:3)]


    cat("Sum\n")
    cl <- makeCluster(11, "SOCK")
    registerDoSNOW(cl)
    test.result.Sum <- list()
    for(k in 1:3) {
        n.test.Sum <- ncol(D.test.Sum[[k]])
        ##ptm <- proc.time()
        test.result.Sum[[k]] <- foreach(j = 1:n.test.Sum, .packages = c("mclust", "aod")) %dopar% {
            x <- D.test.Sum[[k]][, j]
            tryCatch(test.binary(y, x, z, G = 1:5, modelNames = "V"), error = function(ex) NULL)
        }
    }
    ##ttt <- proc.time() - ptm
    ##print(ttt)
    stopCluster(cl)

    cat("WSum\n")
    cl <- makeCluster(11, "SOCK")
    registerDoSNOW(cl)
    test.result.WSum <- list()
    for(k in 1:3) {
        n.test.WSum <- ncol(D.test.WSum[[k]])
        ##ptm <- proc.time()
        test.result.WSum[[k]] <- foreach(j = 1:n.test.WSum, .packages = c("mclust", "aod")) %dopar% {
            x <- D.test.WSum[[k]][, j]
            tryCatch(test.binary(y, x, z, G = 1:5, modelNames = "V"), error = function(ex) NULL)
        }
    }
                                        ##ttt <- proc.time() - ptm
                                        ##print(ttt)
    stopCluster(cl)

    cat("HC\n")
    cl <- makeCluster(11, "SOCK")
    registerDoSNOW(cl)
    test.result.HC <- list()
    for(k in 1:3) {
        n.test.HC <- ncol(D.test.HC[[k]])
        ##ptm <- proc.time()
        test.result.HC[[k]] <- foreach(j = 1:n.test.HC, .packages = c("mclust", "aod")) %dopar% {
            x <- D.test.HC[[k]][, j]
            tryCatch(test.binary(y, x, z, G = 1:5, modelNames = "V"), error = function(ex) NULL)
        }
    }
    ##ttt <- proc.time() - ptm
    ##print(ttt)
    stopCluster(cl)

    save(test.result.Sum, test.result.WSum, test.result.HC, file = paste("test_D_results/chr_", i, "_test_results_D.rdata", sep = ""))
    gc()
}















load(paste("../detrend/chr_", i, "_detrend.rdata", sep = ""))


density.model <- function(m) {
    function(x) {
        d <- 0
        for(i in 1:length(m$parameters$mean)) {
            d <- d + m$parameters$pro[i] * dnorm(x, m$parameters$mean[i], sqrt(m$parameters$variance$sigmasq[i]))
        }
        return(d)
    }
}

plot.model <- function(m, ...) {
    curve(density.model(m)(x), col = 2, ...)
    for(i in 1:length(m$parameters$mean)) {
        f <- function(x) m$parameters$pro[i] * dnorm(x, m$parameters$mean[i], sqrt(m$parameters$variance$sigmasq[i]))
        curve(f, add = TRUE, col = 3)
    }
}
