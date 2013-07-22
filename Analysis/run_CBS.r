library(DNAcopy)
library(doSNOW)

for(i in 22:1) {
    cat(i, "\n")
    load(paste("../chr_", i, ".rdata", sep = ""))
    sample.names <- colnames(data.sort)[-(1:3)]
    rm(data.sort)
    gc()

    load(paste("../detrend/chr_", i, "_detrend.rdata", sep = ""))


    ptm <- proc.time()
    cl <- makeCluster(11)
    registerDoSNOW(cl)

    cbs.result <- foreach(j = 1:nrow(data.detrend), .packages = c("DNAcopy")) %dopar% {
        dd <- matrix(data.detrend[j, ], ncol = 1)
        cna.data <- CNA(genomdat = dd, chrom = probe$chromosome, maploc = probe$position, data.type = "logratio", sampleid = sample.names[j], presorted = TRUE)
        segment(cna.data, alpha = 0.01)
    }

    stopCluster(cl)
    ttt <- proc.time() - ptm
    print(ttt)

    save(cbs.result, file = paste("cbs_results/chr_", i, "_cbs.rdata", sep = ""))
}

##nna <- !apply(is.na(data.sort), 1, any)
##cna.data.orig <- CNA(genomdat = data.sort[nna, 4:8], chrom = probe$chromosome, maploc = probe$position, data.type = "logratio", sampleid = sample.names[1:5], presorted = TRUE)

##cbs.result.orig <- segment(cna.data.orig, alpha = 0.01)
