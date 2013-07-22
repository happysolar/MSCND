source("mis_prob.r")


#########################
## Summarize U results ##
#########################

summary.U.Sum <- list()
summary.U.WSum <- list()
summary.U.HC <- list()


for(i in 1:22) {
    print(i)
    load(paste("test_results/chr_", i, "_test_results_U.rdata", sep = ""))
    pval.Sum <- sapply(test.result.Sum, function(x) ifelse(is.null(x), NA, x$p.value))
    pval.WSum <- sapply(test.result.WSum, function(x) ifelse(is.null(x), NA, x$p.value))
    pval.HC <- sapply(test.result.HC, function(x) ifelse(is.null(x), NA, x$p.value))

    par.Sum <- lapply(test.result.Sum, function(x) x$clst.para)
    par.WSum <- lapply(test.result.WSum, function(x) x$clst.para)
    par.HC <- lapply(test.result.HC, function(x) x$clst.para)

    modes.Sum <- list()
    for(j in 1:length(par.Sum)) {
        x <- par.Sum[[j]]
        if(is.null(x)) {
            modes.Sum[[j]] <- NA
        } else {
            modes.Sum[[j]] <- modes(x$mean, x$variance$sigmasq, x$pro)
        }
    }
    n.modes.Sum <- sapply(modes.Sum, function(x) ifelse(all(is.na(x)), NA, x$nmode))



    mis.Sum <- rep(NA, length(par.Sum))
    for(j in 1:length(par.Sum)) {
        x <- par.Sum[[j]]
        mis.Sum[j] <- ifelse(is.null(x), NA, total.mis.prob(x$mean, x$variance$sigmasq, x$pro))
    }

    modes.WSum <- list()
    for(j in 1:length(par.WSum)) {
        x <- par.WSum[[j]]
        if(is.null(x)) {
            modes.WSum[[j]] <- NA
        } else {
            modes.WSum[[j]] <- modes(x$mean, x$variance$sigmasq, x$pro)
        }
    }
    n.modes.WSum <- sapply(modes.WSum, function(x) ifelse(all(is.na(x)), NA, x$nmode))



    mis.WSum <- rep(NA, length(par.WSum))
    for(j in 1:length(par.WSum)) {
        x <- par.WSum[[j]]
        mis.WSum[j] <- ifelse(is.null(x), NA, total.mis.prob(x$mean, x$variance$sigmasq, x$pro))
    }


    modes.HC <- list()
    for(j in 1:length(par.HC)) {
        x <- par.HC[[j]]
        if(is.null(x)) {
            modes.HC[[j]] <- NA
        } else {
            modes.HC[[j]] <- modes(x$mean, x$variance$sigmasq, x$pro)
        }
    }
    n.modes.HC <- sapply(modes.HC, function(x) ifelse(all(is.na(x)), NA, x$nmode))



    mis.HC <- rep(NA, length(par.HC))
    for(j in 1:length(par.HC)) {
        x <- par.HC[[j]]
        mis.HC[j] <- ifelse(is.null(x), NA, total.mis.prob(x$mean, x$variance$sigmasq, x$pro))
    }

    load(paste("scan_results/chr_", i, "_region_test.rdata", sep = ""))

    ##summary.U.Sum <- list()
    summary.U.Sum[[i]] <- data.frame(chr = i, ch.regions.Sum, pval = pval.Sum, nmode = n.modes.Sum, mis = mis.Sum)

    ##summary.U.WSum <- list()
    summary.U.WSum[[i]] <- data.frame(chr = i, ch.regions.WSum, pval = pval.WSum, nmode = n.modes.WSum, mis = mis.WSum)

    ##summary.U.HC <- list()
    summary.U.HC[[i]] <- data.frame(chr = i, ch.regions.HC, pval = pval.HC, nmode = n.modes.HC, mis = mis.HC)
}

save(summary.U.Sum, summary.U.WSum, summary.U.HC, file = "summary/summary_U.rdata")




#########################
## Summarize D results ##
#########################

summary.D.Sum <- list()
summary.D.WSum <- list()
summary.D.HC <- list()


for(i in 1:22) {
    print(i)
    load(paste("test_D_results/chr_", i, "_test_results_D.rdata", sep = ""))
    pval.Sum <- lapply(test.result.Sum, function(y) sapply(y, function(x) ifelse(is.null(x), NA, x$p.value)))
    pval.WSum <- lapply(test.result.WSum, function(y) sapply(y, function(x) ifelse(is.null(x), NA, x$p.value)))
    pval.HC <- lapply(test.result.HC, function(y) sapply(y, function(x) ifelse(is.null(x), NA, x$p.value)))

    par.Sum <- lapply(test.result.Sum, function(y) lapply(y, function(x) x$clst.para))
    par.WSum <- lapply(test.result.WSum, function(y) lapply(y,  function(x) x$clst.para))
    par.HC <- lapply(test.result.HC, function(y) lapply(y, function(x) x$clst.para))

    modes.Sum <- list()
    for(k in 1:3) {
        temp <- list()
        for(j in 1:length(par.Sum[[k]])) {
            x <- par.Sum[[k]][[j]]
            if(is.null(x)) {
                temp[[j]] <- NA
            } else {
                temp[[j]] <- modes(x$mean, x$variance$sigmasq, x$pro)
            }
        }
        modes.Sum[[k]] <- temp
    }
    n.modes.Sum <- lapply(modes.Sum, function(y) sapply(y, function(x) ifelse(all(is.na(x)), NA, x$nmode)))



    mis.Sum <- list()
    for(k in 1:3) {
        temp <- rep(NA, length(par.Sum[[k]]))
        for(j in 1:length(par.Sum[[k]])) {
            x <- par.Sum[[k]][[j]]
            temp[j] <- ifelse(is.null(x), NA, total.mis.prob(x$mean, x$variance$sigmasq, x$pro))
        }
        mis.Sum[[k]] <- temp
    }


    modes.WSum <- list()
    for(k in 1:3) {
        temp <- list()
        for(j in 1:length(par.WSum[[k]])) {
            x <- par.WSum[[k]][[j]]
            if(is.null(x)) {
                temp[[j]] <- NA
            } else {
                temp[[j]] <- modes(x$mean, x$variance$sigmasq, x$pro)
            }
        }
        modes.WSum[[k]] <- temp
    }
    n.modes.WSum <- lapply(modes.WSum, function(y) sapply(y, function(x) ifelse(all(is.na(x)), NA, x$nmode)))



    mis.WSum <- list()
    for(k in 1:3) {
        temp <- rep(NA, length(par.WSum[[k]]))
        for(j in 1:length(par.WSum[[k]])) {
            x <- par.WSum[[k]][[j]]
            temp[j] <- ifelse(is.null(x), NA, total.mis.prob(x$mean, x$variance$sigmasq, x$pro))
        }
        mis.WSum[[k]] <- temp
    }



    modes.HC <- list()
    for(k in 1:3) {
        temp <- list()
        for(j in 1:length(par.HC[[k]])) {
            x <- par.HC[[k]][[j]]
            if(is.null(x)) {
                temp[[j]] <- NA
            } else {
                temp[[j]] <- modes(x$mean, x$variance$sigmasq, x$pro)
            }
        }
        modes.HC[[k]] <- temp
    }
    n.modes.HC <- lapply(modes.HC, function(y) sapply(y, function(x) ifelse(all(is.na(x)), NA, x$nmode)))



    mis.HC <- list()
    for(k in 1:3) {
        temp <- rep(NA, length(par.HC[[k]]))
        for(j in 1:length(par.HC[[k]])) {
            x <- par.HC[[k]][[j]]
            temp[j] <- ifelse(is.null(x), NA, total.mis.prob(x$mean, x$variance$sigmasq, x$pro))
        }
        mis.HC[[k]] <- temp
    }

    load(paste("scan_D_results/chr_", i, "_chp_test.rdata", sep = ""))

    hh <- c(10, 20, 30)

    chps.Sum <- lapply(1:3, function(x) data.frame(chr = i, chp = chp.max.Sum[[x]], hh = hh[x], pval = pval.Sum[[x]], nmode = n.modes.Sum[[x]], mis = mis.Sum[[x]]))
    summary.D.Sum[[i]] <- do.call(rbind, chps.Sum)

    chps.WSum <- lapply(1:3, function(x) data.frame(chr = i, chp = chp.max.WSum[[x]], hh = hh[x], pval = pval.WSum[[x]], nmode = n.modes.WSum[[x]], mis = mis.WSum[[x]]))
    summary.D.WSum[[i]] <- do.call(rbind, chps.WSum)

    chps.HC <- lapply(1:3, function(x) data.frame(chr = i, chp = chp.max.HC[[x]], hh = hh[x], pval = pval.HC[[x]], nmode = n.modes.HC[[x]], mis = mis.HC[[x]]))
    summary.D.HC[[i]] <- do.call(rbind, chps.HC)
}

save(summary.D.Sum, summary.D.WSum, summary.D.HC, file = "summary/summary_D.rdata")


##elig.Sum <- list()
##for(k in 1:3) {
##    elig.Sum[[k]] <- which(n.modes.Sum[[k]] > 1)
##}
##k=1
##    sapply(test.result.Sum[[k]][elig.Sum[[k]]], function(x) ifelse(is.null(x), NA, x$p.value))



