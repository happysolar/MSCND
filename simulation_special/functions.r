higher.criticism <- function(p)
{
    n <- length(p)
    pp <- sort(p)
    HC <- ((1:n)/n - pp)/sqrt(pp * (1 - pp)) * sqrt(n)
    return(max(abs(HC)))
}

Max <- function(p)
{
    return(min(p))
}

a.rOP <- function(p, u = 1/4, l = 1/2)
{
    n <- length(p)
    pp <- sort(p)
    li <- floor(n * (1 - l))
    ui <- ceiling(n * (1 - u))
    slopes <- (1 - pp[li:ui])/(n + 1 - (li:ui))
    rs <- ceiling(n + 1 - 1/slopes)
    r <- max(c(rs, 1))
    rOP <- pp[r]
    p.rOP <- pbeta(rOP, r, n + 1 - r)
    return(p.rOP)
}

a.rOP.2 <- function(p, u = 1/4, l = 1/2)
{
    n <- length(p)
    pp <- sort(p)
    li <- floor(n * (1 - l))
    ui <- ceiling(n * (1 - u))
    slopes <- (1 - pp[li:ui])/(n + 1 - (li:ui))
    rs <- ceiling(n + 1 - 1/slopes)
    r <- max(c(rs, 1))
    rOP <- pp[r]
    return(abs(rOP-r/n)/sqrt(rOP*(1-rOP)/n))
}

a.rOP.3 <- function(p, u = 1/4, l = 1/2)
{
    n <- length(p)
    pp <- sort(p)
    li <- floor(n * (1 - l))
    ui <- ceiling(n * (1 - u))
    slopes <- (1 - pp[li:ui])/(n + 1 - (li:ui))
    rs <- ceiling(n + 1 - 1/slopes)
    r <- max(c(rs, 1))
    rOP <- pp[r]
    return(abs(rOP-r/n)/sqrt(rOP*(1-rOP)/n))
}

min.error <- function(label, stat, higher = FALSE) {
    if(higher == FALSE) stat <- -stat
    or <- order(stat)
    n <- length(label)
    err <- 1
    for(i in 0:n) {
        if(i == 0) {
            err <- min(err, mean(label == 0))
        } else if(i == n) {
            err <- min(err, mean(label == 1))
        } else {
            err.1 <- mean(label[or[1:i]] == 1)
            err.2 <- mean(label[or[(i+1):n]] == 0)
            err <- min(err, err.1 + err.2)
        }
    }
    return(err)
}


min.error.2 <- function(label, stat, higher = FALSE) {
    if(higher == FALSE) stat <- -stat
    or <- order(stat)
    n <- length(label)
    num.h0=sum(label==0)
    num.h1=sum(label==1)
    label.order <- label[or]
    err <- 1
    for(i in 2:n) {
      err1 <- sum(label.order[i:n]==0)/num.h0
      err2 <- sum(label.order[1:(i-1)]==1)/num.h1
      err=min(err,err1+err2)
    }
    return(err)
}

