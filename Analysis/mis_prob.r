total.mis.prob <- function(means, variances, probs) {
  n.class <- length(means)
  roots <- NULL
  for (i in 1:(n.class - 1)) {
    for (j in (i + 1):n.class) {
      roots <- c(roots, meet.points(means[c(i, j)], variances[c(i, j)], probs[c(i, j)]))
    }
  }
  if (is.null(roots)) {
      prob <- 1 - max(probs)
      return(prob)
  }
  roots <- sort(roots)
  n.roots <- length(roots)
  mid.points <- colMeans(rbind(c(min(roots) - 1, roots), c(roots, max(roots) + 1)))
  highest <- rep(NA, n.roots + 1)
  log.prob.C <- rep(NA, n.roots + 1)
  for (i in 1:(n.roots + 1)) {
    hh <- which.max(log(probs) + dnorm(mid.points[i], mean = means, sd = sqrt(variances), log = TRUE))
    highest[i] <- hh
    if (i == 1) {
      log.prob.C[i] <- log(probs[hh]) + pnorm(roots[1], mean = means[hh], sd = sqrt(variances[hh]), log.p = TRUE)
    }
    if (i == (n.roots + 1)) {
      log.prob.C[i] <- log(probs[hh]) + pnorm(roots[n.roots], mean = means[hh], sd = sqrt(variances[hh]), lower.tail = FALSE, log.p = TRUE)
    }
    if (i < (n.roots + 1) && i > 1) {
      ph <- log(probs[hh]) + pnorm(roots[i], mean = means[hh], sd = sqrt(variances[hh]), log.p = TRUE)
      pl <- log(probs[hh]) + pnorm(roots[i - 1], mean = means[hh], sd = sqrt(variances[hh]), log.p = TRUE)
      log.prob.C[i] <- ph + log(1 - exp(pl - ph))
    }
  }
  prob <-  1 - exp(max(log.prob.C)) * sum(exp(log.prob.C - max(log.prob.C)))
  return(prob)
}

meet.points <- function(mu, sigma2, p) {
  A <- sigma2[2] - sigma2[1]
  B <- -2 * (mu[1] * sigma2[2] - mu[2] * sigma2[1])
  C <- mu[1]^2 * sigma2[2] - mu[2]^2 * sigma2[1] - 2 * sigma2[1] * sigma2[2] * log(p[1]/p[2] * sqrt(sigma2[2]/sigma2[1]))
  delta <- B^2- 4 * A * C
  if (delta < 0) return(NULL)
  roots <- (-B + c(-1, 1) * sqrt(delta))/(2 * A)
  return(roots)
}

mis.prob <- function(mu, sigma2, p) {
  roots <- meet.points(mu, sigma2, p)
  if (is.null(roots)) return(min(p))
  roots.mean <- mean(roots)
  higher <- which.max(log(p) + dnorm(roots.mean, mean = mu, sd = sqrt(sigma2), log = TRUE))
  lower <- 3 - higher
  log.prob1 <- log(p[higher]) + pnorm(roots[1], mean = mu[higher], sd = sqrt(sigma2[higher]), log.p = TRUE)
  log.prob2 <- log(p[lower]) + pnorm(roots[1], mean = mu[lower], sd = sqrt(sigma2[lower]), log.p = TRUE)
  log.prob3 <- log(p[lower]) + pnorm(roots[2], mean = mu[lower], sd = sqrt(sigma2[lower]), log.p = TRUE)
  log.prob4 <- log(p[higher]) + pnorm(roots[2], mean = mu[higher], sd = sqrt(sigma2[higher]), lower.tail = FALSE  , log.p = TRUE)
  log.probs <- c(log.prob1, log.prob2, log.prob3, log.prob4)
  prob <- exp(max(log.probs)) * sum(exp(log.probs - max(log.probs)) * c(1, -1, 1, 1))
  return(prob)
}

modes <- function(means, vars, pro, L = 1000) {
    from <- min(means) - 0.1
    to <- max(means) + 0.1
    n <- length(means)
    den <- function(x) {
        d <- 0
        for(i in 1:n) {
            d <- d + pro[i] * dnorm(x, means[i], sqrt(vars[i]))
        }
        return(d)
    }
    xs <- seq(from, to, length.out = L)
    ys <- den(xs)
    cond1 <- ys[2:(L - 1)] > ys[1:(L - 2)]
    cond2 <- ys[2:(L - 1)] > ys[3:L]
    cond <- cond1 & cond2
    return(list(den = den, nmode = sum(cond), modes.x = xs[2:(L - 1)][cond], modes.y = ys[2:(L - 1)][cond]))
}
