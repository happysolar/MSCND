library("mclust")
library("aod")
#i <- 22
#load(file = paste("scan_results/chr_", i, "_region_test.rdata", sep = ""))

#N <- dim(U.test.Sum)[1]
#W <- dim(U.test.Sum)[2]

#num.clusters.U.Sum <- rep(NA, W)

#for (j in 1:W) {
#  signal <- U.test.Sum[, j]
#  cluster <- Mclust(signal, modelNames = "V", G = 1:5)
#}

test.binary <- function(y, x, z, ...) {
  N <- length(x)
  cluster.results <- Mclust(x, ...)
  if (cluster.results$G == 1) return(NULL)
  probs <- cluster.results$z
  clst.para <- cluster.results$parameters
  dummy <- probs[, -1, drop = FALSE]
  data <- data.frame(y = y, x = dummy, z = z)
  glm.results <- glm(y ~ ., family = "binomial", data = data)
  estimate <- coef(glm.results)[2:(cluster.results$G)]
  sd <- sqrt(diag(vcov(glm.results)))[2:(cluster.results$G)]
  p.value <- wald.test(b = coef(glm.results), Sigma = vcov(glm.results), Terms = 2:(cluster.results$G))$result$chi2[3]
  return(list(estimate = estimate, sd = sd, p.value = p.value, clst.para = clst.para))
}
