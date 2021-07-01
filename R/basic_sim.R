library(cmdstanr)
library(posterior)
library(matrixStats)
library(HDInterval)
library(mvnfast)
library(parallel)
library(data.table)
library(optiRum) # for CJ.dt
library(qs) # for qsave and qread
library(optparse)

option_list = list(
  make_option(c("-c", "--cores"), type="integer", default=20,
              help="number of cores to use", metavar="character"),
  make_option(c("-n", "--nsim"), type="integer", default=2000,
              help="number of simulations to run under each configuration", metavar="character")
);

opt <- parse_args(OptionParser(option_list = option_list))
num_cores <- opt$cores
num_sims  <- opt$nsim

logreg <- cmdstan_model("Stan/logistic_regression.stan")

# Possible Designs
Xtrt <- cbind(1, rbind(0, rbind(cbind(diag(1, 2), 0), 1)))
Xeff <- cbind(1, c(-1, 1, -1, 1), c(-1, -1, 1, 1), c(1, -1, -1, 1))
Xef2 <- cbind(1, 0.5*c(-1, 1, -1, 1), 0.5*c(-1, -1, 1, 1), 0.25*c(1, -1, -1, 1))
colnames(Xtrt) <- colnames(Xeff) <- colnames(Xef2) <- c("1", "x1", "x2", "x1x2")

# Contrasts of interest in terms of cell means
C <- matrix(
  c(-0.5,  0.5, -0.5, 0.5,
    -0.5, -0.5,  0.5, 0.5,
       1,   -1,   -1,   1,
      -1,    1,    0,   0,
      -1,    0,    1,   0,
       0,    0,   -1,   1,
       0,   -1,    0,   1), 7, 4, byrow = T)
rownames(C) <- c(
  "Main effect x1", "Main effect x2",
  "Interaction x1:x2", "x1 alone", "x2 alone",
  "x1 when given with x2 (vs x2 alone)",
  "x2 when given with x1 (vs x1 alone)")

# Priors considered
mu0non <- rep(0, 4)
mu0 <- c(qlogis(1/5), rep(0, 3))
Sigma0non <- diag(rep(10, 4)^2)
Sigma0 <- diag(c(1.8, 1, 1, 1)^2)

# Simulation configurations
design_configs <- data.table::data.table(
  X = list(Xef2, Xef2[, 1:3], Xef2, Xef2[, 1:3]),
  C = list(C, C, C, C),
  mu0 = list(mu0, mu0[1:3], mu0non, mu0non[1:3]),
  Sigma = list(Sigma0, Sigma0[1:3, 1:3], Sigma0, Sigma0[1:3, 1:3])
)
sim_configs <- data.table::CJ(
  sims = num_sims,
  n_each = 25,
  p = list(c(.1,.1,.1,.1), c(.1,.4,.1,.4), c(.1, .4, .2, .4), c(.1, .4, .2, .6)),
  sorted = F
)
configs <- CJ.dt(
  sim_configs,
  design_configs
)

# Run simulation under given configuration
run_sim <- function(config) {
  sims <- config[, sims][1]
  n_each <- config[, n_each][[1]]
  mu0 <- config[, mu0][[1]]
  Sigma0 <- config[, Sigma][[1]]
  p <- config[, p][[1]]
  X <- config[, X][[1]]
  C <- config[, C][[1]]

  K <- ncol(X)
  N <- nrow(X)
  n <- rep(n_each, N)
  dat <- list(N = N, K = K, n = n, X = X, mu = mu0, Sigma = Sigma0)

  res <- mclapply(1:sims, function(x) {
    dat$y <- rbinom(N, n, p)
    fit <- logreg$sample(
      data = dat,
      iter_warmup = 500,
      iter_sampling = 5000,
      refresh = 0,
      show_messages = F,
      chains = 1)
    beta <- as_draws_matrix(fit$draws("beta"))
    mu <- beta %*% t(X)
    p <- plogis(mu)
    contr <- mu %*% t(C)
    list(mu_beta = colMeans(beta),
         med_beta = colMedians(beta),
         sig_beta = var(beta),
         mu_mu = colMeans(mu),
         sig_mu = var(mu),
         mu_p = colMeans(p),
         med_p = colMedians(p),
         mu_contr = colMeans(contr),
         sig_contr = var(contr),
         ci_contr = HDInterval::hdi(contr),
         gt0_contr = colMeans(contr > 0))
  }, mc.cores = num_cores)
  names(res) <- paste0("sim", 1:sims)
  return(res)
}

# Generate sim results
res <- lapply(1:nrow(configs), function(x) run_sim(configs[x]))

# Save results
qsave(configs, file = "~/out_files/mental_health_sims/basic_configs.qs")
qsave(res, file = "~/out_files/mental_health_sims/basic_sim_results.qs")
