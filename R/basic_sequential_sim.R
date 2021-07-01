library(simstudy)
library(data.table)
library(varapproxr) # remotes::install_github("jatotterdell/varapproxr")
library(rctrandr)   # remotes::install_github("jatotterdell/rctrandr")
library(sn)
library(truncdist)
library(Hmisc)

n_list <- function (...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names)
    FALSE
  else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  }
  else {
    names(out)[!has_name] <- nms[!has_name]
  }
  return(out)
}

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

sim_accrual <- function(n, accrual_rate = 1.5, tte = 8) {
  ds <- data.table(id = 1:n, gap = rexp(n, 1.5))
  ds[, t0 := cumsum(gap)]
  ds[, tt := t0 + tte]
  return(ds)
}

sim_outcome <- function(n, ntrt, p) {
  ds <- CJ(
    id = 1:n,
    trt = 1:ntrt,
    sorted = FALSE
  )
  setkey(ds, id)
  for(i in 1:ntrt) {
    r <- rlogis(n)
    cp <- qlogis(cumsum(c(0, p[i, ])))
    ds[trt == i, y := findInterval(r, cp) - 1]
  }
  return(ds)
}

sim_data <- function(n, ntrt, p) {
  ds1 <- data.table(id = 1:n, gap = rexp(n, 1.5))
  ds1[, t0 := cumsum(gap)]
  ds1[, tt := t0 + 8]

  ds2 <- CJ(id = 1:n, trt = 1:ntrt, sorted = TRUE)
  setkey(ds2, id)
  for(i in 1:ntrt) {
    r <- rlogis(n)
    cp <- qlogis(cumsum(c(0, p[i, ])))
    ds2[trt == i, y := findInterval(r, cp) - 1]
  }
  return(list(ds1, ds2))
}

trial_as_dt <- function(res) {
  Reduce(function(x, y) merge(x, y, on = .(analysis, variable), all = TRUE),
         lapply(res, function(d) as.data.table(d, keep.rownames = "analysis")))
}

sim_trial <- function(
  trial_par = list(
    looks = c(30, 60, 90, 120, 150),
    epsilon = 0.95,
    p = {
      pt <- ptrunc(0:55, spec = "logis", a = 0, b = 56, 20, 4)
      p <- diff(c(0, pt, 1))
      rbind(p, p, p)
    },
    mu0 = c(20, 0, 0),
    Sigma0 = diag(5, 3)^2,
    MID = 5,
    NonInfMarg = 3.3,
    drop_rule = "mid"
  )
) {

  looks <- trial_par$looks
  epsilon <- trial_par$epsilon
  p <- trial_par$p
  mu0 <- trial_par$mu0
  Sigma0 <- trial_par$Sigma0
  drop_rule <- trial_par$drop_rule
  MID <- trial_par$MID
  NonInfMarg <- trial_par$NonInfMarg

  K <- length(looks)
  N <- max(looks)
  palloc = rep(1 / 3, 3)
  P <- length(palloc)

  dat <- sim_accrual(N)
  dat[, trt := NA_integer_]
  trtdat <- sim_outcome(N, P, p)

  # Basic design matrix
  Xdes <- cbind(1, rbind(0, diag(1, P - 1)))
  colnames(Xdes) <- paste0("b_mu_", 1:P - 1)
  rownames(Xdes) <- paste0("mu", 1:P)

  # Contrast matrix
  Cdes <- rbind(cbind(-1, diag(1, P - 1)), c(0, 1, -1))
  rownames(Cdes) <- c("trt1-trt0", "trt2-trt0", "trt1-trt2")

  t_seq <- sapply(looks, function(a) dat[, tt][a])
  n_enr <- sapply(looks, function(a) nrow(dat[t0 <= tt[a]]))
  n_new <- diff(c(0, n_enr))
  id_obs <- cbind(c(1, looks[-K] + 1), looks)
  id_enr <- cbind(c(1, n_enr[-length(n_enr)] + 1), n_enr)

  # Storage
  trtlabs <- paste0("trt", 0:(P - 1))
  ss <- expand.grid(0:2, 0:2)
  ssalt <- expand.grid(1:3, 1:3)
  varlabs <- apply(ss[ss[, 2] <= ss[, 1], ], 1, paste0, collapse = "")
  ctrlabs <- apply(ssalt[ssalt[, 2] <= ssalt[, 1], ], 1, paste0, collapse = "")
  n_trt_enr <- matrix(0, K, P, dimnames = list(analysis = 1:K, treatment = paste0("n_enr_", trtlabs)))
  n_trt_obs <- matrix(0, K, P, dimnames = list(analysis = 1:K, treatment = paste0("n_obs_", trtlabs)))
  y_trt_obs <- matrix(0, K, P, dimnames = list(analysis = 1:K, treatment = paste0("y_obs_", trtlabs)))
  b_mu  <- matrix(0, K, P, dimnames = list(analysis = 1:K, parameter = colnames(Xdes)))
  b_sigma <- matrix(0, K, P*(P+1)/2, dimnames = list(analysis = 1:K, parameter = paste0("b_sigma_", varlabs)))
  m_mu <- b_mu
  colnames(m_mu) <- gsub("b", "m", colnames(b_mu))
  m_sigma <- b_sigma
  colnames(m_sigma) <- gsub("b", "m", colnames(b_sigma))
  c_mu <-  matrix(0, K, nrow(Cdes), dimnames = list(analysis = 1:K, parameter = paste0("c_mu_", 1:nrow(Cdes))))
  c_sigma <- matrix(0, K, 3*(3+1)/2, dimnames = list(analysis = 1:K, parameter = paste0("c_sigma_", ctrlabs)))
  p_ctr <- matrix(0, K, 4, dimnames = list(analysis = 1:K, probability = paste0("p_", c(rownames(Cdes), "joint"))))
  drp0 <- matrix(0, K, 1, dimnames = list(analysis = 1:K, drop = "drp0"))

  for(l in 1:length(looks)) {
    # trtenr <- sample.int(P, n_new[l], replace = TRUE, prob = palloc)
    trtenr <- rctrandr:::mass_weighted_urn_rand(palloc, n_new[l], alpha = 3)$trt
    dat[between(id, id_enr[l, 1], id_enr[l, 2]), trt := trtenr]
    n_trt_enr[l, ] <- dat[between(id, 1, id_enr[l, 2]), .N, keyby = trt][["N"]]
    n_trt_obs[l, ] <- dat[between(id, 1, id_obs[l, 2]), .N, keyby = trt][["N"]]
    moddat <- trtdat[dat[between(id, 1, id_obs[l, 2])], on = .(id, trt)]
    y_trt_obs[l, ] <- moddat[between(id, 1, id_obs[l, 2]), .(m = mean(y)), keyby = trt][["m"]]

    X <- model.matrix( ~ factor(trt), data = moddat)[,]
    y <- moddat$y
    fit <- vb_lm(X, y, mu0, Sigma0, 1e-2, 1e-2)
    mu <- Xdes %*% fit$mu
    Sig <- Xdes %*% fit$Sigma %*% t(Xdes)
    ctr_mu <- c(Cdes %*% mu)
    ctr_sig <- Cdes %*% Sig %*% t(Cdes)
    b_mu[l, ] <- c(fit$mu)
    b_sigma[l, ] <- fit$Sigma[lower.tri(fit$Sigma, diag = TRUE)]
    m_mu[l, ] <- mu
    m_sigma[l, ] <- Sig[lower.tri(Sig, diag = T)]
    c_mu[l, ] <- ctr_mu
    c_sigma[l, ] <- ctr_sig[lower.tri(ctr_sig, diag = T)]
    p_ctr[l, ] <- c(pnorm(c(-MID, -MID, NonInfMarg), ctr_mu, sqrt(diag(ctr_sig))),
                    mvtnorm::pmvnorm(-Inf, 0, ctr_mu[1:2], sigma = ctr_sig[1:2,1:2]))
    ### Actions

    # Drop options:
    # - as soon as one superior by MID
    # - as soon as both superior by MID
    # - as soon as jointly superior by MID
    if(drop_rule == "one") {
      if(any(p_ctr[l, 1:2] > epsilon)) {
        drp0[l, ] <- 1
        palloc[1] <- 0
        palloc[-1] <- palloc[-1] / sum(palloc[-1])
      } else if(l > 1) {
        drp0[l, ] <- drp0[l - 1, ]
      }
    } else if (drop_rule == "two") {
        if(all(p_ctr[l, 1:2] > epsilon)) {
          drp0[l, ] <- 1
          palloc[1] <- 0
          palloc[-1] <- palloc[-1] / sum(palloc[-1])
        } else if(l > 1) {
          drp0[l, ] <- drp0[l - 1, ]
        }
    } else if (drop_rule == "joint") {
        if(p_ctr[l, 4] > epsilon) {
          drp0[l, ] <- 1
          palloc[1] <- 0
          palloc[-1] <- palloc[-1] / sum(palloc[-1])
        } else if(l > 1) {
          drp0[l, ] <- drp0[l - 1, ]
        }
    }
  }
  out <- n_list(
    n_trt_enr,
    n_trt_obs,
    y_trt_obs,
    b_mu,
    b_sigma,
    m_mu,
    m_sigma,
    c_mu,
    c_sigma,
    p_ctr,
    drp0
  )
  return(trial_as_dt(out))
}

sim_replicate <- function(nsim, cores = 14, ...) {
  rbindlist(
    parallel::mclapply(
      X = 1 : nsim,
      FUN = function(x) sim_trial(...),
      mc.cores = cores),
    idcol = "trial"
  )
}
