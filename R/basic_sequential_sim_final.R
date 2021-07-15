library(simstudy)
library(data.table)
library(varapproxr) # remotes::install_github("jatotterdell/varapproxr")
library(rctrandr)   # remotes::install_github("jatotterdell/rctrandr")
library(sn)
library(truncdist)
library(Hmisc)
library(parallel)


#' n_list
#'
#' Create a named list from objects (from loo::n_list)
#'
#' @param ... Objects to include in the list
#' @return A named list
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


#' scenario_list
#'
#' Combine a list of scenario parameters into a combined list
#'
#' @param ... Scenario parameters
#' @return A list
scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}


#' sim_accrual
#'
#' Simulate participant accrual into a trial
#'
#' @param n The number of participants to enrol
#' @param accrual_rate The rate of accrual assuming exponential gap times
#' @tte The time between randomisation and the primary outcome
#' @return A data.table of accrual data
sim_accrual <- function(n, accrual_rate = 1.5, tte = 8) {
  ds <- data.table(id = 1:n, gap = rexp(n, 1.5))
  ds[, t0 := cumsum(gap)]
  ds[, tt := t0 + tte]
  return(ds)
}


#' sim_outcome
#'
#' Simulate outcome data
#'
#' @param n The number of participant outcomes to simulate
#' @param ntrt The number of treatments
#' @param p The probability of respones y
#' @return a data.table of outcomes
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


#' trial_as_dat
#'
#' Convert a list of trial results into a data.table
#'
#' @param res A named list of trial output objects
#' @return The list converted to data.table with one row per analysis
trial_as_dt <- function(res) {
  Reduce(function(x, y) merge(x, y, on = .(analysis, variable), all = TRUE),
         lapply(res, function(d) as.data.table(d, keep.rownames = "analysis")))
}


#' Simulate a trial for Phase II SAD.
#'
#' @param trial_par A list providing the trial parameters.
#' - looks: the sample size at which analyses occur
#' - epsilon: a vector of decision thresholds,
#'    1 - superiority,
#'    2 - futility,
#'    3 - inferiority,
#'    4 - highly effective
#' - p: a matrix of dimension 3 by 57 giving probability of y = 0:56 in each treatment arm
#' - mu0: prior mean (intercept, effect1, effect2)
#' - Sigma0: prior covariance
#' - Delta: reference effects,
#'    1 - futility (half the large MD),
#'    2 - treatment relative efficacy (large MD),
#'    3 - highly effective (large MD)
#' - stop_rule: 1 stop if
#' - drop_rule: 1 drop control if one treatment superior, 2 drop only if both superior
#' @return A data.table consisting of trial results at each interim analysis
#' @examples
#' sim_trial()
sim_trial <- function(
  trial_par = list(
    looks = c(30, 60, 90, 120, 150),
    epsilon = c(0.95, 0.05, 0.05, 0.85),
    p = {
      pt <- ptrunc(0:55, spec = "logis", a = 0, b = 56, 20, 4)
      p <- diff(c(0, pt, 1))
      rbind(p, p, p)
    },
    mu0 = c(18, 0, 0),
    Sigma0 = diag(c(10^2, 6^2, 6^2)),
    Delta = c(2.5, 5, 5, 0),
    stop_rule = 1,
    drop_rule = 1
  )
) {

  looks <- trial_par$looks
  epsilon <- trial_par$epsilon
  p <- trial_par$p
  mu0 <- trial_par$mu0
  Sigma0 <- trial_par$Sigma0
  Delta <- trial_par$Delta
  stop_rule <- trial_par$stop_rule
  drop_rule <- trial_par$drop_rule

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
  Cdes <- rbind(cbind(-1, diag(1, P - 1)), c(0, 1, -1), c(0, -1, 1))
  colnames(Cdes) <- colnames(Xdes)
  rownames(Cdes) <- c("t1-t0", "t2-t0", "t1-t2", "t2-t1")

  t_seq <- sapply(looks, function(a) dat[, tt][a])
  n_enr <- sapply(looks, function(a) nrow(dat[t0 <= tt[a]]))
  n_new <- diff(c(0, n_enr))
  id_obs <- cbind(c(1, looks[-K] + 1), looks)
  id_enr <- cbind(c(1, n_enr[-length(n_enr)] + 1), n_enr)

  # Storage
  trtlabs <- paste0("trt", 0:(P - 1))
  ss <- expand.grid(0:2, 0:2)
  ssalt <- expand.grid(1:nrow(Cdes), 1:nrow(Cdes))
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
  c_sigma <- matrix(0, K, nrow(Cdes)*(nrow(Cdes)+1)/2, dimnames = list(analysis = 1:K, parameter = paste0("c_sigma_", ctrlabs)))

  # Specific probabilities to monitor
  p_mon <- matrix(0, K, nrow(Cdes) + 6, dimnames = list(
    analysis = 1:K,
    probability = c(paste0("p_", c(rownames(Cdes))),
                    "p_t1-t0-d1", "p_t2-t0-d1",
                    "p_t1-t2-d2", "p_t2-t1-d2",
                    "p_t1-t0-d3", "p_t2-t0-d3")))
  i_drp <- matrix(0, K, 3, dimnames = list(analysis = 1:K, treatment = c("drp0", "drp1", "drp2")))
  i_fut <- matrix(0, K, 2, dimnames = list(analysis = 1:K, treatment = c("fut1", "fut2")))
  i_eff <- matrix(0, K, 2, dimnames = list(analysis = 1:K, treatment = c("eff1", "eff2")))
  i_inf <- matrix(0, K, 2, dimnames = list(analysis = 1:K, treatment = c("inf1", "inf2")))
  i_heff <- matrix(0, K, 2, dimnames = list(analysis = 1:K, treatment = c("heff1", "heff2")))
  i_stp <- matrix(0, K, 1, dimnames = list(analysis = 1:K, stop = "stp"))
  status <- matrix(0, K, 1, dimnames = list(analysis = 1:K, status = 'status'))

  final <- FALSE
  stopped <- FALSE
  for(l in 1:length(looks)) {

    if(l == length(looks)) final <- TRUE

    if(stopped) {

      final <- TRUE
      findat <- dat[!is.na(trt)]
      n_trt_enr[l, ] <- findat[, .N, keyby = trt][["N"]]
      n_trt_obs[l, ] <- findat[, .N, keyby = trt][["N"]]
      moddat <- trtdat[findat, on = .(id, trt)]
      y_trt_obs[l, ] <- moddat[, .(m = mean(y)), keyby = trt][["m"]]

    } else {

      trtenr <- rctrandr:::mass_weighted_urn_rand(palloc, n_new[l], alpha = 3)$trt
      dat[between(id, id_enr[l, 1], id_enr[l, 2]), trt := trtenr]
      n_trt_enr[l, ] <- dat[between(id, 1, id_enr[l, 2]), .N, keyby = trt][["N"]]
      n_trt_obs[l, ] <- dat[between(id, 1, id_obs[l, 2]), .N, keyby = trt][["N"]]
      moddat <- trtdat[dat[between(id, 1, id_obs[l, 2])], on = .(id, trt)]
      y_trt_obs[l, ] <- moddat[between(id, 1, id_obs[l, 2]), .(m = mean(y)), keyby = trt][["m"]]

    }

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
    p_mon[l, ] <- c(pnorm(-Delta[4], ctr_mu[1:2], sqrt(diag(ctr_sig)[1:2])), # TRT better than PBO by Delta[4] (effectiveness)
                    pnorm(0, ctr_mu[3:4], sqrt(diag(ctr_sig)[3:4])),         # TRT better than other TRT (relative effectiveness)
                    pnorm(-Delta[1], ctr_mu[1:2], sqrt(diag(ctr_sig)[1:2])), # TRT better than PBO by Delta[1] (futility)
                    pnorm(-Delta[2], ctr_mu[3:4], sqrt(diag(ctr_sig)[3:4])), # TRT better than other TRT by Delta[2] (relative efficacy)
                    pnorm(-Delta[3], ctr_mu[1:2], sqrt(diag(ctr_sig)[1:2]))) # TRT better than PBO by Delta[3] (highly effective)
    ### Actions
    # Drop control if either treatment superior (rule 1) or only if both superior (rule 2)
    # Drop a treatment if futile (unlikely to improve more than half MID)
    # Drop a treatment if other treatment superior
    if(l == 1) {
      i_eff[l, ] <- p_mon[l, 1:2] > epsilon[1]
      i_fut[l, ] <- p_mon[l, 5:6] < epsilon[2]
      i_inf[l, ] <- p_mon[l, 3:4] < epsilon[3]
      i_heff[l, ] <- p_mon[l, 9:10] > epsilon[4]
      if(drop_rule == 1) {
        i_drp[l, ] <- c(
          any(p_mon[l, 1:2] > epsilon[1]), # drop PBO if either treatment superior
          (p_mon[l, 5:6] < epsilon[2]) | (p_mon[l, 3:4] < epsilon[3]) # treatment futile or inferior to other treatment
        )
      } else if(drop_rule == 2) {
        i_drp[l, ] <- c(
          all(i_eff[l, ] == 1), # drop PBO if both treatments superior
          (p_mon[l, 5:6] < epsilon[2]) | (p_mon[l, 3:4] < epsilon[3]) # treatment futile or inferior to other treatment
        )
      }

    } else {
      i_eff[l, ] <- (p_mon[l, 1:2] > epsilon[1]) | i_eff[l-1, ]
      i_fut[l, ] <- (p_mon[l, 5:6] < epsilon[2]) | i_fut[l-1, ]
      i_inf[l, ] <- p_mon[l, 3:4] < epsilon[3] | i_inf[l-1, ]
      i_heff[l, ] <- (p_mon[l, 9:10] > epsilon[4]) | i_heff[l-1, ]

      if(drop_rule == 1) {
        i_drp[l, ] <- c(
          any(p_mon[l, 1:2] > epsilon[1]), # either treatment superior
          p_mon[l, 5:6] < epsilon[2] | (p_mon[l, 3:4] < epsilon[3]) # treatment futile or inferior to other treatment
        ) | i_drp[l - 1, ] # already been dropped then stays dropped
      } else if(drop_rule == 2) {
        i_drp[l, ] <- c(
          all(i_eff[l, ] == 1), # both treatments superior
          p_mon[l, 5:6] < epsilon[2] | (p_mon[l, 3:4] < epsilon[3]) # treatment futile or inferior to other treatment
        ) | i_drp[l - 1, ] # already been dropped then stays dropped
      }
    }

    # Stop rule 1:
    # - both treatments futile
    # - both treatments highly effective
    # - one treatment effective and better than the other
    # Stop rule 2:
    # - both treatments futile
    # - both treatments effective
    # - one treatment effective and better than the other
    if(stop_rule == 1) {
      i_stp[l, ] <- all(i_fut[l, ] == 1) | all(i_heff[l, ] == 1) | (i_eff[l, 1] & i_inf[l, 2]) | (i_eff[l, 2] & i_inf[l, 1])
    } else if(stop_rule == 2) {
      i_stp[l, ] <- all(i_fut[l, ] == 1) | all(i_eff[l, ] == 1) | (i_eff[l, 1] & i_inf[l, 2]) | (i_eff[l, 2] & i_inf[l, 1])
    }
    stopped <- i_stp[l, ]

    status[l, ] <- ifelse(
      final == TRUE, "final analysis",
      ifelse(all(i_heff[l, ] == 1), "both highly effective",
      ifelse(all(i_fut[l, ] == 1), "both futile",
      ifelse((i_eff[l, 1] & i_inf[l, 2]) | (i_eff[l, 2] & i_inf[l, 1]), "one effective and better", "continue")
    )))

    if(final) break
    # Stop if only one treatment active OR if both treatments highly effective (>6 reduction HAM-A)
    palloc  <- palloc * (1 - i_drp[l, ]) / sum(palloc * (1 - i_drp[l, ]))


  }

  out <- n_list(
    n_trt_enr[1:l, , drop = F],
    n_trt_obs[1:l, , drop = F],
    y_trt_obs[1:l, , drop = F],
    b_mu[1:l, , drop = F],
    b_sigma[1:l, , drop = F],
    m_mu[1:l, , drop = F],
    m_sigma[1:l, , drop = F],
    c_mu[1:l, , drop = F],
    c_sigma[1:l, , drop = F],
    p_mon[1:l, , drop = F],
    i_eff[1:l, , drop = F],
    i_fut[1:l, , drop = F],
    i_inf[1:l, , drop = F],
    i_heff[1:l, , drop = F],
    i_drp[1:l, , drop = F],
    i_stp[1:l, , drop = F],
    status[1:l, , drop = F]
  )
  return(trial_as_dt(out))
}


#' sim_replicate
#'
#' Replicate a `sim_trial` multiple times for a given scenario
#'
#' @param nsim The number of simulations to repeat
#' @param cores The number of cores to use
#' @param ... The sim_trial scenario
#' @return A data.table of combined simulated trials
#' sim_replicate_mc(10)
sim_replicate_mc <- function(nsim, cores = 14, ...) {
  rbindlist(
    parallel::mclapply(
      X = 1 : nsim,
      FUN = function(z) sim_trial(...),
      mc.cores = cores),
    idcol = "trial"
  )
}

sim_replicate <- function(nsim, ...) {
  rbindlist(
    lapply(
      X = 1 : nsim,
      FUN = function(z) sim_trial(...)),
    idcol = "trial"
  )
}
