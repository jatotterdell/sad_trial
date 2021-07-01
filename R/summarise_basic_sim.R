library(matrixStats)
library(parallel)
library(data.table)
library(optiRum) # for CJ.dt
library(qs) # for qsave and qread

configs <- qread(file = "~/out_files/mental_health_sims/basic_configs.qs")
res     <- qread(file = "~/out_files/mental_health_sims/basic_sim_results.qs")

X <- configs[1, X][[1]]
C <- configs[1, C][[1]]

ORtrue <- lapply(configs[, p], \(x) exp(MASS::ginv(cbind(1, rbind(0, diag(1, 3)))) %*% qlogis(x)))

# Summarising functions
exp_relative_ci <- function(res) {
  apply(apply(exp(simplify2array(lapply(res, `[[`, "ci_contr"))), 2:3, \(x) x[2]/x[1]), 1, mean)
}
exp_ci_width <- function(res) {
  apply(apply(simplify2array(lapply(res, `[[`, "ci_contr")), 2:3, diff), 1, mean)
}
exp_decide_contr_gt0 <- function(res, thres = 0.975) {
  rowMeans(simplify2array(lapply(res, `[[`, "gt0_contr")) > thres)
}
exp_par <- function(res, par = "mu_beta") {
  rowMeans(simplify2array(lapply(res, \(x) `[[`(x, par))))
}
exp_par_var <- function(res, par = "sig_beta") {
  rowMeans(simplify2array(lapply(res, \(x) sqrt(diag(`[[`(x, par))))))
}

p <- do.call(rbind, configs$p)

# Power
round(do.call(rbind, lapply(res, exp_decide_contr_gt0, thres = 0.95)), 3)

rbindlist(lapply(res, function(x) setDT(as.list(exp_par_var(x)))[]), fill = T)
e_par <- rbindlist(lapply(res, function(x) setDT(as.list(exp_par(x)))[]), fill = T)
e_par <- rbindlist(lapply(res, function(x) setDT(as.list(exp_par(x, "med_beta")))[]), fill = T)

cbind(
  Est = as.numeric(e_par[2, ]),
  True = drop(MASS::ginv(X) %*% qlogis(p[2, ])))

round(do.call(rbind, lapply(res, exp_decide_contr_gt0, thres = 0.95)), 3)
do.call(rbind, lapply(res, exp_relative_ci))
do.call(rbind, lapply(res, exp_ci_width))
do.call(rbind, lapply(res, exp_par_var))
do.call(rbind, lapply(res, exp_par))
round(do.call(rbind, lapply(res, exp_par, "mu_contr")), 2)
