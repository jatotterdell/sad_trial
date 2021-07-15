# run_basic_sequential_sim.R

source("R/basic_sequential_sim.R")

# Data Generation setup 1
pt <- ptrunc(0:56, spec = "logis", a = -0.5, b = 56, 17.5, 3)
p0 <- diff(c(0, pt))
m0 <- sum(p0*0:56)
s0 <- sqrt(sum(p0*(0:56 - m0)^2))
or <- c(1/16, 1/8, 1/4, 1/2)
p1 <- pomodm(p = p0, odds.ratio = or[1])
m1 <- sum(p1*0:56)
s1 <- sqrt(sum(p1*(0:56 - m1)^2))
p2 <- pomodm(p = p0, odds.ratio = or[2])
m2 <- sum(p2*0:56)
p3 <- pomodm(p = p0, odds.ratio = or[3])
m3 <- sum(p3*0:56)
p4 <- pomodm(p = p0, odds.ratio = or[4])
m4 <- sum(p4*0:56)


p_alt00 <- rbind(p0, p0, p0)
p_alt11 <- rbind(p0, p1, p1)
p_alt22 <- rbind(p0, p2, p2)
p_alt33 <- rbind(p0, p3, p3)
p_alt44 <- rbind(p0, p4, p4)
p_alt01 <- rbind(p0, p0, p1)
p_alt21 <- rbind(p0, p2, p1)
p_alt31 <- rbind(p0, p3, p1)
p_alt41 <- rbind(p0, p4, p1)
p_alt02 <- rbind(p0, p0, p2)
p_alt12 <- rbind(p0, p1, p2)
p_alt32 <- rbind(p0, p3, p2)
p_alt42 <- rbind(p0, p4, p2)
p_alt03 <- rbind(p0, p0, p3)
p_alt13 <- rbind(p0, p1, p3)
p_alt23 <- rbind(p0, p2, p3)
p_alt43 <- rbind(p0, p4, p3)
p_alt04 <- rbind(p0, p0, p4)
p_alt14 <- rbind(p0, p1, p4)
p_alt24 <- rbind(p0, p2, p4)
p_alt34 <- rbind(p0, p3, p4)
p_alt10 <- rbind(p0, p1, p0)
p_alt20 <- rbind(p0, p2, p0)
p_alt30 <- rbind(p0, p3, p0)
p_alt40 <- rbind(p0, p4, p0)

looks <- list(seq(30, 150, 30), seq(30, 240, 30))
p <- list(
  p_alt01, p_alt11, p_alt21, p_alt31,
  p_alt02, p_alt12, p_alt22, p_alt32,
  p_alt03, p_alt13, p_alt23, p_alt33,
  p_alt00, p_alt10, p_alt20, p_alt30)
epsilon <- list(0.95, 0.975)
mu0 <- list(c(18, 0, 0))
Sigma0 <- list(diag(c(10^2, 3.9^2, 3.9^2)))
drop_rule <- list("one")
MID <- list(0)
NonInfMarg <- list(3.3)
scenarios_sup <- scenario_list(n_list(p, epsilon, mu0, Sigma0, drop_rule, MID, NonInfMarg, looks))

epsilon <- list(0.8, 0.85)
MID <- list(5)
scenarios_mid <- scenario_list(n_list(p, epsilon, mu0, Sigma0, drop_rule, MID, NonInfMarg, looks))

scenarios <- c(scenarios_sup, scenarios_mid)

scenarios_dt <- rbindlist(lapply(scenarios, function(x)
  as.data.table(lapply(x, function(a) ifelse(length(a)>1, list(a), I(a))))))[, scenario := 1:.N]
scenarios_dt[, effect1 := rep(
  c("null", "very large", "large", "moderate",
    "null", "very large", "large", "moderate",
    "null", "very large", "large", "moderate",
    "null", "very large", "large", "moderate"), 8)]
scenarios_dt[, effect2 := rep(c(
  "very large", "very large", "very large", "very large",
  "large", "large", "large", "large",
  "moderate", "moderate", "moderate", "moderate",
  "null", "null", "null", "null"), 8)]

res_sup <- rbindlist(pbapply::pblapply(scenarios_sup, function(x) sim_replicate(1000, trial_par = x, cores = 15)), idcol = "scenario")
res_mid <- rbindlist(pbapply::pblapply(scenarios_mid, function(x) sim_replicate(1000, trial_par = x, cores = 15)), idcol = "scenario")
res_mid[, scenario := scenario + max(res_sup$scenario)]

res <- rbind(res_sup, res_mid)
res_sum <- res[, lapply(.SD, mean), keyby = .(scenario, analysis)]

saveRDS(scenarios_dt, file = "~/out_files/mental_health_sims/basic_sequential_scenarios2.rds")
saveRDS(res, file = "~/out_files/mental_health_sims/basic_sequential_sims2.rds")
saveRDS(res_sum, file = "~/out_files/mental_health_sims/basic_sequential_summary2.rds")

rm(res, res_sum)
